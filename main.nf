#!/home/zeemeeuw/miniconda3/envs/joint/bin/nextflow


// mkdir -p ~/.nextflow/lsp/v25.10
// wget https://github.com/nextflow-io/language-server/releases/download/v25.10.2/language-server-all.jar -O ~/.nextflow/lsp/v25.10/v25.10.2.jar

// nextflow run main.nf -resume -bg


params.input_csv = '/data/xrz/capint/nftide-survirus/samplesheet.csv'


params.HBVfa_csv = "/data/xrz/capint/nftide-survirus/meta_all_assembled_fa.csv"

// prebuilt host bwa index (shared by all samples, no need to build again)
params.host_bwa_dir = "/data/xrz/ref/hg38/hg38_bwa"
params.host_bwa_prefix = "hg38.fa"

// Existing Conda environment and SurVirus installation.  The environment must
// contain bwa, samtools, Python 2, pysam, pyfaidx, NumPy, and sdust.
params.survirus_conda = "${System.getenv('HOME')}/miniforge3/envs/survirus"
params.survirus_dir = "${System.getenv('HOME')}/SurVirus"
params.survirus_wgs = false
params.filter_min_host_pbs = 0.8
params.filter_min_pairs = 2
params.filter_split_min_pairs = 1
params.filter_min_split_reads = 1
params.filter_min_seq_length = 30
params.filter_mask_mode = 'flag'
params.filter_max_masked_fraction = 0.8
params.filter_pairing_rescue = true
params.filter_pair_min_distance = -50
params.filter_pair_max_distance = 1000
params.filter_deduplicate = true
params.filter_dedup_fraction = 0.8



process MERGE_FQ {
    tag "Merging fastq files of ${meta.id}..."
    
    input:
    tuple val(meta), path(r1s) , path(r2s)
    
    output:
    tuple val(meta), path("${meta.id}_merged_R1.fq.gz"), path("${meta.id}_merged_R2.fq.gz"), emit: merged_fq
    
    script:
    """
    cat ${r1s.join(' ')} > ${meta.id}_merged_R1.fq.gz
    cat ${r2s.join(' ')} > ${meta.id}_merged_R2.fq.gz
    """
}


process CUTADAPT {
    tag "cutadapt on ${meta.id}"
    conda params.survirus_conda

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_cutadapt_R1.fq.gz"), path("${meta.id}_cutadapt_R2.fq.gz"), emit: trimmed_reads
    tuple val(meta), path("${meta.id}_cutadapt.log"), emit: cutadapt_log

    script:
    """
    # Cutadapt 1.18 runs under the Python-2 SurVirus environment and does not
    # support parallel execution.
    cutadapt \
        -j 1 -m 2:2 \
        -a "CTGTCTCTTATACACATCT" \
        -A "CTGTCTCTTATACACATCT" \
        --pair-filter=any \
        --overlap 1 \
        -o ${meta.id}_cutadapt_R1.fq.gz \
        -p ${meta.id}_cutadapt_R2.fq.gz \
        ${r1} \
        ${r2} > "${meta.id}_cutadapt.log"

    """
}


process CHECK_SAMPLES {
    tag "Checking that every sample has an HBV reference genome"

    input:
    val missing

    script:
    """
    if [ -n "${missing}" ]; then
        echo "ERROR: no HBV reference genome found for sample(s): ${missing}" >&2
        exit 1
    fi
    echo "all samples have an HBV reference genome"
    """
}


process EXTRACT_HBV_FA {
    tag "Extracting HBV genome of ${meta.id}"

    input:
    tuple val(meta), val(hbv_seq)

    output:
    tuple val(meta), path("${meta.id}_HBV.fa"), emit: hbv_fa

    script:
    """
    printf '>HBV\\n${hbv_seq}\\n' > ${meta.id}_HBV.fa
    """
}


process BUILD_HOST_FAI {
    tag "Indexing host FASTA once"
    conda params.survirus_conda

    input:
    path host_reference

    output:
    path "${params.host_bwa_prefix}.fai", emit: host_fai

    script:
    """
    samtools faidx ${host_reference}
    """
}


process RUN_SURVIRUS {
    tag "SurVirus on ${meta.id}"
    conda params.survirus_conda

    input:
    tuple val(meta), path(r1), path(r2), path(hbv_fa)
    path host_reference
    path host_bwa_index
    path host_fai

    output:
    tuple val(meta), path("${meta.id}_survirus"), emit: survirus_workdir

    script:
    def sample_dir = "${meta.id}_survirus"
    // hbv_fa is staged by Nextflow with its original name (${meta.id}_HBV.fa),
    // so use a distinct local filename for the writable virus reference.
    def virus_fa = "${meta.id}_virus.fa"
    def merged_fa = "${meta.id}_host_plus_HBV.fa"
    def wgs_arg = params.survirus_wgs ? '--wgs' : ''
    """
    set -euo pipefail

    # surveyor.py is invoked through Python, so it only needs to be readable;
    # the source file is not necessarily marked executable.
    test -r "${params.survirus_dir}/surveyor.py"
    command -v bwa
    command -v samtools
    command -v sdust

    # host_reference and its prebuilt BWA sidecar files are staged inputs.
    # host_fai is generated once by BUILD_HOST_FAI and staged as hg38.fa.fai.
    cp ${hbv_fa} ${virus_fa}
    samtools faidx ${virus_fa}
    cat ${host_reference} ${virus_fa} > ${merged_fa}
    samtools faidx ${merged_fa}

    # SurVirus maps reads to both the HBV-only and concatenated references.
    bwa index ${virus_fa}
    bwa index ${merged_fa}

    mkdir -p ${sample_dir}
    python "${params.survirus_dir}/surveyor.py" \\
        "${r1},${r2}" \\
        ${sample_dir} \\
        ${host_reference} \\
        ${virus_fa} \\
        ${merged_fa} \\
        --fq \\
        --threads ${task.cpus} \\
        --bwa bwa \\
        --samtools samtools \\
        --dust sdust \\
        ${wgs_arg}

    """
}

process CUSTOM_FILTER_SURVIRUS {
    tag "Custom filtering ${meta.id}"
    conda params.survirus_conda
    cpus 1

    input:
    tuple val(meta), path(survirus_workdir)
    path filter_script
    path alignment_source

    output:
    tuple val(meta), path("${meta.id}_custom_filter"), emit: filtered

    script:
    """
    set -euo pipefail
    gcc -O2 -fPIC -c "${params.survirus_dir}/libs/ssw.c" -o ssw.o
    g++ -std=c++11 -O2 -fPIC -shared -I"${params.survirus_dir}/libs" \\
        "${alignment_source}" "${params.survirus_dir}/libs/ssw_cpp.cpp" ssw.o -o libsurvirus_alignment.so
    python "${filter_script}" --workdir "${survirus_workdir}" \\
        --outdir "${meta.id}_custom_filter" --ssw-library ./libsurvirus_alignment.so \\
        --min-host-pbs ${params.filter_min_host_pbs} \\
        --min-pairs ${params.filter_min_pairs} \\
        --split-min-pairs ${params.filter_split_min_pairs} \\
        --min-split-reads ${params.filter_min_split_reads} \\
        --min-seq-length ${params.filter_min_seq_length} \\
        --mask-mode '${params.filter_mask_mode}' \\
        --max-masked-fraction ${params.filter_max_masked_fraction} \\
        --pairing-rescue ${params.filter_pairing_rescue} \\
        --pair-min-distance ${params.filter_pair_min_distance} \\
        --pair-max-distance ${params.filter_pair_max_distance} \\
        --deduplicate ${params.filter_deduplicate} \\
        --dedup-fraction ${params.filter_dedup_fraction}
    """
}

// Refilter existing outputs without repeating mapping/indexing:
// nextflow run main.nf --filter_workdir /path/N10_survirus --filter_sample N10
params.filter_workdir = null
params.filter_sample = 'sample'

workflow {

    main:
    if (params.filter_workdir) {
        ch_filter_input = channel.of(tuple([id: params.filter_sample], file(params.filter_workdir, checkIfExists: true)))
        ch_publish_fastqs = channel.empty()
        ch_publish_logs = channel.empty()
        ch_publish_workdirs = channel.empty()
    } else {
    log.info """\
      nftide-hivid (HBV integration breakpoint detection)
      ===================================
      projectDir             :  ${projectDir}
      workingDir             :  ${workflow.outputDir}
      input csv              :  ${params.input_csv}
      HBV meta csv           :  ${params.HBVfa_csv}
      host bwa index         :  ${params.host_bwa_dir}/${params.host_bwa_prefix}
      SurVirus installation  :  ${params.survirus_dir}
    """.stripIndent()

    ch_read_pairs = channel.fromPath(params.input_csv)
    .splitCsv(header:true)
    .map { row -> 
        [
            row.sample,
            row
        ]
    }
    .groupTuple()
    .map { _sample, rows -> 
        rows.withIndex().collect { row, index ->
            row + [rep: index + 1]
        }
    }
    .flatMap { item -> item }
   .map { row -> 

        [
            [
                id: row.sample,
                rep: row.rep,

            ], 
            [
                file(row.fastq_1, checkIfExists: true), 
                file(row.fastq_2, checkIfExists: true)
            ]
        ]
    }
    .map{meta, files -> [meta.subMap(['id']), files]}
    .groupTuple()
    .map { meta, filePairs ->
        [ meta, filePairs.collect { pair -> pair[0] }, filePairs.collect { pair -> pair[1] }]
    }

    MERGE_FQ(ch_read_pairs)
    CUTADAPT(MERGE_FQ.out.merged_fq)

    // HBV reference genome per sample (from meta_all_assembled_fa.csv)
    ch_hbv = channel.fromPath(params.HBVfa_csv)
        .splitCsv(header: true)
        .map { row -> [ [id: row.sample], row.sequence ] }

    // sanity check: every sample must have an HBV reference genome
    // (read both CSVs directly; simple comma-separated format)
    def _samples = new File(params.input_csv).readLines().drop(1).collect { line -> line.split(',')[0].trim() }
    def _hbv_samples = new File(params.HBVfa_csv).readLines().drop(1).collect { line -> line.split(',')[0].trim() }
    def missing_ids = _samples - _hbv_samples
    CHECK_SAMPLES(missing_ids ? missing_ids.join(' ') : '')

    // Extract the HBV reference once for each sample.  RUN_SURVIRUS then
    // builds the sample-specific HBV and host+HBV BWA indices.
    EXTRACT_HBV_FA(ch_hbv)

    // The host BWA index is reused.  Only its FASTA index is generated here,
    // once, because the supplied host directory does not contain hg38.fa.fai.
    ch_host_reference = channel.value(
        file("${params.host_bwa_dir}/${params.host_bwa_prefix}", checkIfExists: true)
    )
    def host_bwa_index_files = ['amb', 'ann', 'bwt', 'pac', 'sa'].collect { suffix ->
        file("${params.host_bwa_dir}/${params.host_bwa_prefix}.${suffix}", checkIfExists: true)
    }
    ch_host_bwa_index = channel.value(host_bwa_index_files)
    BUILD_HOST_FAI(ch_host_reference)
    ch_host_fai = BUILD_HOST_FAI.out.host_fai.collect()

    // Explicitly join by sample ID: metadata maps are not used as join keys.
    ch_trimmed_with_id = CUTADAPT.out.trimmed_reads.map { meta, r1, r2 ->
        tuple(meta.id, meta, r1, r2)
    }
    ch_hbv_with_id = EXTRACT_HBV_FA.out.hbv_fa.map { meta, hbv_fa ->
        tuple(meta.id, hbv_fa)
    }
    ch_survirus_input = ch_trimmed_with_id
        .join(ch_hbv_with_id)
        .map { sample_id, meta, r1, r2, hbv_fa -> tuple(meta, r1, r2, hbv_fa) }

    RUN_SURVIRUS(ch_survirus_input, ch_host_reference, ch_host_bwa_index, ch_host_fai)
    ch_filter_input = RUN_SURVIRUS.out.survirus_workdir
    ch_publish_fastqs = CUTADAPT.out.trimmed_reads
    ch_publish_logs = CUTADAPT.out.cutadapt_log
    ch_publish_workdirs = RUN_SURVIRUS.out.survirus_workdir
    }
    CUSTOM_FILTER_SURVIRUS(ch_filter_input,
        file("${projectDir}/bin/filter_survirus.py", checkIfExists: true),
        file("${projectDir}/bin/survirus_alignment.cpp", checkIfExists: true))




    publish:
    out_trimmed_fastqs = ch_publish_fastqs
    out_cutadapt_logs = ch_publish_logs
    out_survirus_workdirs = ch_publish_workdirs
    out_custom_filter = CUSTOM_FILTER_SURVIRUS.out.filtered


}

output {
    out_custom_filter {
        path { meta, _folder -> "${meta.id}/survirus" }
    }
    out_trimmed_fastqs {
        path { meta, _f1, _f2 -> "${meta.id}/fastqs" }
    }
    out_cutadapt_logs {
        path { meta, _f1 -> "${meta.id}/fastqs" }
    }
    out_survirus_workdirs {
        path { meta, _workdir -> "${meta.id}/survirus" }
    }




}
