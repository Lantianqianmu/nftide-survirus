# nftide-survirus Pipeline #

Nextflow pipeline for detecting HBV integration from paired-end HBV-probe enriched libraries using sample-specific HBV genomes and SurVirus.

## System requirements ##

The pipeline runs on Linux with Bash, Conda or Mamba, and GCC/G++. SurVirus uses Python 2.7; Nextflow uses Java. The custom filtering process compiles a small alignment library with `gcc` and `g++`, so both commands must be available inside the task environment or through the system PATH.

Plan storage and memory for a complete host+HBV FASTA and BWA index **for each sample**, plus intermediate FASTQs and BAMs. The pipeline reuses the host-only BWA index but builds each merged index separately. A fixed minimum memory requirement has not been benchmarked for this workflow.

The current `nextflow.config` allocates 16 CPUs per process by default and allows up to six concurrent tasks per process, including six SurVirus tasks. Cutadapt and custom filtering use one CPU each. Reduce concurrency on smaller servers. Published results use copy mode, and `cleanup = false` retains Nextflow work files for resuming and diagnosis.

## Software dependencies ##

Dependencies | Environment / version used in local verification
------------- | -------------
Nextflow | `placseq`; 25.10.2
Java / OpenJDK | `placseq`; Java 17
Python | `survirus`; 2.7.15
NumPy | `survirus`; 1.16.5
PyFaidx | `survirus`; 0.7.0
PySam | `survirus`; 0.20.0
Cutadapt | `survirus`; 1.18
BWA | `survirus`; 0.7.18
Samtools | `survirus`; 1.18, using HTSlib 1.17
sdust | `survirus`; required for low-complexity masking
SurVirus | Compiled installation at `~/SurVirus`
GCC / G++ | System or task PATH; C++11 support required

## Installation ##

(1) Use **two separate Conda environments**:

Environment | Purpose | How it is used
------------- | ------------- | -------------
`placseq` | Launch Nextflow and provide Java | Activate it in the terminal before starting the workflow.
`survirus` | Run Python-2 SurVirus and its bioinformatics dependencies | Nextflow activates it automatically for processing tasks.

SurVirus's `surveyor.py` contains Python-2 syntax and requires its legacy Python dependencies, whereas Nextflow needs Java. Separate environments isolate those dependency stacks and prevent updates to the workflow launcher from disturbing SurVirus. Two environments are the supported setup for this pipeline, not an intrinsic Nextflow requirement. These names describe the existing local setup; another environment containing compatible Nextflow and Java can launch the pipeline.

The existing environment paths are:

```text
/home/xrz/miniforge3/envs/placseq
/home/xrz/miniforge3/envs/survirus
```

Nextflow uses `conda.enabled = true` in `nextflow.config`. The Cutadapt, host FASTA indexing, SurVirus, and custom filtering processes declare `conda params.survirus_conda`. The launch shell can therefore stay in `placseq` throughout execution. Do not switch the launch shell to `survirus` to run these tasks manually.

(2) Check the existing installations:

```bash
conda activate placseq
java -version
nextflow -version

conda run -n survirus python -c 'import numpy, pyfaidx, pysam, cutadapt'
conda run -n survirus samtools --version
conda run -n survirus cutadapt --version
conda run -n survirus python /home/xrz/SurVirus/surveyor.py --help
conda run -n survirus bash -c 'command -v bwa; command -v sdust; command -v gcc; command -v g++'
```

On another server, provision the two environments separately with the dependencies above, compile SurVirus, and set `--survirus_conda` and `--survirus_dir` to their actual paths. Keep the SurVirus source `libs/` directory as well as its compiled executables: the custom filter builds its alignment helper from those sources. Record working environments with `conda list -n placseq --explicit` and `conda list -n survirus --explicit` for reproducibility.

(3) Prepare the host FASTA and its BWA index. For prefix `hg38.fa`, the host directory must contain:

```text
hg38.fa
hg38.fa.amb
hg38.fa.ann
hg38.fa.bwt
hg38.fa.pac
hg38.fa.sa
```

For a new host reference, run `bwa index /path/to/hg38.fa` using the SurVirus environment. The workflow generates the host FASTA `.fai` once in its task directory. Do not supply only the BWA sidecar files: the host FASTA itself is required.

## Overview ##

This pipeline merges FASTQ files sharing a sample name, trims adapters, and extracts that sample's HBV sequence from a CSV. It then concatenates the host FASTA with the sample-specific HBV FASTA and builds the HBV-only and merged BWA and FASTA indexes.

SurVirus runs on the trimmed paired FASTQs in `--fq` mode. It identifies candidate host–virus junctions, constructs breakpoint sequences, remaps supporting reads, masks low-complexity sequence with sdust, and generates its original filtered results. An independent `CUSTOM_FILTER_SURVIRUS` process applies configurable criteria to **all candidates in `results.remapped.txt`**, preserving accepted calls, rejected calls, and the reasons for each decision.

## Usage ##

(1) Prepare `samplesheet.csv`. The CSV file __must__ contain these three columns:

`sample`: Sample/library name, used in output paths. Multiple rows with the same name are merged before trimming. Use simple names containing letters, numbers, underscores or hyphens; avoid spaces and shell-special characters.  
`fastq_1`: Path to the gzipped read-1 FASTQ.  
`fastq_2`: Path to the matching gzipped read-2 FASTQ.

```csv
sample,fastq_1,fastq_2
N10,/data/reads/N10_lane1_R1.fq.gz,/data/reads/N10_lane1_R2.fq.gz
N10,/data/reads/N10_lane2_R1.fq.gz,/data/reads/N10_lane2_R2.fq.gz
```

(2) Prepare `meta_all_assembled_fa.csv`. The CSV file __must__ contain:

`sample`: Exactly matching sample name from the samplesheet; provide one reference row per sample.  
`sequence`: The HBV nucleotide sequence itself, on a single CSV line, without a FASTA header. This is not a path to a FASTA file.

```csv
sample,sequence
N10,ACGTACGTACGTACGT
```

The short sequence above illustrates the format only; use the actual assembled HBV genome for analysis. Each emitted viral FASTA uses the contig name `HBV`; that name must not already identify a host contig. HBV coordinates consequently refer to the supplied sample-specific sequence, which may differ between samples.

(3) Activate the Nextflow environment and execute:

```bash
conda activate placseq
cd /data/xrz/capint/nftide-survirus

nextflow run main.nf \
  -output-dir /data/xrz/capint/output_survirus \
  --input_csv samplesheet.csv \
  --HBVfa_csv meta_all_assembled_fa.csv \
  --host_bwa_dir /data/xrz/ref/hg38/hg38_bwa \
  --host_bwa_prefix hg38.fa \
  --survirus_conda /home/xrz/miniforge3/envs/survirus \
  --survirus_dir /home/xrz/SurVirus \
  --filter_mask_mode flag \
  -with-report survirus_report.html \
  -with-timeline survirus_timeline.html \
  -bg -resume
```

`-resume` reuses eligible completed tasks; retain both `work/` and `.nextflow/`. `-bg` runs in the background. Inspect `.nextflow.log` for progress and failures. Nextflow options use one hyphen; pipeline parameters use two.

### General parameters ###

Parameter | Default in the current pipeline | Description
------------- | ------------- | -------------
`-output-dir` | `/data/xrz/capint/output_survirus` | Published output directory; a Nextflow option.
`--input_csv` | `/data/xrz/capint/nftide-survirus/samplesheet.csv` | FASTQ samplesheet.
`--HBVfa_csv` | `/data/xrz/capint/nftide-survirus/meta_all_assembled_fa.csv` | Sample-specific HBV sequences.
`--host_bwa_dir` | `/data/xrz/ref/hg38/hg38_bwa` | Host FASTA and BWA sidecars directory.
`--host_bwa_prefix` | `hg38.fa` | Host FASTA filename and index prefix, including extension.
`--survirus_conda` | `$HOME/miniforge3/envs/survirus` | Existing Conda environment directory, not merely its name.
`--survirus_dir` | `$HOME/SurVirus` | Compiled SurVirus installation and source directory; capitalization matters on Linux.
`--survirus_wgs` | `false` | Passes SurVirus's `--wgs` when true. The pipeline always uses FASTQ mode; the installed SurVirus uses its FASTQ insert-size estimator, so this switch does not select its BAM sampling path.
`--filter_workdir` | unset | Existing SurVirus result directory; when set, executes custom filtering only.
`--filter_sample` | `sample` | Output sample label in filtering-only mode; set explicitly, e.g. `N10`.

Resource settings are in `nextflow.config`, not exposed as custom `--cpus` parameters. Cutadapt uses `-j 1` because its Python-2-compatible release cannot run parallel workers. The adapter sequence is `CTGTCTCTTATACACATCT` for both reads; minimum length is `2:2`, minimum overlap is 1, and pair filtering is `any`. These trimming settings are currently fixed in `main.nf`.

### Custom filtering parameters ###

Parameter | Default passed by Nextflow | Description
------------- | ------------- | -------------
`--filter_min_host_pbs` | `0.8` | Minimum HOST_PBS; equality passes. NaN/infinite values fail.
`--filter_min_pairs` | `2` | Supporting pairs sufficient without split-read evidence.
`--filter_split_min_pairs` | `1` | Alternative minimum supporting pairs when enough split reads exist.
`--filter_min_split_reads` | `1` | Split reads required for the alternative support rule.
`--filter_min_seq_length` | `30` | Minimum length, in bases, of each host and virus breakpoint sequence.
`--filter_mask_mode` | `flag` | `flag`: annotate masking without rejecting for it. `reject`: reject if either masked fraction exceeds the threshold.
`--filter_max_masked_fraction` | `0.8` | Masked-fraction threshold; exactly 80% passes and >80% is flagged/rejected.
`--filter_pairing_rescue` | `true` | Enable pairing and recovery of rejected partners.
`--filter_pair_min_distance` | `-50` | Minimum signed host breakpoint separation for pairing, in bases.
`--filter_pair_max_distance` | `1000` | Maximum signed host breakpoint separation for pairing, in bases.
`--filter_deduplicate` | `true` | Remove later calls matching an earlier retained call on both sides.
`--filter_dedup_fraction` | `0.8` | Minimum locally aligned fraction of the shorter sequence, on each side, to classify a duplicate.

These are the defaults in **`main.nf`**. The standalone Python script defaults to `--mask-mode reject`; Nextflow explicitly passes its own default, `flag`.

(1) **Initial quality and support filtering.** A candidate must have finite `HOST_PBS >= 0.8`, both breakpoint sequences at least 30 bases long, and satisfy:

```text
supporting_pairs >= 2
OR
(supporting_pairs >= 1 AND split_reads >= 1)
```

This is the intended “two pairs or 1 + 1” rule. A split read can come from a supporting pair, so the two counts are not necessarily independent DNA molecules. `HOST_PBS` in the remapped file is SurVirus's normalized alignment score assembled from junction-supporting alignments; it should not be interpreted as a mapping probability or simply host percent identity. Empty sequences fail the length test. **Coverage is reported but never used for rejection or rescue by the custom filter.**

(2) **Low-complexity masking.** The script reads the host and virus sdust BED files and divides the union of masked intervals by the corresponding sequence length. Intervals use zero-based, half-open BED coordinates (`end-start`); overlaps are counted once. Both `flag` and `reject` modes calculate masking. Missing BEDs are errors; present but empty BEDs mean no masked bases.

With the pipeline default `flag`, repetitive candidates can pass the quality filters and receive annotations. `HOST_HAS_MASKED_BASES` and `VIRUS_HAS_MASKED_BASES` indicate any masked bases, while `REPEAT_FLAGS` identifies sides above the configured fraction. These annotations describe **sdust low-complexity sequence**, not overlap with an external genomic repeat annotation such as RepeatMasker.

(3) **Pairing rescue.** Pairing requires the same host chromosome and virus contig, opposite host strands, opposite virus strands, and a signed host distance within the inclusive configured interval. The distance is the reverse-strand host position minus the forward-strand host position; the position is `start` for a minus-strand breakpoint and `end` for a plus-strand breakpoint.

The filter pairs initially passing calls first. Each remaining unpaired passing call can rescue one nearest eligible rejected call. Ties follow candidate-file order. Rescue can override initial failures, including score, support, length, and masking; those failures remain recorded. Rescued calls do not initiate rescue chains. After deduplication, a rescued call is rejected if its original qualifying sponsor was removed. Use `--filter_pairing_rescue false` if every retained call must independently pass all initial criteria.

(4) **Deduplication.** The filter compares candidates in input order against earlier retained calls using striped Smith–Waterman local alignment. It uses match score 1, mismatch penalty 4, gap-open penalty 6, and gap-extension penalty 1, with reverse complementation when the corresponding breakpoint strands differ. Both host and virus sides must satisfy:

```text
(aligned_query_end - aligned_query_begin + 1) / shorter_sequence_length >= 0.8
```

This is an aligned-span threshold, not 80% sequence identity. Matching uses sequence similarity rather than a chromosome-distance window, so different mapped locations may be collapsed. The first qualifying call is kept; the script does not rank duplicates by support. Set `--filter_deduplicate false` to inspect all passing/rescued candidates. The Nextflow process currently compiles the alignment helper even when deduplication is disabled.

### Filtering existing results ###

Change thresholds without rebuilding indexes or rerunning SurVirus:

```bash
conda activate placseq
cd /data/xrz/capint/nftide-survirus

nextflow run main.nf \
  --filter_workdir /path/to/N10_survirus \
  --filter_sample N10 \
  --filter_mask_mode reject \
  --filter_min_pairs 2 \
  --filter_split_min_pairs 1 \
  --filter_min_split_reads 1 \
  --filter_min_seq_length 30 \
  -output-dir /path/to/refiltered_results \
  -resume
```

The directory must contain `results.remapped.txt`, `host_bp_seqs.fa`, `virus_bp_seqs.fa`, `host_bp_seqs.masked.bed`, and `virus_bp_seqs.masked.bed`. Candidate IDs must match the FASTA record IDs. This mode publishes only the custom filtering directory; it does not republish original FASTQs or the full SurVirus directory. Use separate output directories when comparing parameter sets.

## Expected output ##

Go to `-output-dir`. Each sample has the following structure:

```text
<sample>/
├── fastqs/
│   ├── <sample>_cutadapt_R1.fq.gz
│   ├── <sample>_cutadapt_R2.fq.gz
│   └── <sample>_cutadapt.log
└── survirus/
    ├── <sample>_survirus/
    │   ├── results.txt
    │   ├── results.remapped.txt
    │   ├── results.t1.txt
    │   ├── results.remapped.t1.txt
    │   ├── results.discarded.txt
    │   ├── results.alternative.txt
    │   ├── host_bp_seqs.fa
    │   ├── virus_bp_seqs.fa
    │   ├── host_bp_seqs.masked.bed
    │   ├── virus_bp_seqs.masked.bed
    │   ├── config.txt
    │   ├── contig_map
    │   ├── log.txt
    │   ├── bam_0/
    │   └── readsx/
    └── <sample>_custom_filter/
        ├── accepted.tsv
        ├── rejected.tsv
        ├── all_candidates.tsv
        └── summary.json
```

The entire `<sample>_survirus` directory is published, so additional SurVirus intermediates such as sampled FASTQs, SAMs, and breakpoint BAMs may also appear. `bam_0/` contains retained reads and intermediate alignments for the merged paired FASTQ input. `readsx/` contains candidate-associated BAMs where generated.

### Original SurVirus files ###

File | Description
------------- | -------------
`results.txt` | Candidate calls before breakpoint-consensus remapping; internal whitespace-delimited format without a header.
`results.remapped.txt` | Candidates with support/score information updated by remapping. Input to the custom filter; not a final accepted-call list.
`results.t1.txt` | Stock SurVirus filtering of `results.txt`.
`results.remapped.t1.txt` | Stock SurVirus filtering of `results.remapped.txt`. Uses the original SurVirus criteria, independent of custom parameters.
`results.discarded.txt` | Rejections from stock filtering of `results.txt`. It is not a rejection audit of the custom filter or of all remapped candidates; duplicates removed later are not comprehensively represented here.
`results.alternative.txt` | Alternative host mappings produced by SurVirus; not an additional custom-filter stage.
`host_bp_seqs.fa`, `virus_bp_seqs.fa` | Candidate breakpoint sequences, with headers linking to candidate IDs and orientation.
`host_bp_seqs.masked.bed`, `virus_bp_seqs.masked.bed` | sdust intervals relative to the corresponding breakpoint sequences, not host-genome BED coordinates.
`config.txt` | SurVirus runtime settings, including threads, read length and clipping settings.
`contig_map` | Mapping between reference contig names and SurVirus's internal IDs.
`log.txt` | SurVirus remapper log. Full task stdout/stderr are in the Nextflow task directory.

The raw candidate field order is:

```text
id host_breakpoint virus_breakpoint reads good_pairs split_reads score
host_pbs virus_pbs reads_w_dups unique_reads_w_dups host_cov virus_cov
```

These are 13 fields on one line per candidate. The custom filter uses `good_pairs` as `SUPPORTING_PAIRS`. Breakpoint strings retain SurVirus's `contig:strand:start:end` representation. The custom TSV does not convert coordinates to BED or VCF; verify coordinate conventions before exporting to another format. In SurVirus's compact printed breakpoint, minus uses `start` and plus uses `end`.

### Custom filtering files ###

Start with **`<sample>_custom_filter/accepted.tsv`** for calls retained under the custom parameters. This includes both independently passing calls and calls recovered by pairing rescue. Use `all_candidates.tsv` to audit every decision.

File | Description
------------- | -------------
`accepted.tsv` | Final retained candidates; `STATUS` is `accepted` or `rescued`.
`rejected.tsv` | Final rejected candidates, including sequence duplicates, with rejection reasons.
`all_candidates.tsv` | Every input candidate, in original order, with the same columns as the other TSVs.
`summary.json` | Total candidates, initial passes, final retained count, final status counts and parameters used. Zero-count statuses may be absent from the status dictionary.

All TSVs have a header even when there are no records. Their columns are:

Column | Meaning
------------- | -------------
`ID` | Original SurVirus candidate ID; stable across these output tables.
`HOST_BREAKPOINT`, `VIRUS_BREAKPOINT` | Unmodified raw SurVirus breakpoint strings.
`SUPPORTING_PAIRS` | Remapped `good_pairs` count.
`SPLIT_READS` | Remapped split-read count; may overlap supporting-pair evidence.
`HOST_PBS` | Remapped normalized alignment score used by the score filter.
`COVERAGE` | `(host_cov + virus_cov)/2`; reported only. SurVirus derives these coverages relative to its maximum insert-size window, not total sequencing depth.
`HOST_SEQ_LENGTH`, `VIRUS_SEQ_LENGTH` | Breakpoint FASTA sequence lengths in bases.
`HOST_MASKED_FRACTION`, `VIRUS_MASKED_FRACTION` | Fraction of bases covered by the union of sdust intervals. Empty sequences report 0 but fail the length criterion.
`REPEAT_FLAGS` | Semicolon-separated `HOST_LOW_COMPLEXITY` / `VIRUS_LOW_COMPLEXITY` when fractions exceed the threshold; blank otherwise.
`STATUS` | Final `accepted`, `rescued`, or `rejected` status.
`REJECTION_REASON` | Final rejection reasons, separated by semicolons; blank for retained calls. Duplicate rejection is reported as `DUPLICATE_SEQUENCE`.
`INITIAL_FILTER_FAILURES` | Original score, support, length and mask failures, preserved even after rescue or duplicate removal.
`PAIRED_WITH` | Partner ID if both paired calls remain retained; blank otherwise. ID 0 is a valid partner.
`DUPLICATE_OF` | Earlier candidate ID responsible for sequence deduplication; blank if not a duplicate. This records the deduplication decision before final rescue-sponsor cleanup.
`RESCUED_BY` | Original passing sponsor ID for a rescued candidate; preserved even if that candidate is subsequently rejected.
`HOST_HAS_MASKED_BASES`, `VIRUS_HAS_MASKED_BASES` | 1 if any bases on that side are masked, otherwise 0; independent of the >80% threshold.

Rejection codes are `LOW_HOST_PBS`, `INVALID_HOST_PBS`, `INSUFFICIENT_SUPPORT`, `SHORT_HOST_SEQUENCE`, `SHORT_VIRUS_SEQUENCE`, `HOST_LOW_COMPLEXITY`, `VIRUS_LOW_COMPLEXITY`, `DUPLICATE_SEQUENCE`, and `RESCUE_PARTNER_REMOVED`. Multiple initial failures can coexist. In `flag` mode, low-complexity annotations do not become initial rejection reasons. Malformed inputs, missing files or duplicate IDs stop the task rather than silently dropping records.

### Work files and reports ###

Merged pre-trimming FASTQs, sample HBV/host+HBV FASTAs and their indexes, and the compiled custom alignment helper reside in Nextflow's `work/` task directories. They are outside `<sample>_survirus` and are **not included in that published folder**. Published trimmed FASTQs and the full SurVirus intermediate directory can consume substantial additional disk space because output mode is `copy`.

`-with-report` and `-with-timeline` write reports at the specified paths. A failed task's work directory contains `.command.sh`, `.command.out`, `.command.err`, and `.exitcode`; check these and the Cutadapt/SurVirus logs when diagnosing a failure.

## Verification ##

The independent custom filtering process was verified through Nextflow 25.10.2 on N10's 411 candidates with the corrected alignment-length calculation:

Mask mode | Initially passing | Final retained | Final rejected
------------- | ------------- | ------------- | -------------
`reject` | 297 | 267 | 144
`flag` | 341 | 295 | 116

No N10 candidates were rescued in these runs. Many additional retained candidates use the one-pair/one-split support rule; passing these computational filters is not independent experimental validation.

The audit tables are in `tests/n10_validation/fixed_reject/` and `tests/n10_validation/fixed_flag/`. Targeted tests in `tests/test_custom_filter.py` cover thresholds, overlapping masking intervals, empty sequences, pairing rescue and deduplication. The synthetic caller smoke test can be run with:

```bash
bash tests/survirus_smoke_test.sh
```

## Troubleshooting ##

`java: command not found`: Activate `placseq` before launching Nextflow. Calling Nextflow by its absolute path does not automatically add that environment's Java to PATH. Alternatively:

```bash
PATH=/home/xrz/miniforge3/envs/placseq/bin:$PATH \
  /home/xrz/miniforge3/envs/placseq/bin/nextflow run main.nf -resume -bg
```

`Running in parallel is not supported on Python 2`: Keep Cutadapt at `-j 1`; other SurVirus tasks can still use multiple threads.

`libcrypto.so.1.0.0` missing from Samtools: Check the packages in `survirus`; this previously occurred with an incompatible old Samtools build. Resolve compatible packages within that environment, then recheck Python 2 and its modules. Do not substitute a symlink to a different OpenSSL ABI.

Custom filter compilation failure: Verify `gcc`, `g++`, and `${survirus_dir}/libs/ssw.c` / `ssw_cpp.cpp` are available. The SurVirus installation needs its source files as well as its binaries.
