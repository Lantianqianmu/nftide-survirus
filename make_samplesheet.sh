#!/bin/bash
# Produces samplesheet.csv from a directory containing sample subfolders.
# Each subfolder may contain paired _1.fq.gz / _2.fq.gz files.
# Output: samplesheet.csv with columns: folder_name,read1,read2

OUTPUT_FILE="samplesheet.csv"

# Get the directory where the script is run (should contain the sample folders)
BASE_DIR="$(pwd)"

# Write header
echo "sample,fastq_1,fastq_2" > "$OUTPUT_FILE"

# Iterate over each subdirectory in the base directory
for sample_dir in "$BASE_DIR"/*/; do
    # Skip if not a directory
    [ -d "$sample_dir" ] || continue
    
    # Get the folder name
    folder_name="$(basename "$sample_dir")"
    
    # Find all paired _1.fq.gz / _2.fq.gz files
    # Collect _1.fq.gz files
    mapfile -t r1_files < <(find "$sample_dir" -maxdepth 1 -name '*_1.fq.gz' -type f | sort)
    
    for r1 in "${r1_files[@]}"; do
        # Derive the expected _2.fq.gz path by replacing _1.fq.gz with _2.fq.gz
        r2="${r1/_1.fq.gz/_2.fq.gz}"
        
        # Check if the paired file exists
        if [ -f "$r2" ]; then
            # Get absolute paths
            r1_abs="$(realpath "$r1")"
            r2_abs="$(realpath "$r2")"
            
            # Write to samplesheet
            echo "${folder_name},${r1_abs},${r2_abs}" >> "$OUTPUT_FILE"
        fi
    done
done

echo "samplesheet.csv generated with $(tail -n +2 "$OUTPUT_FILE" | wc -l) entries."