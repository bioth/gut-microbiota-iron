#!/bin/bash

read_length_distribution() {
    # Check if the path variable is set
    if [ -z "$path" ]; then
        echo "Error: 'path' variable is not set."
        return 1
    fi
    for file in "$path"/*.fastq; do
        # Extracting file name without extension
        filename=$(basename -- "$file")
        filename_no_ext="${filename%.*}"
        echo "$filename_no_ext" >> "$path/reads_counts_distribution.txt"
        awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' "$file" >> "$path/reads_counts_distribution.txt"
    done
}