#!/bin/bash

# Source directory containing compressed files
source_dir="../../16s_data/2060250822"

# Destination directory for uncompressed files
destination_dir="../../16s_data/uncompressed_files"

# Check if source directory exists and if it contains any .gz files
if [ ! -d "$source_dir" ]; then
    echo "Source directory does not exist: $source_dir"
    exit 1
fi

gz_files=("$source_dir"/*.gz)
if [ ${#gz_files[@]} -eq 0 ]; then
    echo "No .gz files found in source directory: $source_dir"
    exit 1
fi

# Check if destination directory exists, create it if not
mkdir -p "$destination_dir"

# Loop through all compressed files in the source directory
for file in "${gz_files[@]}"; do
    # Extract filename without extension
    filename=$(basename "$file" .gz)
    # Uncompress the file and save to destination directory
    gunzip -c "$file" > "$destination_dir/$filename"
done
