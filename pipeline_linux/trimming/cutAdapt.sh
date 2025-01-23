#!/bin/bash

# Changing directory to access the FASTQ files
cd ../../../Microbiota_18/2106590411

# Check if destination directory exists, create it if not
mkdir -p "../trimmed_fastq"

# Supprimer les fichiers de logs existants s'ils existent
rm -f ../trimmed_fastq/combined_r1_cutadapt_logs.txt
rm -f ../trimmed_fastq/combined_r2_cutadapt_logs.txt

# Loop through each R1 FASTQ file
for file in *R1.fastq; do

    # Extracting file name without extension
    filename=$(basename -- "$file")
    filename_no_ext="${filename%.*}"
    
    # Trimming non-biological bases on R1 FASTQ files / -u 18
    cutadapt -g file:../primers/r1_primers.fasta --cores=0 --discard-untrimmed --no-indels -o "../trimmed_fastq/${filename_no_ext}_trimmed.fastq" "$file" > "../trimmed_fastq/${filename}_cutadapt_log.txt"

     # Combine cutadapt logs into a single file
    cat "../trimmed_fastq/${filename}_cutadapt_log.txt" >> ../trimmed_fastq/combined_r1_cutadapt_logs.txt
    
    # Remove individual cutadapt log files
    rm "../trimmed_fastq/${filename}_cutadapt_log.txt"

done

# Loop through each R2 FASTQ file
for file in *R2.fastq; do

    # Extracting file name without extension
    filename=$(basename -- "$file")
    filename_no_ext="${filename%.*}"

    # Trimming non-biological bases on R2 FASTQ files / -u 15 
    cutadapt -g file:../primers/r2_primers.fasta --cores=0 --discard-untrimmed --no-indels -o "../trimmed_fastq/${filename_no_ext}_trimmed.fastq" "$file" > "../trimmed_fastq/${filename}_cutadapt_log.txt"

    # Combine cutadapt logs into a single file
    cat "../trimmed_fastq/${filename}_cutadapt_log.txt" >> ../trimmed_fastq/combined_r2_cutadapt_logs.txt
    
    # Remove individual cutadapt log files
    rm "../trimmed_fastq/${filename}_cutadapt_log.txt"

done
