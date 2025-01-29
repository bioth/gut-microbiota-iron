#!/bin/bash

# Changing directory to access the FASTQ files
cd ../../../Microbiota_18/2106590411

# Check if destination directory exists, create it if not
mkdir -p "../trimmed_fastq"

# Delete existing log files if there are any
rm -f ../trimmed_fastq/combined_cutadapt_logs.txt

# Define adapter sequences (modify these accordingly)
FWD_PRIMER="^GGMTTAGATACCCBDGTA" # 
REV_PRIMER="^GGGTYKCGCTCGTTR" # 

for r1_file in *R1.fastq; do

    # Extracting file name without extension
    filename=$(basename -- "$r1_file")
    filename_no_ext="${filename%_R1.fastq}"
    r2_file="${filename_no_ext}_R2.fastq"

    # Run cutadapt in paired-end mode
    cutadapt --cores=0 --discard-untrimmed --no-indels --pair-adapters -g "${FWD_PRIMER};e=0" -G "${REV_PRIMER};e=0" -o "../trimmed_fastq/${filename_no_ext}_R1_trimmed.fastq" -p "../trimmed_fastq/${filename_no_ext}_R2_trimmed.fastq" $r1_file $r2_file > "../trimmed_fastq/${filename}_cutadapt_log.txt"
    
    # Combine cutadapt logs into a single file
    cat "../trimmed_fastq/${filename}_cutadapt_log.txt" >> ../trimmed_fastq/combined_cutadapt_logs.txt
    
    # Remove individual cutadapt log files
    rm "../trimmed_fastq/${filename}_cutadapt_log.txt"
    
done

        #> $LOG_FILE 2>&1  # Redirect output and errors to the log file
        #--report=minimal --discard-untrimmed --cores 0
