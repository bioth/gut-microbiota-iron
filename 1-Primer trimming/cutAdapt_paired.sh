#!/bin/bash

# Changing directory to access the FASTQ files
#cd ../../../Microbiota_18/2106590411
cd ../../../Microbiota_19/1068099565

# Check if destination directory exists, create it if not
mkdir -p "../trimmed_fastq"

# Delete existing log files if there are any
rm -f ../trimmed_fastq/combined_cutadapt_logs.txt

# Define adapter sequences (modify these accordingly)
FWD_PRIMER="^GGMTTAGATACCCBDGTA" # 
REV_PRIMER="^GGGTYKCGCTCGTTR" # 

for r1_file in *R1.fastq.gz; do

    # Extracting file name without extension
    filename=$(basename -- "$r1_file")
    filename_no_ext="${filename%_R1.fastq.gz}"
    r2_file="${filename_no_ext}_R2.fastq.gz"

    # Run cutadapt in paired-end mode / --no-indels / ;e=0
    cutadapt --cores=0 --discard-untrimmed  --pair-adapters -g "${FWD_PRIMER}" -G "${REV_PRIMER}" -o "../trimmed_fastq/${filename_no_ext}_R1_trimmed.fastq.gz" -p "../trimmed_fastq/${filename_no_ext}_R2_trimmed.fastq.gz" $r1_file $r2_file > "../trimmed_fastq/${filename}_cutadapt_log.txt"
    
    # Combine cutadapt logs into a single file
    cat "../trimmed_fastq/${filename}_cutadapt_log.txt" >> ../trimmed_fastq/combined_cutadapt_logs.txt
    
    # Remove individual cutadapt log files
    rm "../trimmed_fastq/${filename}_cutadapt_log.txt"
    
done

        #> $LOG_FILE 2>&1  # Redirect output and errors to the log file
        #--report=minimal --discard-untrimmed --cores 0
