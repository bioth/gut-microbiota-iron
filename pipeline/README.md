# 1- Remove the primers at the end of the reads (trimming of non-biological bases):

## First uncompress the fastq files with: 
bash uncompress.sh

## Secondly use cutadapt to remove the primers: 
bash cutadapt.sh

# 2- Filtering of the reads using DADA2:

## dada2Tests.r

