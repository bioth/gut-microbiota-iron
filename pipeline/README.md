# 1- Remove the primers at the end of the reads (trimming of non-biological bases):

## First uncompress the fastq files with: 
bash uncompress.sh

## Secondly create files listing the forward and reverse primers:
py create_primers_fasta_file.py

## Thirdly use cutadapt to remove the primers (it will automatically select reverse and forward primers accordingly): 
bash cutadapt.sh

# 2- Filtering of the reads using DADA2:

## dada2Tests.r

