# Welcome to the gut microbiota analysis pipeline!
Here is a walkthrough of the pipeline as well as the instructions to run it properly. Pipeline is designed to function with MiSeq V3-V4 16S sequencing data (see example below).

# 1- Prepare you data folder:
You data folder should contain a file with your compressed reads, and a folder called metadata which contains the metadata as well as the MiSeqReadSet excel file which contains all the information about the primers and adapters.

<p align="center">
  <img src="https://github.com/bioth/gut-microbiota-iron/blob/main/pipeline/photos/folder_format.png?raw=true" height="100" />
</p>


# 2- Remove the primers at the end of the reads (trimming of non-biological bases):

## First uncompress the fastq files with: 
bash uncompress.sh

## Secondly create files listing the forward and reverse primers:
py create_primers_fasta_file.py

## Thirdly use cutadapt to remove the primers (it will automatically select reverse and forward primers accordingly): 
bash cutadapt.sh

# 3- Filtering of the reads using DADA2:

## dada2Tests.r

