# Welcome to the gut microbiota analysis pipeline!
<p align="center">
  <img src="https://github.com/bioth/gut-microbiota-iron/blob/main/pipeline/photos/pipeline.png?raw=true"/>
</p>

Here is a walkthrough of the pipeline as well as the instructions to run it properly. Pipeline is designed to function with MiSeq V3-V4 16S sequencing data (see example below).

# The requirements:
Running this pipeline cannot be done without the required versions and packages:

- A bash command line (not powershell)
- Python (3.0 or later)
- pandas (python library)
- R (4.3.0 or later)
- R tools (build packages)
- BiocManager (version 3.18)
- DADA2 (version 3.18)
- ggplot2 (version 3.4.4 or later)
- phyloseq (version 1.46.0 or later)

# 1- Prepare you data folder and your bash environment:
You data folder should contain a file with your compressed reads, and a folder called metadata which contains the metadata as well as the MiSeqReadSet excel file which contains all the information about the primers and adapters.

<p align="center">
  <img src="https://github.com/bioth/gut-microbiota-iron/blob/main/pipeline/photos/data_folder_format.png?raw=true" height="100" />
</p>

<p align="center">
  <img src="https://github.com/bioth/gut-microbiota-iron/blob/main/pipeline/photos/metadata_folder_format.png?raw=true" height="100" />
</p>

Now open a bash command line, cd into your data folder, cd into the folder containing the reads, and type the following command:
DATA_FOLDER=$(pwd)

# 2- Uncompress the files:
bash uncompress.sh $DATA_FOLDER

# 3- Remove the primers at the end of the reads (trimming of non-biological bases):
## First create files listing the forward and reverse primers:
py create_primers_fasta_file.py

## Secondly use cutadapt to remove the primers (it will automatically select reverse and forward primers accordingly): 
bash cutadapt.sh $DATA_FOLDER

# 4- Filtering of the reads using DADA2:

## dada2Tests.r

