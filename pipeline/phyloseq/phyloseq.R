#loading libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

library(phyloseq)
library(Biostrings)
library(ggplot2)

theme_set(theme_bw())

#set working directory
setwd("D:/CHUM_git/16s_data/asv_table/")
asv_table <- read.csv("seqtab.nochim2.csv", sep = ";")
samples.out <- rownames(asv_table)
samples.out <- gsub("^((?:[^_]*_){2}).*", "\\1", basename(samples.out))
substr(samples.out,1,53)

"MI.M05812_0020.001.FLD0003.15172T14_R1_filt.fastq.gz"
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
