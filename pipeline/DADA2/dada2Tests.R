library("dada2"); packageVersion("dada2")

#Setting working directory
setwd("D:/CHUM_git/16s_data/trimmed_fastq/")

#Listing files in the current working directory
list.files()

#Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
r1_fastq <- sort(list.files(pattern="_R1_trimmed.fastq", full.names = TRUE))
r2_fastq <- sort(list.files(pattern="_R2_trimmed.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample_names <- gsub("^((?:[^_]*_){2}).*", "\\1", basename(r1_fastq))
# Remove the last character
sample_names <- substr(sample_names, 1, nchar(sample_names) - 1)


plotQualityProfile(r1_fastq[1:2])


#Quality filtering of the reads

# Place filtered files in filtered/ subdirectory
path = "../dada2_filtered_and_trimmed"

filtR1 <- file.path(path, paste0(sample_names, "_R1_filt.fastq.gz"))
filtR2 <- file.path(path, paste0(sample_names, "_R2_filt.fastq.gz"))
names(filtR1) <- sample_names
names(filtR2) <- sample_names

out <- filterAndTrim(r1_fastq, filtR1, r2_fastq, filtR2, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
out
