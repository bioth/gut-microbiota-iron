library("dada2"); packageVersion("dada2")

#Setting working directory
setwd("D:/CHUM_git/16s_data_test/trimmed_fastq/")

#Listing files in the current working directory
list.files()

#Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
r1_fastq <- sort(list.files(pattern="_R1_trimmed.fastq", full.names = TRUE))
r2_fastq <- sort(list.files(pattern="_R2_trimmed.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample_names <- gsub("^((?:[^_]*_){2}).*", "\\1", basename(r1_fastq))
# Remove the last character
sample_names <- substr(sample_names, 1, nchar(sample_names) - 1)


plotQualityProfile(r1_fastq)
plotQualityProfile(r2_fastq)


#Quality filtering of the reads, merging of reverse and forward reads
# Place filtered files in filtered/ subdirectory
path = "../dada2_filtered_and_trimmed"

filtR1 <- file.path(path, paste0(sample_names, "_R1_filt.fastq.gz"))
filtR2 <- file.path(path, paste0(sample_names, "_R2_filt.fastq.gz"))
names(filtR1) <- sample_names
names(filtR2) <- sample_names

out <- filterAndTrim(r1_fastq, filtR1, r2_fastq, filtR2, truncLen=c(220,200),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
out


#error rates estimation
#change working directory
setwd("D:/CHUM_git/16s_data_test/dada2_filtered_and_trimmed/")
list.files()

#load filtered files
filtR1 <- sort(list.files(pattern="_R1_filt.fastq.gz", full.names = TRUE))
filtR2 <- sort(list.files(pattern="_R2_filt.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample_names <- gsub("^((?:[^_]*_){2}).*", "\\1", basename(filtR1))
# Remove the last character
sample_names <- substr(sample_names, 1, nchar(sample_names) - 1)

errR1 <- learnErrors(filtR1, multithread=TRUE)
errR2 <- learnErrors(filtR2, multithread=TRUE)
plotErrors(errR1, nominalQ=TRUE)
plotErrors(errR2, nominalQ=TRUE)


#sample inference (core DADA2 algorithm)
dadaFs <- dada(filtR1, err=errR1, multithread=TRUE)
dadaRs <- dada(filtR2, err=errR2, multithread=TRUE)

#merge paired reads
merge_data <- mergePairs(dadaFs, filtR1, dadaRs, filtR2, verbose=TRUE)






#comparing size of the files untrimmed

#Setting working directory
setwd("D:/CHUM_git/16s_data_test/uncompressed_files/")

#Listing files in the current working directory
list.files()

#Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
r1_fastq <- sort(list.files(pattern="R1.fastq", full.names = TRUE))
r2_fastq <- sort(list.files(pattern="R2.fastq", full.names = TRUE))

r1_fastq[1]


plotQualityProfile(r1_fastq)
plotQualityProfile(r2_fastq)
