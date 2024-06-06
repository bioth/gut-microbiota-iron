library("dada2"); packageVersion("dada2")
library(ShortRead)
library(ggplot2)
library(reshape2)

#Setting working directory and selecting sample ids of interest
setwd("D:/CHUM_git/16s_data/metadata/")
metadata <- read.csv("metadata.csv", sep = ";")
ids <- metadata[1:180,5]

#loading files
setwd("../trimmed_fastq/")

#Listing files in the current working directory
list.files()

#Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
r1_fastq <- sort(list.files(pattern="T.+_R1_trimmed.fastq", full.names = TRUE))
# Create a pattern that matches any of the IDs
pattern <- paste(ids, collapse = "|")
# Find the file names that contain any of the IDs
r1_fastq <- r1_fastq[grepl(pattern, r1_fastq)]

r2_fastq <- sort(list.files(pattern="_R2_trimmed.fastq", full.names = TRUE))
# Find the file names that contain any of the IDs
r2_fastq <- r2_fastq[grepl(pattern, r2_fastq)]

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample_names <- gsub("^((?:[^_]*_){2}).*", "\\1", basename(r1_fastq))
# Remove the last character
sample_names <- substr(sample_names, 1, nchar(sample_names) - 1)

#Quality filtering of the reads
# Place filtered files in filtered/ subdirectory
path = "../dada2_filtered_and_trimmed"

filtR1 <- file.path(path, paste0(sample_names, "_R1_filt.fastq.gz"))
filtR2 <- file.path(path, paste0(sample_names, "_R2_filt.fastq.gz"))
names(filtR1) <- sample_names
names(filtR2) <- sample_names

out <- filterAndTrim(r1_fastq, filtR1, r2_fastq, filtR2, truncLen=c(225,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, matchIDs=TRUE) # On Windows set multithread=FALSE
#Check how many reads made it out
head(out)
print(mean(out[,2]/out[,1])*100)
print(min(out[,2]/out[,1])*100)
print(max(out[,2]/out[,1])*100)


#error rates estimation
#change working directory
setwd("D:/CHUM_git/16s_data/dada2_filtered_and_trimmed/")
list.files()

#load filtered files
filtR1 <- sort(list.files(pattern="_R1_filt.fastq.gz", full.names = TRUE))
filtR1 <- filtR1[grepl(pattern, filtR1)]
filtR2 <- sort(list.files(pattern="_R2_filt.fastq.gz", full.names = TRUE))
filtR2 <- filtR2[grepl(pattern, filtR2)]

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample_names <- gsub("^((?:[^_]*_){2}).*", "\\1", basename(filtR1))
# Remove the last character
sample_names <- substr(sample_names, 1, nchar(sample_names) - 1)

errR1 <- learnErrors(filtR1, multithread=FALSE)
errR2 <- learnErrors(filtR2, multithread=FALSE)
plotErrors(errR1, nominalQ=TRUE)
plotErrors(errR2, nominalQ=TRUE)


#running the inference algorithm => based on trimmed/filtered fastq and error rates
dadaFs <- dada(filtR1, err=errR1, multithread=FALSE)
dadaRs <- dada(filtR2, err=errR2, multithread=FALSE)

#merge paired reads 
mergers <- mergePairs(dadaFs, filtR1, dadaRs, filtR2, verbose=TRUE)

#construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#what percentage of chimears over the total dataset
sum(seqtab.nochim)/sum(seqtab)


#trying to export the data
setwd("../")
rownames(seqtab.nochim)[1]
write.table(seqtab.nochim, sep = ";", file = "seqtab.nochim_run2.csv", col.names = TRUE)



#tracking what reads made it through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names
head(track)
track

mean(track[,6]/track[,1])*100
max(track[,6]/track[,1])*100
min(track[,6]/track[,1])*100
