library("dada2"); packageVersion("dada2")
library(ggplot2)

#microbiota 17
#subset working
#Setting working directory and selecting loading files of interest
setwd("~/Documents/CHUM_git/Microbiota_17/trimmed_fastq/")

#Listing files in the current working directory
list.files()

#Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
r1_fastq <- sort(list.files(pattern="R1_trimmed.fastq", full.names = TRUE))
r2_fastq <- sort(list.files(pattern="R2_trimmed.fastq", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample_names <- gsub(".*\\.(.*?)_R1_trimmed\\.fastq$", "\\1", basename(r1_fastq))

#Printing aggregate quality profile plot for all of the R1 files
plotQualityProfile(r1_fastq[1:132], aggregate = TRUE)
ggsave(dpi = 300, filename = "../r_console_output/quality_profileR1.png", height = 6, width = 10)

#Printing aggregate quality profile plot for all of the R2 files
plotQualityProfile(r2_fastq[1:132], aggregate = TRUE)
ggsave(dpi = 300, filename = "../r_console_output/quality_profileR2.png", height = 6, width = 10)

# Start redirecting output to a file
run <- stringr::str_replace_all(as.character(Sys.time()), "[ :]", "-")
dataset <- "Microbiota17"
run <- paste("../r_console_output/", run, dataset,".txt", sep = '')
sink(run)


#Quality filtering of the reads
# Place filtered files in filtered/ subdirectory
# Specify the name of the new directory
new_dir <- "dada2_filtered_and_trimmed"

# Check if the directory already exists
if (!file.exists(new_dir)) {
  # If it doesn't exist, create it
  dir.create(new_dir)
  cat("Directory created:", new_dir)
} else {
  cat("Directory already exists:", new_dir)
}
path = "../dada2_filtered_and_trimmed"

filtR1 <- file.path(path, paste0(sample_names, "_R1_filt.fastq.gz"))
filtR2 <- file.path(path, paste0(sample_names, "_R2_filt.fastq.gz"))
names(filtR1) <- sample_names
names(filtR2) <- sample_names

out <- filterAndTrim(r1_fastq, filtR1, r2_fastq, filtR2, truncLen=c(FALSE,FALSE),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, matchIDs=TRUE) # On Windows set multithread=FALSE

#Check how many reads made it out
print(out)
print("Mean number of reads that passed filter:")
print(mean(out[,2]/out[,1])*100)
print("Min:")
print(min(out[,2]/out[,1])*100)
print("Max:")
print(max(out[,2]/out[,1])*100)

#dereplication of the fastq files
derepR1 <- derepFastq(filtR1)
derepR2 <- derepFastq(filtR2)
names(derepR1) <- sample_names
names(derepR2) <- sample_names


#error rates estimation
errR1 <- learnErrors(derepR1, multithread=TRUE, randomize = TRUE)
errR2 <- learnErrors(derepR2, multithread=TRUE, randomize = TRUE)
plotErrors(errR1, nominalQ=TRUE)
ggsave(dpi = 300, filename = "../r_console_output/error_plotR1.png", height = 10, width = 10)
plotErrors(errR2, nominalQ=TRUE)
ggsave(dpi = 300, filename = "../r_console_output/error_plotR2.png", height = 10, width = 10)


#running the inference algorithm => based on trimmed/filtered fastq and error rates
dadaR1 <- dada(derepR1, err=errR1, multithread=TRUE, pool = "pseudo")
dadaR2 <- dada(derepR2, err=errR2, multithread=TRUE, pool = "pseudo")

#merge paired reads 
mergers <- mergePairs(dadaR1, derepR1, dadaR2, derepR2, verbose=TRUE)

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
setwd("../asv_table/")
write.table(seqtab.nochim, sep = ";", file = "seqtab.nochim_run.csv", col.names = TRUE)



#tracking what reads made it through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFsSamuel, getN), sapply(dadaRsSamuel, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names
head(track)
track

mean(track[,4]/track[,1])*100
max(track[,4]/track[,1])*100
min(track[,4]/track[,1])*100



#extracting taxonomical information
#Loading ASV table previously generated
setwd("D:/CHUM_git/Microbiota_17/asv_table/")
asv_table <- read.csv("seqtab.nochim_run2.csv", sep = ";")

#transforming asv_table into matrix so that it can be used by dada2 taxonomic assignment algorithm
asv_table <- as.matrix(seqtab.nochim)
taxa <- assignTaxonomy(asv_table, "../../training_set/silva_nr99_v138.1_wSpecies_train_set.fa.gz")
taxa.print <- taxa #removing rownames for display
rownames(taxa.print) <- NULL
head(taxa.print)


#save taxa matrix so that we can use it later
write.table(taxa, sep = ";", "../taxonomy/taxa_annotation.csv", col.names = TRUE)


sink()


