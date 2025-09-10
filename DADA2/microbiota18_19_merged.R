library("dada2"); packageVersion("dada2")
library(ShortRead)
library(data.table)

# Function checking if a dir exists and creating it otherwise
existingDirCheck <- function(path){
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Directory created: ", path)
  } else {
    message("Directory already exists: ", path)
  }
  
}

# Listing files in the current working directory
setwd("~/Documents/CHUM_git/Microbiota_18/trimmed_fastq/")
list.files()

pattern <- paste(c("T0","T35","T49"), collapse = "|") # Select only diet timepoints

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
r1_fastq_M18 <- sort(list.files(pattern="R1_trimmed.fastq", full.names = TRUE)) 
r2_fastq_M18 <- sort(list.files(pattern="R2_trimmed.fastq", full.names = TRUE))

# Keep only samples of interest
r1_fastq_M18 <- r1_fastq_M18[grepl(pattern, r1_fastq_M18)]
r2_fastq_M18 <- r2_fastq_M18[grepl(pattern, r2_fastq_M18)]

# Inspect if we should apply similar trimming parameters 
plotQualityProfile(r1_fastq_M18, aggregate = TRUE)
plotQualityProfile(r2_fastq_M18, aggregate = TRUE)

# Listing files in the current working directory
setwd("~/Documents/CHUM_git/Microbiota_19/trimmed_fastq/")
list.files()

pattern <- paste(c("T0","T35","T49"), collapse = "|") # Select only diet timepoints

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
r1_fastq_M19 <- sort(list.files(pattern="R1_trimmed.fastq", full.names = TRUE)) 
r2_fastq_M19 <- sort(list.files(pattern="R2_trimmed.fastq", full.names = TRUE))

# Keep only samples of interest
r1_fastq_M19 <- r1_fastq_M19[grepl(pattern, r1_fastq_M19)]
r2_fastq_M19 <- r2_fastq_M19[grepl(pattern, r2_fastq_M19)]

# Inspect if we should apply similar trimming parameters 
plotQualityProfile(r1_fastq_M19, aggregate = TRUE)
plotQualityProfile(r2_fastq_M19, aggregate = TRUE)



# Find optimal trimming parameters for Microbiota_19
#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample_names <- gsub(".*\\.(.*?)_R1_trimmed\\.fastq.gz$", "\\1", basename(r1_fastq_M19))

# Quality filtering of the reads
# Place filtered files in filtered/ subdirectory
existingDirCheck("../../Microbiota_18_19_merged/dada2_filtered_and_trimmed")
path = "../../Microbiota_18_19_merged/dada2_filtered_and_trimmed"

filtR1 <- file.path(path, paste0(sample_names, "_R1_filt.fastq.gz"))
filtR2 <- file.path(path, paste0(sample_names, "_R2_filt.fastq.gz"))
names(filtR1) <- sample_names
names(filtR2) <- sample_names

# Filter and trim
out <- filterAndTrim(r1_fastq_M19, filtR1, r2_fastq_M19, filtR2, truncLen=c(250,230),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, matchIDs=TRUE) 

#Check how many reads made it out
print(out)
print("Mean number of reads that passed filter:")
print(mean(out[,2]/out[,1])*100)
print("Min:")
print(min(out[,2]/out[,1])*100)
print("Max:")
print(max(out[,2]/out[,1])*100)

# Load asv table for each microbiota 18 and microbiota 19 runs
asv_table_M18 <- as.matrix(read.csv("~/Documents/CHUM_git/Microbiota_18_19_merged/asv_table/asv_table_M18.csv", sep = ";"))
asv_table_M18 <- as.data.frame(fread("~/Documents/CHUM_git/Microbiota_18_19_merged/asv_table/asv_table_M18.csv", sep = ";"))
rownames(asv_table_M18) <- asv_table_M18[,1]  # Use the first column as row names
asv_table_M18 <- asv_table_M18[,-1]  # Drop the first column

asv_table_M19 <- as.matrix(read.csv("~/Documents/CHUM_git/Microbiota_18_19_merged/asv_table/asv_table_M19.csv", sep = ";"))
asv_table_M19 <- as.matrix(read.csv("~/Documents/CHUM_git/Microbiota_18_19_merged/asv_table/asv_table_M19_2.csv", sep = ";"))
asv_table_M19 <- as.data.frame(fread("~/Documents/CHUM_git/Microbiota_18_19_merged/asv_table/asv_table_M19_2.csv", sep = ";"))
rownames(asv_table_M19) <- asv_table_M19[,1]  # Use the first column as row names
asv_table_M19 <- asv_table_M19[,-1]  # Drop the first column

asv_table_M18 <- as.matrix(asv_table_M18)
asv_table_M19 <- as.matrix(asv_table_M19)

# Inspect distribution of sequence lengths
print("Distribution of sequence lengths:")
table(nchar(getSequences(asv_table_M18)))
table(nchar(getSequences(asv_table_M19)))

# Merge the ASV table by summing similar sequences ids
asv_table <- mergeSequenceTables(asv_table_M18, asv_table_M19, repeats="sum")

# Removing chimeras
seqtab.nochim <- removeBimeraDenovo(asv_table, method="consensus", multithread=TRUE, verbose=TRUE,
                                    allowOneOff=FALSE, minFoldParentOverAbundance=8)

seqtab.nochim <- collapseNoMismatch(seqtab.nochim)

ncol(seqtab.nochim)
# What percentage of chimeras over the total dataset
print("Percentage of chimeras over the total dataset:")
1-sum(seqtab.nochim)/sum(asv_table)

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/CHUM_git/training_set/silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)
taxa_w_species <- addSpecies(taxtab = taxa, refFasta = "~/Documents/CHUM_git/training_set/silva_species_assignment_v138.1.fa.gz")

existingDirCheck("~/Documents/CHUM_git/Microbiota_18_19_merged/taxonomy/")
write.table(seqtab.nochim, sep = ";", file = "~/Documents/CHUM_git/Microbiota_18_19_merged/asv_table/merged_asv_table.csv", col.names = TRUE)
write.table(taxa_w_species, sep = ";", file = "~/Documents/CHUM_git/Microbiota_18_19_merged/taxonomy/merged_taxa_annotation.csv", col.names = TRUE)
