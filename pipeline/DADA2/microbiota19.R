library("dada2"); packageVersion("dada2")
library(ggplot2)
library(dplyr)
library(readxl)

# Function checking if a dir exists and creating it otherwise
existingDirCheck <- function(path){
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Directory created: ", path)
  } else {
    message("Directory already exists: ", path)
  }
  
}

# Setting working directory and selecting loading files of interest
setwd("~/Documents/CHUM_git/Microbiota_19/trimmed_fastq/")

# Loading metdata
meta <- read_excel("../metadata/dissection.xlsx")
samples <- substring(meta$ID, 1, 5)

# Listing files in the current working directory
list.files()

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
r1_fastq <- sort(list.files(pattern="R1_trimmed.fastq", full.names = TRUE)) 
r2_fastq <- sort(list.files(pattern="R2_trimmed.fastq", full.names = TRUE))

r1_fastq <- r1_fastq[grepl(paste0(paste0(samples, "s*"), collapse = "|"), r1_fastq)]
r2_fastq <- r2_fastq[grepl(paste0(paste0(samples, "s*"), collapse = "|"), r2_fastq)]

# Plot quality profiles
{
  plotQualityProfile(r1_fastq, aggregate = TRUE)
  plotQualityProfile(r2_fastq, aggregate = TRUE)
}



# Transforming asv_table into matrix so that it can be used by dada2 taxonomic assignment algorithm
asv_table <- as.matrix(read.csv("~/Documents/CHUM_git/Microbiota_19/asv_table/asv_table.csv", sep = ";"))
taxa <- assignTaxonomy(asv_table, "~/Documents/CHUM_git/training_set/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread = TRUE)
taxa.print <- taxa #removing rownames for display
rownames(taxa.print) <- NULL
head(taxa.print)

# Save taxa matrix so that we can use it later
existingDirCheck("~/Documents/CHUM_git/Microbiota_19/taxonomy")
write.table(taxa, sep = ";", file = "~/Documents/CHUM_git/Microbiota_19/taxonomy/taxa_annotation.csv", col.names = TRUE)
