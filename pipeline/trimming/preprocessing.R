library("dada2")
library(readxl)

setwd("D:/CHUM_git/16s_data/")
#checking data



#In uncompressed files, they are 468 files, removing the first two, that's 466 fastq
#for both directions, meaning 233 uniques IDs


metadata <- read.csv("metadata.csv", header = TRUE, sep = ";")
length(unique(metadata$SampleID))
duplicated(metadata$SampleID)
#metadata.csv contains 230 unique samples IDS, removing last one which corresponds to nothing
#that's 229 

#thus how do we explain the differences between the number of samples in metadata and
#the number of samples associated with the fastq files = let's check the informations
#about the things sequenced

data <- read_excel("Microbiota_10_MiSeqReadSet_2020-12-08.xlsx")
length(unique(data$Nom)) #229 too if you remove negative control
data$Nom[duplicated(data$Nom)] #"15195T21" "14465"    "15172T14"
length(unique(data$`Identifiant de lecture`)) #233


#everything is good, don't have to perform demultiplexing so we can move on to get rid
#of the non-biological bases
