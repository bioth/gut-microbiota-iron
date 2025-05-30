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
unique(data$`Séquence de l'amorce sens`)
unique(data$`Séquence de l'amorce antisens`)
unique(data$`Adaptateur Read 1 (NOTE: Usage restreint par le Disclaimer Illumina visible dans la page Projet de Nanuq)`)
unique(data$`Adaptateur Read 2 (NOTE: Usage restreint par le Disclaimer Illumina visible dans la page Projet de Nanuq)`)
for (i in 1:ncol(data)){
  if (nrow(unique(data[,i]))>=2){
    print(nrow(unique(data[,i])))
    print(unique(data[,i]))
  }
}
#everything is good, don't have to perform demultiplexing so we can move on to get rid
#of the non-biological bases. Thibault C. suggested that we might use cutAdapt as different adapters
#were used in the dataset but it does not seem to be the case
