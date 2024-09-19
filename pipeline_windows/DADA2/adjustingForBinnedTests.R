library("dada2"); packageVersion("dada2")
library(ggplot2)


#Setting working directory
setwd("D:/CHUM_git/Microbiota_17/dada2_filtered_and_trimmed/")

#list filtered trimmed fastq files (passed the first two steps of the pipeline)
list.files()

#printing files names
filtR1 <- sort(list.files(pattern="_R1_filt.fastq.gz", full.names = TRUE))
filtR2 <- sort(list.files(pattern = "_R2_filt.fastq.gz", full.names = TRUE))

#create subset of filtR1 and filtR2 to have results faster
filtR1_subset <- filtR1[1:3]
filtR2_subset <- filtR2[1:3]

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
errR1 <- learnErrors(filtR1_subset, multithread=FALSE)
errR2 <- learnErrors(filtR2_subset, multithread=FALSE)
plotErrors(errR1, nominalQ=TRUE)
plotErrors(errR2, nominalQ=TRUE)