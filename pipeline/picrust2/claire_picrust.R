library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(readxl)
library(openxlsx)

source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/utilities.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/picrust2/picrust2_utilities.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/other scripts/dataManipFunctions.R")

# Metadata handling
{
  #loading metadata of interest
  metadata <- read.csv("~/Documents/CHUM_git/Microbiota_17/metadata/metadata.csv", sep = ";")
  metadata <- metadata[metadata$student == "Claire",]
  metadata <- metadata[metadata$week == "10",]
  
  # Adding id col as rownames too
  rownames(metadata) <- metadata$sample_id
  
  # Putting diet as factor
  metadata$diet <- factor(metadata$diet, levels = c("50", "500"))
}

# Analysis of picrust2 results at 10 weeks / tryptophan metabolism
setwd("~/Documents/CHUM_git/Microbiota_17/claire_picrust/picrust2/picrust2_out_pipeline/")
ko_df <- read.table("KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE) # Load KO annotations
colnames(ko_df)[2:ncol(ko_df)] <- substring(colnames(ko_df)[2:ncol(ko_df)], 2)
pattern <- paste(metadata$sample_id, collapse = "|")
indexes <-  grep(pattern, colnames(ko_df)) 
ko_df <- ko_df[,c(1,indexes)] # Keep only samples for 10 weeks
# K01667 (tryptophan synthase alpha), K11819 (kynureninase)
ko_df <- ko_df[ko_df$function. %in% c("K01667","K11819"),] # Look at KOs of interest
ko_df <- ko_df[ko_df$function. %in% c("K01695"),]
ko_df <- ko_df[ko_df$function. %in% c("K00453","K00459","K01465","K00480","K01440"),]
rownames(ko_df) <- ko_df[,1]
ko_df <- ko_df[,-1]
relative_trypto_ko <- apply(ko_df, 2, function(x) x / sum(x))
trypto_total <- colSums(ko_df)
ko_df <- pivot_longer(ko_df, cols = c(2:21), names_to = "sample_id", values_to = "K01667") # Put in long format 
ko_df <- merge(ko_df, metadata, by = "sample_id")


ironBoxplot(ko_df, measure = "K01667", group = "diet",
            title = "Tryptophan metabolism at 10 weeks",
            y_axis_title = "prediceted abundance", custom_colors = c("blue","red"))

verifyStatsAssumptions(ko_df, group = "diet", measure = "K01667")
wilcox.test(K01667 ~ diet, data = ko_df)
t.test(K01667 ~ diet, data = ko_df, var.equal = TRUE)





