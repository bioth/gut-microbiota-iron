library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(car)
library(ggprism)
library(patchwork)
library(readxl)
library(openxlsx)
library(compositions) # For CLR transformation

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
ko_df <- read.table("KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE) # Load KO annotation
ko_df <- read.table("EC_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE) # Load KO annotationss
colnames(ko_df)[2:ncol(ko_df)] <- substring(colnames(ko_df)[2:ncol(ko_df)], 2)
pattern <- paste(metadata$sample_id, collapse = "|")
indexes <-  grep(pattern, colnames(ko_df)) 
ko_df <- ko_df[,c(1,indexes)] # Keep only samples for 10 weeks
# K01667 (tryptophan synthase alpha), K11819 (kynureninase)
ko_df <- ko_df[ko_df$function. %in% c("K01667","K11819"),] # Look at KOs of interest
ko_df <- ko_df[ko_df$function. %in% c("K01695"),]
ko_df <- ko_df[ko_df$function. %in% c("K00453","K00459","K01465","K00480","K01440"),]
ko_df <- ko_df[ko_df$function. %in% c("K01667","K01593","K00466","K00453","K00444",
                                      "K00451","	K01555","K00832","K00141","K00016"),] # exclided KO0274 => serotonin pathway
ko_df <- ko_df[ko_df$function. %in% c("K00832","K00466","K01667","K00463","K01593"),]
ko_df <- ko_df[ko_df$function. %in% c("EC:2.6.1.57","EC:1.13.12.3","EC:4.1.99.1","EC:1.13.11.52","EC:4.1.1.28"),]

test = pivot_longer(ko_df, cols = -function., names_to = "sample_id", values_to = "abundance")
test = merge(test, metadata, by = "sample_id")

ggplot(test, aes(x = diet, y = abundance, color = diet))+
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), alpha = 0.5) +
  stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..-0.25, yend=..y..), size = 1)+ #adding horizontal bars representing means
  stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..+0.25, yend=..y..), size = 1)+ 
  stat_summary(fun.data="mean_cl_normal", geom="errorbar", width=0.2, size = 0.7) + #adding SEM error bars
  facet_wrap(~ function., scales = "free_y")+
  scale_color_manual(values = c("blue","red"))
  labs(x = "", y = "Abundance", title = "Abundance per KO across samples")
  
  # Lactate dehydrogenase (fldH)	K00016	IPYA → ILA
  # Monoamine oxidase (MAO)	K00274	Serotonin → 5-HIAA (tentative) uhuh
  # IDO / TDO	K00453 / K00444	Tryptophan → Kynurenine
  # Tryptophanase (TNA)	K01667	Tryptophan → Indole
  # Kynurenine aminotransferase (KAT)	K00832	Kynurenine → Kynurenic acid
  
verifyStatsAssumptions(test[test$function. == "K00832",], group = "diet", measure = "abundance")
t.test(abundance ~ diet, data = test[test$function. == "K00832",], var.equal = TRUE)

verifyStatsAssumptions(test[test$function. == "K01667",], group = "diet", measure = "abundance")
t.test(abundance ~ diet, data = test[test$function. == "K01667",], var.equal = TRUE)

rownames(ko_df) <- ko_df[,1]
ko_df <- ko_df[,-1]

clr_data <- as.data.frame(clr(t(ko_df + 1e-6))) # CLR transform with pseudocount
clr_data$score <- rowMeans(clr_data) # Summary score per sample
clr_data$sample_id <- rownames(clr_data)
clr_data <- merge(clr_data, metadata, by = "sample_id")


ironBoxplot(clr_data, measure = "score", group = "diet",
            title = "Tryptophan metabolism at 10 weeks",
            y_axis_title = "prediceted abundance", custom_colors = c("blue","red"))

verifyStatsAssumptions(clr_data, group = "diet", measure = "score")
wilcox.test(K01667 ~ diet, data = ko_df)
t.test(score ~ diet, data = clr_data, var.equal = TRUE)




# Load pathways abundances
pwy <- read.table("pathways_out/path_abun_unstrat.tsv", sep = "\t", header = TRUE) # Load pathways abundances
colnames(pwy)[2:ncol(pwy)] <- substring(colnames(pwy)[2:ncol(pwy)], 2)
trpPwys <- c("TRPSYN-PWY","TRPCAT-PWY","TRPIAACAT-PWY","PWY-5655","TRYPTOPHAN-RXN")
pwy <- pwy[pwy$pathway %in% trpPwys,]
test = pivot_longer(pwy, cols = -pathway, names_to = "sample_id", values_to = "abundance")
test = merge(test, metadata, by = "sample_id")
verifyStatsAssumptions(test, group = "diet", measure = "abundance")
t.test(abundance ~ diet, data = test, var.equal = FALSE)
wilcox.test(K01667 ~ diet, data = ko_df)

ironBoxplot(test, measure = "abundance", group = "diet",
            title = "Tryptophan metabolism at 10 weeks",
            y_axis_title = "prediceted abundance", custom_colors = c("blue","red"))



