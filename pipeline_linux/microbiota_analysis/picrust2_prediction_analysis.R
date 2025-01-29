library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)

# Installing ggpicrust2 and other required libraries
{
# Had to install GSL via "sudo apt-get install libgsl-dev" from command line, for Aldex2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("phyloseq", "ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "BiocManager", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

# install.packages("devtools")
devtools::install_github("cafferychen777/ggpicrust2")
}

# Analysis of picrust2 results for Samuel
setwd("~/Documents/CHUM_git/Microbiota_17/picrust2/input/picrust2_out_pipeline/")

pathways <- read.table("pathways_out/path_abun_unstrat.tsv.gz", sep = "\t", header = TRUE)
ko_metagenome <- list(read.table("KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE),
                      read.table("KO_metagenome_out/seqtab_norm.tsv.gz", sep = "\t", header = TRUE),
                      read.table("KO_metagenome_out/weighted_nsti.tsv.gz", sep = "\t", header = TRUE))
                    
View(ko_metagenome[[1]]) # ko numbers
View(ko_metagenome[[2]]) # osef
View(ko_metagenome[[3]]) # osef


ec_metagenome <- list(read.table("EC_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE),
                      read.table("EC_metagenome_out/seqtab_norm.tsv.gz", sep = "\t", header = TRUE),
                      read.table("EC_metagenome_out/weighted_nsti.tsv.gz", sep = "\t", header = TRUE))

View(ec_metagenome[[1]]) # ec numbers
View(ec_metagenome[[2]]) # osef
View(ec_metagenome[[3]]) # osef


marker <- read.table("marker_predicted_and_nsti.tsv.gz", sep = "\t", header = TRUE)
View(marker)

ko_predicted <- read.table("KO_predicted.tsv.gz", sep = "\t", header = TRUE, row.names = 1)

# Remove the X at the start of the names
colnames(ko_metagenome[[1]]) = substring(colnames(ko_metagenome[[1]]), first = 2, 6)
colnames(ko_metagenome[[1]])[1] = "KO_id"

# Metadata handling
{
  metadata <- read.csv("~/Documents/CHUM_git/Microbiota_17/metadata/metadata.csv", sep = ";")
  
  #adding first col as rownames too
  rownames(metadata) <- metadata$sample_id
  
  #transforming week col from num to character
  metadata$week <- as.factor(metadata$week)
  metadata$diet <- as.factor(metadata$diet)
  
  # Replaces Samuel names by shorter versions
  metadata$treatment <- gsub(".*Putrescine.*", "Putrescine", metadata$treatment)
  metadata$treatment <- gsub(".*Vehicle.*", "Vehicle", metadata$treatment)
  metadata$genotype <- gsub(".*IL-22.*", "IL-22ra1-/-", metadata$genotype)
  
  # Creates gg_group specific to Claire
  metadata$gg_group[metadata$student == "Claire"] <- 
    paste(metadata$week[metadata$student == "Claire"], 
          metadata$diet[metadata$student == "Claire"], 
          sep = ":")
  
  # Creates gg_group specific to Samuel
  metadata$gg_group[metadata$student == "Samuel"] <- 
    paste(metadata$genotype[metadata$student == "Samuel"], 
          metadata$treatment[metadata$student == "Samuel"], 
          sep = ":")
}

# Only Samuel's metadata
meta <- metadata[metadata$student == "Samuel",]
meta <- meta[,-c(2,4)]
rownames(meta) <- 1:nrow(meta)

# Testing ggpicrust2 workflow
# Run ggpicrust2 with input data
{
results_data_input <- ggpicrust2(#file <-"KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
                                data = ko_metagenome[[1]],
                                 metadata = meta,
                                 group = "gg_group", # For example dataset, group = "Environment"
                                 pathway = "KO",
                                 daa_method = "ALDEx2",
                                 ko_to_kegg = TRUE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")

results_data_input[[1]]$plot
results_data_input[[2]]$plot
results_data_input[[1]]$results

data("metadata")
data("ko_abundance")
} # This is pretty much awful and does not make any sense at all, I guess we are doing it on our own
















# Select SCFA-related KOs
scfa_kos <- c("K00929", "K01034", "K00169")  # Example KOs for butyrate and propionate

# Extract SCFA KO abundances
scfa_abundance <- ko_metagenome[[1]][rownames(ko_metagenome[[1]]) %in% scfa_kos, ]

# Sum SCFA KO abundances for each sample
scfa_total <- as.data.frame(colSums(scfa_abundance))

# data <- scfa_abundance[1,]  %>% 
#   pivot_longer(
#     cols = 1:ncol(scfa_abundance), 
#     names_to = "sample",
#     values_to = "count"
#   )

# Metadata handling
{
metadata <- read.csv("~/Documents/CHUM_git/Microbiota_17/metadata/metadata.csv", sep = ";")

#adding first col as rownames too
rownames(metadata) <- metadata$sample_id

#transforming week col from num to character
metadata$week <- as.factor(metadata$week)
metadata$diet <- as.factor(metadata$diet)

# Replaces Samuel names by shorter versions
metadata$treatment <- gsub(".*Putrescine.*", "Putrescine", metadata$treatment)
metadata$treatment <- gsub(".*Vehicle.*", "Vehicle", metadata$treatment)
metadata$genotype <- gsub(".*IL-22.*", "IL-22ra1-/-", metadata$genotype)

# Creates gg_group specific to Claire
metadata$gg_group[metadata$student == "Claire"] <- 
  paste(metadata$week[metadata$student == "Claire"], 
        metadata$diet[metadata$student == "Claire"], 
        sep = ":")

# Creates gg_group specific to Samuel
metadata$gg_group[metadata$student == "Samuel"] <- 
  paste(metadata$genotype[metadata$student == "Samuel"], 
        metadata$treatment[metadata$student == "Samuel"], 
        sep = ":")
}

# Associate gg_group with each sample_id using metadata information
for(i in 1:nrow(scfa_total)){
  scfa_total$gg_group[i] <- metadata[metadata[rownames(metadata)==rownames(scfa_total)[i],"sample_id"],"gg_group"]
}

# Replacing a colname and adding a sample_id column, + gg_plot as a factor
colnames(scfa_total)[1] <- "scfa_ab"
scfa_total$sample_id <- rownames(scfa_total)
factor(scfa_total$gg_group, levels = c("Wt:Vehicle", "Wt:Putrescine", "IL-22ra1-/-:Vehicle", "IL-22ra1-/-:Putrescine"))

# Visualize with ggplot2
# df <- data.frame(Sample = names(scfa_total), SCFA_Potential = scfa_total)

p <- ggplot(scfa_total, aes(x = gg_group, y = scfa_ab, colors = gg_group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Predicted butyrate + propionate production", y = "Butyrate + propionate")
p
