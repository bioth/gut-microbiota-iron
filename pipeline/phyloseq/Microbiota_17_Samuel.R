#Loading required packages
{
  library(phyloseq)
  library(DESeq2) #for differential abundance analysis
  library(Biostrings)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(ggsignif) #for adding significance bars to a ggplot
  library(car)
  library(ggpubr)
  library(patchwork)
  library(cowplot) #for adding texts to ggplots
  library(vegan) #for beta diversity PCOAs and unifrac
  library(DECIPHER) #for performing multiple sequence alignments
  library(phangorn) #for building trees
  library(ape)
  library(tidyverse)
  library(writexl)
  library(openxlsx)
  library(reshape2)
  library(Hmisc)
  library(plotly) #To plot 3D pcoas
  library(StackbarExtended) # Thibault C. package
  library(ggh4x)
  
}


#Load custom functions for microbiota analysis
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/utilities.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/alpha_diversity_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/beta_diversity_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/correlation_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/relab_analysis_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/taxa_distrib_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/plot_microbiota_extension.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/plot_microbiota_ext.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/deseq2_log2fold_change_analysis.R")

# For microbiota 17 - specific to Samuel's data, and wt only
#set working directory
setwd("~/Documents/CHUM_git/Microbiota_17/")
asv_table <- as.data.frame(fread("asv_table/seqtab.nochim_run_m1.csv", sep = ";"))

# Set the first column as row names and remove it from the data frame
rownames(asv_table) <- asv_table[,1]  # Use the first column as row names
asv_table <- asv_table[,-1]  # Drop the first column

# Loading metadata of interest
{
metadata <- read.csv("metadata/metadata.csv", sep = ";")
rownames(metadata) <- metadata$sample_id
metadata$diet <- as.factor(metadata$diet)
metadata$treatment <- gsub(".*Putrescine.*", "Putrescine", metadata$treatment)
metadata$treatment <- gsub(".*Vehicle.*", "Vehicle", metadata$treatment)

# Creates gg_group specific to Samuel
metadata$gg_group[metadata$student == "Samuel"] <- 
  paste(metadata$genotype[metadata$student == "Samuel"], 
        metadata$treatment[metadata$student == "Samuel"], 
        sep = ":")
}

#load taxonomical assignments
taxa <- as.matrix(fread("taxonomy/taxa_annotation_m1.csv", sep = ";"))
rownames(taxa) <- taxa[,1]  # Use the first column as row names
taxa <- taxa[,-1]  # Drop the first column

# Load phylogenetic tree if possible
tree <- read.tree("~/Documents/CHUM_git/figures/old/samuel/beta_diversity/phylo_tree/phylogenetic_tree.newick")

#creating phyloseq object
ps <- phyloseq(otu_table(asv_table, taxa_are_rows = FALSE),
               sample_data(metadata),
               tax_table(taxa))

#use short names for the asvs (eg ASV21) rather than full dna sequence name
#though we can keep these dna sequences for other purposes (?)
#And adds refseq entry to the phyloseq object
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
names(dna) <- taxa_names(ps)

# Add tree to phyloseq object
ps <- merge_phyloseq(ps, phy_tree(tree))
sum(taxa_sums(ps)) #total number of reads
length(taxa_sums(ps)) #total number of ASVs

ps_samuel <- subset_samples(ps, student == "Samuel" & genotype == "Wt")
sum(taxa_sums(ps_samuel))
length(taxa_sums(ps_samuel))

# Filtering 
{
#function filtering out ASVs for which they were in total less than a threshold count
ps_samuel <- prune_taxa(taxa_sums(ps_samuel) > 0, ps_samuel) #Actual ASVs from Samuel's data
sum(taxa_sums(ps_samuel))
length(taxa_sums(ps_samuel))

ps_samuel <- prune_taxa(taxa_sums(ps_samuel) > 10, ps_samuel)
sum(taxa_sums(ps_samuel))
length(taxa_sums(ps_samuel))

#Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
ps_samuel <- prune_taxa(colSums(otu_table(ps_samuel) > 0) >= (0.05 * nsamples(ps_samuel)), ps_samuel)
sum(taxa_sums(ps_samuel))
length(taxa_sums(ps_samuel))
}

#for Samuel's data, put treatment as factor and define order
sample_data(ps_samuel)$treatment <- factor(sample_data(ps_samuel)$treatment, levels = c("Vehicle", "Putrescine")) # Vehicle as reference
customColors = list(list('black','#A22004'))
# Beta diversity
{
  betaDiversityPairwise(ps_samuel, "treatment", pairs = list(list("Vehicle","Putrescine")), "bray", customColors, font = "Times New Roman", displayPValue = FALSE, dim = c(3.75,3.75), transform = "none", path = "~/Documents/CHUM_git/figures/Samuel_final_wt/beta_diversity/")
  betaDiversityPairwise(ps_samuel, "treatment", pairs = list(list("Vehicle","Putrescine")), "wunifrac", customColors, font = "Times New Roman", displayPValue = FALSE, dim = c(3.75,3.75), transform = "none", path = "~/Documents/CHUM_git/figures/Samuel_final_wt/beta_diversity/")
}

# Differential abundance analysis
{
  
}
                                           