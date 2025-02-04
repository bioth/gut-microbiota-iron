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
  library(ggpattern)
  
}

#To unload all packages
{
  # Get all loaded packages
  loaded_pkgs <- setdiff(names(sessionInfo()$otherPkgs), "base")
  
  # Unload each package
  lapply(loaded_pkgs, function(pkg) {
    detach(paste("package", pkg, sep = ":"), character.only = TRUE, unload = TRUE)
  })
}

#Load custom functions for microbiota analysis
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/utilities.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/alpha_diversity_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/beta_diversity_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/correlation_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/relab_analysis_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/taxa_distrib_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/plot_microbiota_ext.R")

#for microbiota 18
#set working directory
setwd("~/Documents/CHUM_git/Microbiota_18/")
asv_table <- as.data.frame(fread("asv_table/asv_table_server.csv", sep = ";"))

# Set the first column as row names and remove it from the data frame
rownames(asv_table) <- asv_table[,1]  # Use the first column as row names
asv_table <- asv_table[,-1]  # Drop the first column

# Metadata handling
{
#loading metadata of interest
metadata <- read.csv("metadata/metadata.csv", sep = ";")
  
# Remove the non-metadata stuff (liver measures and stuff)
metadata <- metadata[,-c(5:8)]
  
# Remove the letter at the end of id
metadata$id <- substring(metadata$id, 1, 5)

#adding id col as rownames too
rownames(metadata) <- metadata$id

# Remove dead mouse
metadata <- metadata[-46,]

# Extract 16S reads sample ids
samples <- read.xlsx("metadata/Microbiota_18_samples_2025-01-13.xlsx")
samples <- as.data.frame(samples$Nom)
colnames(samples) <- "full_id"
samples$id <- substring(samples$full_id, 1, 5)
samples$timepoint <- substring(samples$full_id, 8, nchar(samples$full_id)) 

# Bind both metadata df to link timepoints with their metadata (diet and treatment)
metadata <- merge(samples, metadata, by = "id")

# Consider timepoint 53 similar as timepoint 54
metadata[metadata$timepoint=="53","timepoint"] <- "54"

# Add week column 
metadata$week <- ifelse(metadata$timepoint == "final", "18", as.character(round(as.numeric(metadata$timepoint)/7, 1)+3))

# Adding gg_group variable (combination of time, diet and treatment)
metadata$gg_group <- 
  paste(metadata$timepoint, 
        metadata$diet,
        metadata$treatment, 
        sep = ":")

# Another gg_group variable (diet and treatment only)
metadata$gg_group2 <- 
  paste(metadata$diet,
        metadata$treatment, 
        sep = ":")

# Put full_id as rownames
rownames(metadata) <- metadata$full_id
}

#load taxonomical assignments
taxa <- as.matrix(fread("taxonomy/taxa_annotation_server.csv", sep = ";"))

# Set the first column as row names and remove it from the data frame
rownames(taxa) <- taxa[,1]  # Use the first column as row names
taxa <- taxa[,-1]  # Drop the first column

{
# Load phylogenetic tree if possible
tree <- read.tree("~/Documents/CHUM_git/Microbiota_18/taxonomy/phylogenetic_tree.newick")
}
  
#creating phyloseq object
ps <- phyloseq(otu_table(asv_table, taxa_are_rows = FALSE),
               tax_table(taxa), sample_data(metadata))


#use short names for the asvs (eg ASV21) rather than full dna sequence name
#though we can keep these dna sequences for other purposes (?)
#And adds refseq entry to the phyloseq object
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
names(dna) <- taxa_names(ps)

#Creating phylogenetic tree
{
  #Rarefaction curve, not necessary when lots of high quality data such as here
  {
    tab = otu_table(ps)
    class(tab) <- "matrix" # as.matrix() will do nothing
    tab <- t(tab)
    rarecurve(tab, step=10000, lwd=2, ylab="ASV",  label=F)
  }
  
  #adding a phylogenetic tree
  # Align sequences
  alignment <- AlignSeqs(DNAStringSet(dna), processors = 20)
  
  # Build a distance matrix
  dist_matrix <- DistanceMatrix(alignment, processors = 20)
  
  #build tree using neighbor-joining method
  tree  <- nj(dist_matrix)
  
  # Export the tree as a Newick file
  write.tree(tree, file = "~/Documents/CHUM_git/Microbiota_18/taxonomy/phylogenetic_tree.newick")
  
  #refinement with maximum likelihood
  {
    # Convert the alignment to a phangorn-compatible format
    alignment_matrix <- as.matrix(alignment)
    
    # Estimate the substitution model (e.g., GTR)
    fitGTR <- pml(initial_tree, data = alignment_matrix)
    
    # Optimize the ML tree
    ml_tree <- optim.pml(fitGTR, model = "GTR", rearrangement = "stochastic")
  }
}

# Add tree to phyloseq object
{
ps <- merge_phyloseq(ps, phy_tree(tree))
}

sum(taxa_sums(ps)) #total number of reads
length(taxa_sums(ps)) #total number of ASVs

View(tax_table(ps))
nrow(tax_table(ps))-sum(is.na(tax_table(ps)[,7])) #how many species detected

# Function filtering out ASVs for which they were in total less than a threshold count
ps_flt <- prune_taxa(taxa_sums(ps) > 10, ps)
sum(taxa_sums(ps_flt))
length(taxa_sums(ps_flt))

# Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
ps_flt <- prune_taxa(colSums(otu_table(ps_flt) > 0) >= (0.05 * nsamples(ps_flt)), ps_flt)
sum(taxa_sums(ps_flt))
length(taxa_sums(ps_flt))

existingDirCheck("../figures/thibault")
sample_data(ps)$gg_group <- factor(sample_data(ps)$gg_group2, levels = c("50:water","500:water","50:dss","500:dss")) # Put gg_group2 as factor
sample_data(ps)$timepoint <- factor(sample_data(ps)$timepoint, levels = c("0","35","49","54","final")) # Put timepoint as factor
sample_data(ps)$treatment <- factor(sample_data(ps)$treatment, levels = c("water","dss")) # Put timepoint as factor
p = alphaDiversityTimeSeries2(ps, "../figures/thibault/", time = "timepoint", group = "gg_group2")
p+scale_fill_manual(values = c("blue","red","blue","red"))+
  scale_color_manual(values = c("blue","red","blue","red"))+
  scale_pattern_manual(values = c("circle","stripe"))+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_blank(),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
    # legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    # legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
  
  
  
#Samuel alpha diversity
customColors = list('blue','red','blue','red')
pairs <- list(list("Wt:Vehicle","Wt:Putrescine"), list("IL-22ra1-/-:Vehicle","IL-22ra1-/-:Putrescine"), list("Wt:Vehicle","IL-22ra1-/-:Vehicle"), list("Wt:Putrescine","IL-22ra1-/-:Putrescine"))
alphaDiversityGgGroup2(ps_samuel, path = "~/Documents/CHUM_git/figures/samuel/new_alpha_diversity/", gg_group = "gg_group", customColors = customColors)

