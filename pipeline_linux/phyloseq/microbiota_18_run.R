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
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/alpha_diversity_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/beta_diversity_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/correlation_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/relab_analysis_graphs_and_stats.R")

#function to add SEM (1.96 for 95% confidence interval)
mean_cl_normal <- function(x, mult = 1.96) { #mult is 1.96 for a 95% confidence interval
  # Calculate the mean of the input vector x
  mean_val <- mean(x, na.rm = TRUE)
  
  # Calculate the standard error of the mean
  se_val <- sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
  
  # Return a data frame with the mean (y), and the lower (ymin) and upper (ymax) bounds
  data.frame(y = mean_val, ymin = mean_val - mult * se_val, ymax = mean_val + mult * se_val)
}

#Function checking if a dir exists and creating it otherwise
existingDirCheck <- function(path){
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Directory created: ", path)
  } else {
    message("Directory already exists: ", path)
  }
  
}

#Function to remove "/" and "-" characters from a string + lowercase
clean_string <- function(input_string) {
  result <- tolower(gsub("[/-]", "", input_string))
  return(result)
}

# Transform otu_table of a phyloseq object based on the chosen transformation
transformCounts <- function(ps, transformation = "rel_ab", log_base = 10) {
  if (transformation == "rel_ab") {
    
    # Check orientation of OTU table and transpose if necessary
    if (isFALSE(taxa_are_rows(otu_table(ps)))) {
      otu_table(ps) <- t(otu_table(ps))
      message("OTU table was transposed to have ASVs as rows.")
    } else {
      message("ASVs already as rows, no need to transpose OTU table.")
    }
    
    # Apply prop.table for relative abundance calculation and multiply by 100 for percentage
    otu_matrix <- apply(otu_table(ps), 2, prop.table) * 100
    
    # Ensure that the result is an otu_table object
    otu_table(ps) <- otu_table(otu_matrix, taxa_are_rows = TRUE)
    
  } else if (transformation == "log") {
    
    # Check if there are any zero counts to avoid log(0) issues
    if (any(otu_table(ps) == 0)) {
      message("Warning: Zero values detected in the OTU table. Adding a small constant to avoid log(0).")
      otu_table(ps) <- otu_table(ps) + 1e-6
    }
    
    # Apply log transformation (base 10 by default)
    otu_matrix <- log(otu_table(ps), base = log_base)
    
    # Ensure that the result is an otu_table object
    otu_table(ps) <- otu_table(otu_matrix, taxa_are_rows = TRUE)
  } else {
    stop("Transformation method not recognized. Use 'rel_ab' for relative abundance or 'log' for log transformation.")
  }
  
  return(ps) # Return the transformed phyloseq object
}

#for microbiota 18
#set working directory
setwd("~/Documents/CHUM_git/Microbiota_18/")
asv_table <- as.data.frame(fread("asv_table/asv_table_m1.csv", sep = ";"))

# Set the first column as row names and remove it from the data frame
rownames(asv_table) <- asv_table[,1]  # Use the first column as row names
asv_table <- asv_table[,-1]  # Drop the first column

{
#loading metadata of interest
metadata <- read.csv("metadata/metadata.csv", sep = ";")

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

#load taxonomical assignments
taxa <- as.matrix(fread("taxonomy/taxa_annotation_m1.csv", sep = ";"))

# Set the first column as row names and remove it from the data frame
rownames(taxa) <- taxa[,1]  # Use the first column as row names
taxa <- taxa[,-1]  # Drop the first column

{
# Load phylogenetic tree if possible
tree <- read.tree("~/Documents/CHUM_git/figures/samuel/beta_diversity/phylo_tree/phylogenetic_tree.newick")
}
  
#creating phyloseq object
ps <- phyloseq(otu_table(asv_table, taxa_are_rows = FALSE),
               tax_table(taxa)) #sample_data(metadata)


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
  write.tree(tree, file = "~/Documents/CHUM_git/figures/samuel/beta_diversity/phylo_tree/phylogenetic_tree.newick")
  
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
nrow(tax_table(ps))
sum(is.na(tax_table(ps)[,7]))
is.na(tax_table(ps)[,7])

# Function filtering out ASVs for which they were in total less than a threshold count
ps <- prune_taxa(taxa_sums(ps) > 10, ps)
sum(taxa_sums(ps))
length(taxa_sums(ps))

# Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
ps <- prune_taxa(colSums(otu_table(ps) > 0) >= (0.05 * nsamples(ps)), ps)
sum(taxa_sums(ps))
length(taxa_sums(ps))


