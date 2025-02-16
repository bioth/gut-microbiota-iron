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
  library(pairwiseAdonis)
  
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

# Load phylogenetic tree if possible
tree <- read.tree("~/Documents/CHUM_git/Microbiota_18/taxonomy/phylogenetic_tree.newick")
  
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
ps <- merge_phyloseq(ps, phy_tree(tree))

sum(taxa_sums(ps)) #total number of reads
length(taxa_sums(ps)) #total number of ASVs
nrow(tax_table(ps))-sum(is.na(tax_table(ps)[,7])) #how many species detected

# Function filtering out ASVs for which they were in total less than a threshold count
ps_flt <- prune_taxa(taxa_sums(ps) > 10, ps)
sum(taxa_sums(ps_flt))
length(taxa_sums(ps_flt))

# Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
ps_flt <- prune_taxa(colSums(otu_table(ps_flt) > 0) >= (0.05 * nsamples(ps_flt)), ps_flt)
sum(taxa_sums(ps_flt))
length(taxa_sums(ps_flt))

# Put as factors variables that are going to be used
sample_data(ps)$gg_group2 <- factor(sample_data(ps)$gg_group2, levels = c("50:water", "50:dss", "500:water", "500:dss")) # Put gg_group2 as factor
sample_data(ps)$timepoint <- factor(sample_data(ps)$timepoint, levels = c("0","35","49","54","final")) # Put timepoint as factor
sample_data(ps)$treatment <- factor(sample_data(ps)$treatment, levels = c("water","dss")) # Put treatment as factor
sample_data(ps)$diet <- factor(sample_data(ps)$diet, levels = c("50","500")) # Put diet as factor

# Put as factors variables that are going to be used
sample_data(ps_flt)$gg_group2 <- factor(sample_data(ps_flt)$gg_group2, levels = c("50:water", "50:dss", "500:water", "500:dss")) # Put gg_group2 as factor
sample_data(ps_flt)$timepoint <- factor(sample_data(ps_flt)$timepoint, levels = c("0","35","49","54","final")) # Put timepoint as factor
sample_data(ps_flt)$treatment <- factor(sample_data(ps_flt)$treatment, levels = c("water","dss")) # Put treatment as factor
sample_data(ps_flt)$diet <- factor(sample_data(ps_flt)$diet, levels = c("50","500")) # Put diet as factor


# Alpha diveristy
existingDirCheck("../figures/thibault")
graphs = alphaDiversityTimeSeries2(ps, "../figures/thibault/", time = "timepoint", group = "gg_group2", writeData = TRUE)

# Chao1
graphs[[1]]+
  scale_fill_manual(values = c("blue","blue4","red","red4"),
                    labels = c("50 ppm control","50 ppm DSS","500 ppm control","500 ppm DSS"))+
  scale_pattern_manual(values = c("circle","stripe"))+
  guides(pattern = "none")+
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

# Shannon
graphs[[2]]+
  scale_fill_manual(values = c("blue","blue4","red","red4"),
                    labels = c("50 ppm control","50 ppm DSS","500 ppm control","500 ppm DSS"))+
  scale_pattern_manual(values = c("circle","stripe"))+
  guides(pattern = "none")+
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

# InvSimpson
graphs[[3]]+
  scale_fill_manual(values = c("blue","blue4","red","red4"),
                    labels = c("50 ppm control","50 ppm DSS","500 ppm control","500 ppm DSS"))+
  scale_pattern_manual(values = c("circle","stripe"))+
  guides(pattern = "none")+
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

# Stats
alpha_d <- read.xlsx("../figures/thibault/alpha_diversity/alpha_diversity_data.xlsx")
alpha_d$diet <- factor(alpha_d$diet, levels = c("50", "500"))
alpha_d$treatment <- factor(alpha_d$treatment, levels = c("water", "dss"))
alpha_d$timepoint <- factor(alpha_d$timepoint, levels = c("0", "35", "49", "54", "final"))
measures <- c("Chao1", "Shannon", "InvSimpson")

sink("../figures/thibault/alpha_diversity/alpha_diversity_stats.txt")
for(measure in measures){
  
  print(paste("Tests for measure", measure))
  for(timepoint in levels(alpha_d$timepoint)){
    # hist(alpha_d[alpha_d$timepoint == "final","Chao1"], breaks = 20, main = "Histogram of Diversity Index")
    print(paste("At", timepoint))
    print(shapiro.test(alpha_d[alpha_d$timepoint == timepoint,measure])) # Nornality test
    print(bartlett.test(as.formula(paste(measure,"~ interaction(diet, treatment)")), data = alpha_d[alpha_d$timepoint == timepoint,])) # Homoscedasticity test (if data is normally distributed)
    anova_model <- aov(as.formula(paste(measure, "~ diet * treatment")), data = alpha_d[alpha_d$timepoint == timepoint,])
    print(summary(anova_model))
    print(TukeyHSD(anova_model))
    
  }
  
}
sink()


# Beta diversity (for only diets timepoints)
# Bray curtis
betaDiversityTimepoint2Factors(ps, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="bray",
                               transform = "rel_ab", customColors = c("blue","red"),
                               font = "Arial", path = "../figures/thibault/beta_diversity/diet_only/")
# Bray curtis filtered
betaDiversityTimepoint2Factors(ps_flt, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="bray",
                               transform = "rel_ab", customColors = c("blue","red"),
                               font = "Arial", path = "../figures/thibault/beta_diversity/filtered/")

# Weighted Unifrac
betaDiversityTimepoint2Factors(ps, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="wunifrac",
                               customColors = c("blue","red"),
                               font = "Arial", path = "../figures/thibault/beta_diversity/diet_only/")

# Beta diversity (for diet + treatment)
# Bray curtis
betaDiversityTimepoint2Factors(ps, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "gg_group2", distMethod ="bray",
                               customColors = c("blue","blue4","red","red4"),
                               font = "Arial", path = "../figures/thibault/beta_diversity/dss/")

# Weighted Unifrac
betaDiversityTimepoint2Factors(ps, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "gg_group2", distMethod ="wunifrac",
                               customColors = c("blue","blue4","red","red4"),
                               font = "Arial", path = "../figures/thibault/beta_diversity/dss/")


ps_dss <- prune_samples(sample_data(ps)$treatment == "dss", ps)

# Beta diversity (for diet + treatment), but only dss ones
# Bray curtis
betaDiversityTimepoint2Factors(ps_dss, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "gg_group2", distMethod ="bray",
                               customColors = c("blue","blue4","red","red4"),
                               font = "Arial", path = "../figures/thibault/beta_diversity/dss_only/")

# Weighted Unifrac
betaDiversityTimepoint2Factors(ps_dss, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "gg_group2", distMethod ="wunifrac",
                               customColors = c("blue","blue4","red","red4"),
                               font = "Arial", path = "../figures/thibault/beta_diversity/dss_only/")





# Idk
# Weighted Unifrac
betaDiversityTimepoint2Factors(ps_flt, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "gg_group2", distMethod ="wunifrac",
                               customColors = c("blue","blue4","red","red4"),
                               font = "Arial", path = "../figures/thibault/beta_diversity/desperate/")







