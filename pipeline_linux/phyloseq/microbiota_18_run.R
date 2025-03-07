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
  library(caret)
  
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
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/plot_microbiota_extension.R")

#for microbiota 18
#set working directory
setwd("~/Documents/CHUM_git/Microbiota_18/")
asv_table <- as.data.frame(fread("asv_table/asv_table_server.csv", sep = ";"))
asv_table <- as.data.frame(fread("from_server/to_transfer/asv_table/asv_table_m1.csv", sep = ";"))

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
colnames(samples) <- "sample_id"
samples$id <- substring(samples$sample_id, 1, 5)
samples$timepoint <- substring(samples$sample_id, 8, nchar(samples$sample_id)) 

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
rownames(metadata) <- metadata$sample_id
}

#load taxonomical assignments
taxa <- as.matrix(fread("taxonomy/taxa_annotation_server.csv", sep = ";"))
taxa <- as.matrix(fread("from_server/to_transfer/taxonomy/taxa_annotation_m1.csv", sep = ";"))

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
length(taxa_sums(ps))
{
  # Check if there are ASVs that are NAs at the phylum level
  length(tax_table(ps)[is.na(tax_table(ps)[,"Phylum"])])
  
  # Remove them as they probably correspond to chimeric sequences
  ps <- prune_taxa(!is.na(as.data.frame(tax_table(ps))$Phylum), ps)
  length(taxa_sums(ps))
  length(tax_table(ps)[is.na(tax_table(ps)[,"Phylum"])])
  
}

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

# Create four ps objects, one for diets (first 3 timepoints), one for dss (last 2 timepoints)
ps_diet <- prune_samples(sample_data(ps)$timepoint %in% c("0", "35", "49"), ps)
ps_dss_alpha <- prune_samples(sample_data(ps)$timepoint %in% c("49", "54", "final"), ps)
ps_dss_relab_flt <- prune_samples(sample_data(ps_flt)$timepoint %in% c("49", "54", "final"), ps_flt)
ps_dss <- prune_samples(sample_data(ps)$timepoint %in% c("54", "final"), ps)
ps_flt_diet <- prune_samples(sample_data(ps)$timepoint %in% c("0", "35", "49"), ps_flt)
ps_flt_dss <- prune_samples(sample_data(ps)$timepoint %in% c("54", "final"), ps_flt)

{
ps_foo <- ps
rownames(otu_table(ps_foo))
sample_data(ps_foo)[sample_data(ps_foo)$timepoint == 0,c("id","diet")]
sample_data(ps_foo)$diet[sample_data(ps_foo)$timepoint == 0 & sample_data(ps_foo)$id %in% c("33111","33112","33114")] <- c("500","500","500")
sample_data(ps_foo)$diet[sample_data(ps_foo)$timepoint == 0 & sample_data(ps_foo)$id %in% c("33111","33112","33114")] 
sample_data(ps_foo)$diet[sample_data(ps_foo)$timepoint == 0 & sample_data(ps_foo)$id %in% c("33105","10946","10966")] <- c("50","50","50")
sample_data(ps_foo)$diet[sample_data(ps_foo)$timepoint == 0 & sample_data(ps_foo)$id %in% c("33105","10946","10966")] 
ps_foo <- prune_samples(sample_data(ps_foo)$timepoint == "0", ps_foo)
}

# Alpha diveristy for diet only
existingDirCheck("../figures/thibault/diet/")
graphs = alphaDiversityTimeSeries2(ps_diet, "../figures/thibault/diet/", time = "timepoint", group = "diet", writeData = TRUE)

{
  ps_sub <- prune_samples(sample_data(ps_dss_alpha)$timepoint %in% c("49","54"), ps_dss_alpha)
  graphs = alphaDiversityTimeSeries2(ps_sub, "../figures/thibault/icm/", time = "timepoint", group = "gg_group2", writeData = FALSE)
  
  # Chao1
  p = graphs[[1]]+
    scale_fill_manual(values = c("blue","darkblue","red","darkred"),
                      labels = c("50 ppm control","50 ppm DSS","500 ppm control","500 ppm DSS"))+
    scale_pattern_manual(values = c("circle","stripe","circle","stripe"))+
    scale_x_discrete(labels = c("DSS day 0","DSS day 5", "End of\nrecovery"),
                     expand = c(0, 0.5))+
    labs(y = "Species Richness")+
    guides(pattern = "none")+
    theme_minimal()+
    theme(
      plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
      axis.title.x = element_blank(),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
      axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, face = "bold"),  # Adjust x-axis tick label font size
      axis.text.y = element_text(size = 12, face = "bold"),  # Adjust y-axis tick label font size
      axis.title.y = element_text(size = 14, face = "bold"),
      # legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
      # legend.text = element_text(size = 12),  # Adjust legend font size
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
  ggsave(plot = p, filename = "../figures/thibault/icm/chao1.png", width = 5, height = 4, dpi = 300, bg = "white")
  
  
  ps_sub <- prune_samples(sample_data(ps_dss_alpha)$timepoint %in% c("final"), ps_dss_alpha)
  graphs = alphaDiversityTimeSeries2(ps_sub, "../figures/thibault/icm/", time = "timepoint", group = "gg_group2", writeData = FALSE)
  
  # Chao1
  p = graphs[[1]]+
    scale_fill_manual(values = c("blue","darkblue","red","darkred"),
                      labels = c("50 ppm control","50 ppm DSS","500 ppm control","500 ppm DSS"))+
    scale_pattern_manual(values = c("circle","stripe","circle","stripe"))+
    scale_x_discrete(labels = c("End of\nrecovery"),
                     expand = c(0, 0.5))+
    labs(y = "Species Richness")+
    guides(pattern = "none")+
    theme_minimal()+
    theme(
      plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
      axis.title.x = element_blank(),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
      axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, face = "bold"),  # Adjust x-axis tick label font size
      axis.text.y = element_text(size = 12, face = "bold"),  # Adjust y-axis tick label font size
      axis.title.y = element_text(size = 14, face = "bold"),
      # legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
      # legend.text = element_text(size = 12),  # Adjust legend font size
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
  ggsave(plot = p, filename = "../figures/thibault/icm/chao1_2.png", width = 3.5, height = 3.5, dpi = 300, bg = "white")
  
  # Inverse simpson
  p = graphs[[3]]+
    scale_fill_manual(values = c("blue","darkblue","red","darkred"),
                      labels = c("50 ppm control","50 ppm DSS","500 ppm control","500 ppm DSS"))+
    scale_pattern_manual(values = c("circle","stripe","circle","stripe"))+
    scale_x_discrete(labels = c("End of\nrecovery"),
                     expand = c(0, 0.5))+
    labs(y = "Inverse Simpson")+
    guides(pattern = "none")+
    theme_minimal()+
    theme(
      plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
      axis.title.x = element_blank(),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
      axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, face = "bold"),  # Adjust x-axis tick label font size
      axis.text.y = element_text(size = 12, face = "bold"),  # Adjust y-axis tick label font size
      axis.title.y = element_text(size = 14, face = "bold"),
      # legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
      # legend.text = element_text(size = 12),  # Adjust legend font size
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
  ggsave(plot = p, filename = "../figures/thibault/icm/invsim.png", width = 3.5, height = 3.5, dpi = 300, bg = "white")
  
  
}


# Stats
alpha_d <- read.xlsx("../figures/thibault/diet/alpha_diversity/alpha_diversity_data.xlsx")
alpha_d$diet <- factor(alpha_d$diet, levels = c("50", "500"))
alpha_d$timepoint <- factor(alpha_d$timepoint, levels = c("0", "35", "49"))

# Chao1
graphs[[1]]+
  scale_fill_manual(values = c("blue","red"),
                    labels = c("50 ppm","500 ppm"))+
  scale_x_discrete(labels = c("3 weeks\n(Weaning)","5 weeks", "7 weeks"))+
  labs(y = "Chao1 Index")+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_blank(),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, face = "bold"),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12, face = "bold"),  # Adjust y-axis tick label font size
    axis.title.y = element_text(size = 14, face = "bold"),
    # legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    # legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars

print(shapiro.test(alpha_d[alpha_d$timepoint == "0","Chao1"])) # Nornality test
print(leveneTest(Chao1~diet, data = alpha_d[alpha_d$timepoint == "0",]))
wilcox.test(alpha_d$Chao1[alpha_d$timepoint == "0" & alpha_d$diet == "50"], alpha_d$Chao1[alpha_d$timepoint == "0" & alpha_d$diet == "500"])
print(bartlett.test(as.formula(paste(measure,"~ diet")), data = alpha_d[alpha_d$timepoint == timepoint,])) # Homoscedasticity test (if data is normally distributed)
# anova_model <- aov(as.formula(paste(measure, "~ diet")), data = alpha_d[alpha_d$timepoint == timepoint,])
print(t.test(as.formula(paste(measure,"~ diet")), data = alpha_d[alpha_d$timepoint == timepoint,]))

# Shannon
graphs[[2]]+
  scale_fill_manual(values = c("blue","red"),
                    labels = c("50 ppm","500 ppm"))+
  scale_x_discrete(labels = c("3 weeks\n(Weaning)","5 weeks", "7 weeks"))+
  labs(y = "Shannon Index")+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_blank(),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, face = "bold"),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12, face = "bold"),  # Adjust y-axis tick label font size
    axis.title.y = element_text(size = 14, face = "bold"),
    # legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    # legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars

print(shapiro.test(alpha_d[alpha_d$timepoint == "49","Shannon"])) # Nornality test
print(leveneTest(Shannon~diet, data = alpha_d[alpha_d$timepoint == "49",]))
wilcox.test(alpha_d$Shannon[alpha_d$timepoint == "49" & alpha_d$diet == "50"], alpha_d$Shannon[alpha_d$timepoint == "49" & alpha_d$diet == "500"])
print(bartlett.test(as.formula(paste(measure,"~ diet")), data = alpha_d[alpha_d$timepoint == timepoint,])) # Homoscedasticity test (if data is normally distributed)
# anova_model <- aov(as.formula(paste(measure, "~ diet")), data = alpha_d[alpha_d$timepoint == timepoint,])
print(t.test(as.formula(paste(measure,"~ diet")), data = alpha_d[alpha_d$timepoint == timepoint,]))

# InvSimpson
graphs[[3]]+
  scale_fill_manual(values = c("blue","red"),
                    labels = c("50 ppm","500 ppm"))+
  scale_x_discrete(labels = c("3 weeks\n(Weaning)","5 weeks", "7 weeks"))+
  labs(y = "Inverse Simpson")+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_blank(),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, face = "bold"),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12, face = "bold"),  # Adjust y-axis tick label font size
    axis.title.y = element_text(size = 14, face = "bold"),
    # legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    # legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars

print(shapiro.test(alpha_d[alpha_d$timepoint == "0","InvSimpson"])) # Nornality test
print(leveneTest(InvSimpson~diet, data = alpha_d[alpha_d$timepoint == "0",]))
wilcox.test(alpha_d$InvSimpson[alpha_d$timepoint == "0" & alpha_d$diet == "50"], alpha_d$InvSimpson[alpha_d$timepoint == "0" & alpha_d$diet == "500"])
# print(bartlett.test(as.formula(paste(measure,"~ diet")), data = alpha_d[alpha_d$timepoint == timepoint,])) # Homoscedasticity test (if data is normally distributed)
# anova_model <- aov(as.formula(paste(measure, "~ diet")), data = alpha_d[alpha_d$timepoint == timepoint,])
# print(t.test(as.formula(paste(measure,"~ diet")), data = alpha_d[alpha_d$timepoint == timepoint,]))

print(shapiro.test(alpha_d[alpha_d$timepoint == "35","InvSimpson"])) # Nornality test
print(leveneTest(InvSimpson~diet, data = alpha_d[alpha_d$timepoint == "35",]))
wilcox.test(alpha_d$InvSimpson[alpha_d$timepoint == "35" & alpha_d$diet == "50"], alpha_d$InvSimpson[alpha_d$timepoint == "35" & alpha_d$diet == "500"])

print(shapiro.test(alpha_d[alpha_d$timepoint == "49","InvSimpson"])) # Nornality test
print(leveneTest(InvSimpson~diet, data = alpha_d[alpha_d$timepoint == "49",]))
wilcox.test(alpha_d$InvSimpson[alpha_d$timepoint == "49" & alpha_d$diet == "50"], alpha_d$InvSimpson[alpha_d$timepoint == "49" & alpha_d$diet == "500"])

# Stats
alpha_d <- read.xlsx("../figures/thibault/diet/alpha_diversity/alpha_diversity_data.xlsx")
alpha_d$diet <- factor(alpha_d$diet, levels = c("50", "500"))
alpha_d$timepoint <- factor(alpha_d$timepoint, levels = c("0", "35", "49"))
measures <- c("Chao1", "Shannon", "InvSimpson")

sink("../figures/thibault/diet/alpha_diversity/alpha_diversity_stats.txt")
for(measure in measures){
  
  print(paste("Tests for measure", measure))
  for(timepoint in levels(alpha_d$timepoint)){
    # hist(alpha_d[alpha_d$timepoint == "final","Chao1"], breaks = 20, main = "Histogram of Diversity Index")
    print(paste("At", timepoint))
    print(shapiro.test(alpha_d[alpha_d$timepoint == timepoint,measure])) # Nornality test
    print(bartlett.test(as.formula(paste(measure,"~ diet")), data = alpha_d[alpha_d$timepoint == timepoint,])) # Homoscedasticity test (if data is normally distributed)
    # anova_model <- aov(as.formula(paste(measure, "~ diet")), data = alpha_d[alpha_d$timepoint == timepoint,])
    print(t.test(as.formula(paste(measure,"~ diet")), data = alpha_d[alpha_d$timepoint == timepoint,]))
    # print(summary(anova_model))
    # print(TukeyHSD(anova_model))
    
  }
  
}
sink()



# Alpha diveristy for dss_diet only
existingDirCheck("../figures/thibault/diet_dss/")
graphs = alphaDiversityTimeSeries2(ps_dss_alpha, "../figures/thibault/diet_dss/", time = "timepoint", group = "gg_group2", writeData = TRUE)


# Chao1
p = graphs[[1]]+
  scale_fill_manual(values = c("blue","darkblue","red","darkred"),
                    labels = c("50 ppm control","50 ppm DSS","500 ppm control","500 ppm DSS"))+
  scale_pattern_manual(values = c("circle","stripe","circle","stripe"))+
  scale_x_discrete(labels = c("DSS day 0","DSS day 5", "End of\nrecovery"),
                   expand = c(0, 0.5))+
  labs(y = "Species Richness")+
  guides(pattern = "none")+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_blank(),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, face = "bold"),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12, face = "bold"),  # Adjust y-axis tick label font size
    axis.title.y = element_text(size = 14, face = "bold"),
    # legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    # legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
ggsave(plot = p, filename = "../figures/thibault/icm/chao1.png", width = 6, height = 4, dpi = 300, bg = "white")
nrmTest <- function(x,y,log=FALSE){
  
  if(log){
    print(shapiro.test(log(alpha_d[alpha_d$timepoint == x & alpha_d$gg_group2 == levels(alpha_d$gg_group2)[1],y]))) 
    print(shapiro.test(log(alpha_d[alpha_d$timepoint == x & alpha_d$gg_group2 == levels(alpha_d$gg_group2)[2],y]))) 
    print(shapiro.test(log(alpha_d[alpha_d$timepoint == x & alpha_d$gg_group2 == levels(alpha_d$gg_group2)[3],y]))) 
    print(shapiro.test(log(alpha_d[alpha_d$timepoint == x & alpha_d$gg_group2 == levels(alpha_d$gg_group2)[4],y]))) 
  }else{
    print(shapiro.test(alpha_d[alpha_d$timepoint == x & alpha_d$gg_group2 == levels(alpha_d$gg_group2)[1],y])) 
    print(shapiro.test(alpha_d[alpha_d$timepoint == x & alpha_d$gg_group2 == levels(alpha_d$gg_group2)[2],y])) 
    print(shapiro.test(alpha_d[alpha_d$timepoint == x & alpha_d$gg_group2 == levels(alpha_d$gg_group2)[3],y])) 
    print(shapiro.test(alpha_d[alpha_d$timepoint == x & alpha_d$gg_group2 == levels(alpha_d$gg_group2)[4],y])) 
  }
}

nrmTest(49, "Chao1")
nrmTest(54, "Chao1")
nrmTest("final", "Chao1")

# Shannon
graphs[[2]]+
  scale_fill_manual(values = c("blue","blueviolet","red","darkorange"),
                    labels = c("50 ppm control","50 ppm DSS","500 ppm control","500 ppm DSS"))+
  scale_pattern_manual(values = c("circle","stripe"))+
  scale_x_discrete(labels = c("DSS day 0","DSS day 5", "Final timepoint\n(8 weeks recovery)"))+
  labs(y = "Shannon Index")+
  guides(pattern = "none")+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_blank(),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, face = "bold"),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12, face = "bold"),  # Adjust y-axis tick label font size
    axis.title.y = element_text(size = 14, face = "bold"),
    # legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    # legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars

print(shapiro.test(alpha_d[alpha_d$timepoint == 49,"Shannon"])) # Nornality test
nrmTest(49, "Shannon", log = TRUE)
print(leveneTest(Shannon~diet*treatment, data = alpha_d[alpha_d$timepoint == "49",]))
hist(alpha_d[alpha_d$timepoint == "49","Shannon"], breaks = 10, main = "Histogram of Diversity Index")
model <- glm(Shannon ~ diet * treatment, family = Gamma(link = "log"), data = alpha_d[alpha_d$timepoint == 49,])
summary(model)
# print(bartlett.test(as.formula(paste(measure,"~ interaction(diet, treatment)")), data = alpha_d[alpha_d$timepoint == timepoint,])) # Homoscedasticity test (if data is normally distributed)
# anova_model <- aov(as.formula(paste(measure, "~ diet * treatment")), data = alpha_d[alpha_d$timepoint == timepoint,])
# print(summary(anova_model))
# print(TukeyHSD(anova_model))

nrmTest(54, "Shannon")

nrmTest("final", "Shannon")
print(leveneTest(Shannon~diet*treatment, data = alpha_d[alpha_d$timepoint == "final",]))
hist(alpha_d[alpha_d$timepoint == "final","Shannon"], breaks = 10, main = "Histogram of Diversity Index")
model <- glm(Shannon ~ diet * treatment, family = Gamma(link = "log"), data = alpha_d[alpha_d$timepoint == "final",])
summary(model)


# InvSimpson
p = graphs[[3]]+
  scale_fill_manual(values = c("blue","darkblue","red","darkred"),
                    labels = c("50 ppm control","50 ppm DSS","500 ppm control","500 ppm DSS"))+
  scale_pattern_manual(values = c("circle","stripe"))+
  scale_x_discrete(labels = c("DSS day 0","DSS day 5", "End of\nrecovery"),
                   expand = c(0, 0.5))+
  labs(y = "Inverse Simpson")+
  guides(pattern = "none")+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_blank(),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, face = "bold"),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12, face = "bold"),  # Adjust y-axis tick label font size
    axis.title.y = element_text(size = 14, face = "bold"),
    # legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    # legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
ggsave(plot = p, filename = "../figures/thibault/icm/invSim.png", width = 6, height = 4, dpi = 300, bg = "white")


nrmTest("49", "InvSimpson", log = TRUE)
print(leveneTest(InvSimpson~diet*treatment, data = alpha_d[alpha_d$timepoint == "49",]))
hist(alpha_d[alpha_d$timepoint == "49","InvSimpson"], breaks = 10, main = "Histogram of Diversity Index")
model <- glm(InvSimpson ~ diet * treatment, family = Gamma(link = "log"), data = alpha_d[alpha_d$timepoint == "49",])
summary(model)

nrmTest("54", "InvSimpson", log = TRUE)
alpha_d$logSimpson <- log(alpha_d$InvSimpson)
print(bartlett.test(as.formula("logSimpson ~ interaction(diet, treatment)"), data = alpha_d[alpha_d$timepoint == 54,]))
print(leveneTest(InvSimpson~diet*treatment, data = alpha_d[alpha_d$timepoint == "54",]))
hist(log(alpha_d[alpha_d$timepoint == "54","InvSimpson"]), breaks = 10, main = "Histogram of Diversity Index")
anova_model <- aov(as.formula("logSimpson~ diet * treatment"), data = alpha_d[alpha_d$timepoint == 54,])
print(summary(anova_model))
print(TukeyHSD(anova_model))

nrmTest("final", "InvSimpson", log = TRUE)
print(bartlett.test(as.formula("logSimpson ~ interaction(diet, treatment)"), data = alpha_d[alpha_d$timepoint == "final",]))
anova_model <- aov(as.formula("logSimpson~ diet * treatment"), data = alpha_d[alpha_d$timepoint == "final",])
print(summary(anova_model))
print(TukeyHSD(anova_model))

# Stats
alpha_d <- read.xlsx("../figures/thibault/diet_dss/alpha_diversity/alpha_diversity_data.xlsx")
alpha_d$diet <- factor(alpha_d$diet, levels = c("50", "500"))
alpha_d$treatment <- factor(alpha_d$treatment, levels = c("water", "dss"))
alpha_d$timepoint <- factor(alpha_d$timepoint, levels = c("49", "54", "final"))
alpha_d$gg_group2 <- factor(alpha_d$gg_group2, levels = c("50:water","50:dss","500:water","500:dss"))
measures <- c("Chao1", "Shannon", "InvSimpson")

sink("../figures/thibault/diet_dss/alpha_diversity/alpha_diversity_stats.txt")
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
betaDiversityTimepoint2Factors(ps_diet, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="bray",
                               transform = "rel_ab", customColors = c("blue","red"), displaySampleIDs = TRUE,
                               font = "Arial", path = "../figures/thibault/diet/beta_diversity/")
# Bray curtis filtered
betaDiversityTimepoint2Factors(ps_flt_diet, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="bray",
                               transform = "rel_ab", customColors = c("blue","red"),
                               font = "Arial", path = "../figures/thibault/diet/beta_diversity/filtered/")

# Weighted Unifrac
betaDiversityTimepoint2Factors(ps_diet, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="wunifrac", displaySampleIDs = TRUE,
                               customColors = c("blue","red"),
                               font = "Arial", path = "../figures/thibault/diet/beta_diversity/")

# Weighted Unifrac filtered
betaDiversityTimepoint2Factors(ps_flt_diet, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="wunifrac",
                               customColors = c("blue","red"), dim = c(4.5,4.5),
                               font = "Arial", path = "../figures/thibault/icm/",
                               additionnalAes = 
                                 theme(
                                   plot.title = element_text(),
                                   panel.grid.major = element_blank(),  # Add major grid lines
                                   panel.grid.minor = element_blank(),  # Remove minor grid lines
                                   aspect.ratio = 1
                                 ))

# Weighted Unifrac (nice graph at t0 for ICM day DO NOT USE AFTER)
betaDiversityTimepoint2Factors(ps_foo, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="wunifrac", displaySampleIDs = FALSE,
                               customColors = c("blue","red"), dim = c(4.5,4.5),
                               font = "Arial", path = "../figures/thibault/icm/", 
                               additionnalAes = 
                                 theme(
                                   plot.title = element_text(),
                                   panel.grid.major = element_blank(),  # Add major grid lines
                                   panel.grid.minor = element_blank(),  # Remove minor grid lines
                                   aspect.ratio = 1
                                 ))



# Beta diversity (for diet + treatment)
# Bray curtis
betaDiversityTimepoint2Factors(ps_dss, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "gg_group2", distMethod ="bray",
                               customColors = c("blue","blue4","red","red4"), dim = c(5,5),
                               font = "Arial", path = "../figures/thibault/diet_dss/beta_diversity/",
                               additionnalAes = 
                                 theme(
                                   plot.title = element_text(),
                                   panel.grid.major = element_blank(),  # Add major grid lines
                                   panel.grid.minor = element_blank(),  # Remove minor grid lines
                                 ))

# Bray curtis filtered
betaDiversityTimepoint2Factors(ps_flt_dss, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "gg_group2", distMethod ="bray",
                               customColors = c("blue","blue4","red","red4"),
                               font = "Arial", path = "../figures/thibault/diet_dss/beta_diversity/filtered/")

# Weighted Unifrac filtered, at t54 and tFinal
betaDiversityTimepoint2Factors(ps_flt_dss, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "gg_group2", distMethod ="wunifrac", 
                               customColors = c("blue","blueviolet","red","darkorange"), dim = c(5,5),
                               font = "Arial", path = "../figures/thibault/diet_dss/beta_diversity/filtered/", 
                               additionnalAes = 
                               theme(
                                 plot.title = element_text(),
                                 panel.grid.major = element_blank(),  # Add major grid lines
                                 panel.grid.minor = element_blank(),  # Remove minor grid lines
                               ))

# dbRDA method at t54 and for DSS groups
ps_sub <- prune_samples(sample_data(ps_flt_dss)$timepoint == "54" & sample_data(ps_flt_dss)$treatment == "dss", ps_flt_dss)
sample_data(ps_sub)$diet
betaDiversityTimepointsGroupedDbRDA(ps_sub, sample_id = "full_id", varToCompare = "diet", formula = "diet",
                                    transform = "none", distMethod = "wunifrac", customColors = c("blue4","red4"),
                                    font = "Arial", path = "../figures/thibault/dss_only/beta_diversity/dbRDA_t54/")

# dbRDA method at last timepoint and for DSS groups
ps_sub <- prune_samples(sample_data(ps_flt_dss)$timepoint == "final" & sample_data(ps_flt_dss)$treatment == "dss", ps_flt_dss)
sample_data(ps_sub)$diet
print(length(taxa_sums(ps_sub)))
ps_sub <- prune_taxa(taxa_sums(ps_sub) > 10, ps_sub)
ps_sub <- prune_taxa(colSums(otu_table(ps_sub) > 0) >= (0.05 * nsamples(ps_sub)), ps_sub)
print(length(taxa_sums(ps_sub)))
betaDiversityTimepointsGroupedDbRDA(ps_sub, sample_id = "full_id", varToCompare = "diet", formula = "diet",
                                    transform = "none", distMethod = "wunifrac", customColors = c("darkblue","darkred"),
                                    font = "Arial", path = "../figures/thibault/dss_only/beta_diversity/dbRDA/",
                                    additionnalAes = 
                                      theme(
                                        plot.title = element_text(),
                                        panel.grid.major = element_blank(),  # Add major grid lines
                                        panel.grid.minor = element_blank(),  # Remove minor grid lines
                                        aspect.ratio = 1
                                      ))

# dbRDA method at t54 and last timepoint and for DSS groups
ps_sub <- prune_samples(sample_data(ps_flt_dss)$treatment == "dss", ps_flt_dss)
sample_data(ps_sub)$diet
sample_data(ps_sub)$gg_group <- factor(sample_data(ps_sub)$gg_group, levels = c("54:50:dss","54:500:dss","final:50:dss" ,"final:500:dss"))
betaDiversityTimepointsGroupedDbRDA(ps_sub, sample_id = "full_id", varToCompare = "gg_group", formula = "diet*timepoint",
                                    transform = "none", distMethod = "wunifrac", customColors = c("blue","red","blue4","red4"),
                                    font = "Arial", path = "../figures/thibault/dss_only/beta_diversity/dbRDA_2Last/")

# dbRDA method at last timepoint and for all groups
ps_sub <- prune_samples(sample_data(ps_flt_dss)$timepoint == "final", ps_flt_dss)

{
  # Check if there are ASVs that are NAs at the phylum level
  length(tax_table(ps_sub)[is.na(tax_table(ps_sub)[,2])])
  
  # Check if there are ASVs that are NAs at the order level
  length(tax_table(ps_sub)[is.na(tax_table(ps_sub)[,3])])
  
  # Remove them as they probably correspond to chimeric sequences
  ps <- prune_taxa(!is.na(as.data.frame(tax_table(ps))$Phylum), ps)
  ps_sub <- prune_taxa(!is.na(as.data.frame(tax_table(ps_sub))$Order), ps_sub)
  length(taxa_sums(ps_sub))
  length(tax_table(ps)[is.na(tax_table(ps)[,"Phylum"])])
  
}

sum(taxa_sums(ps_sub)) #total number of reads
length(taxa_sums(ps_sub)) #total number of ASVs
nrow(tax_table(ps_sub))-sum(is.na(tax_table(ps_sub)[,7])) #how many species detected

# Function filtering out ASVs for which they were in total less than a threshold count
ps_sub <- prune_taxa(taxa_sums(ps_sub) > 10, ps_sub)
sum(taxa_sums(ps_flt))
length(taxa_sums(ps_flt))

# Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
ps_sub <- prune_taxa(colSums(otu_table(ps_sub) > 0) >= (0.05 * nsamples(ps_sub)), ps_sub)
sum(taxa_sums(ps_sub))
length(taxa_sums(ps_sub))


sample_data(ps_sub)$gg_group2
betaDiversityTimepointsGroupedDbRDA(ps_sub, sample_id = "full_id", varToCompare = "gg_group2", formula = "diet*treatment",
                                    transform = "none", distMethod = "wunifrac", customColors = c("blue","red","blue4","red4"),
                                    font = "Arial", path = "../figures/thibault/dss_only/beta_diversity/dbRDA_Last/")

{
# Only final timepoint, across all groups
ps_final_dss <- prune_samples(sample_data(ps_flt_dss)$timepoint == "final", ps_flt_dss)
betaDiversityTimepointsGroupedDbRDA(ps_final_dss, sample_id = "full_id", varToCompare = "gg_group2", distMethod = "wunifrac", formula = "diet*treatment",
                                  transform = "none", customColors = c("blue","blueviolet","red","darkorange"), 
                                  font = "Arial", path = "../figures/thibault/diet_dss/beta_diversity/filtered/dbRDA_final/")

# Only t54 timepoint, across all groups
ps_54_dss <- prune_samples(sample_data(ps_flt_dss)$timepoint == "54", ps_flt_dss)
betaDiversityTimepointsGroupedDbRDA(ps_54_dss, sample_id = "full_id", varToCompare = "gg_group2", distMethod = "wunifrac", formula = "diet*treatment",
                                    transform = "none", customColors = c("blue","blueviolet","red","darkorange"), 
                                    font = "Arial", path = "../figures/thibault/diet_dss/beta_diversity/filtered/dbRDA_final/")



# Only DSS groups
ps_dss_only <- prune_samples(sample_data(ps_dss)$treatment == "dss", ps_dss)


# Beta diversity (for diet + treatment), but only dss ones
# Bray curtis
betaDiversityTimepoint2Factors(ps_dss_only, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "gg_group2", distMethod ="bray",
                               customColors = c("blue4","red4"),
                               font = "Arial", path = "../figures/thibault/diet_dss/beta_diversity/dss_only/")

# Weighted Unifrac
betaDiversityTimepoint2Factors(ps_dss_only, sample_id = "full_id", timeVariable = "timepoint",
                               varToCompare =  "gg_group2", distMethod ="wunifrac",
                               customColors = c("blue4","red4"),
                               font = "Arial", path = "../figures/thibault/diet_dss/beta_diversity/dss_only/")




# RDA method (constrained analysis), for both diet and treatment explanotory variables
betaDiversityTimepoint2FactorsRDA(ps_dss, sample_id = "full_id", timeVariable = "timepoint",
                                  varToCompare = "gg_group2", formula = "diet * treatment",
                                  transform = "rel_ab", customColors = c("blue","blue4","red","red4"),
                                  font = "Arial", path = "../figures/thibault/diet_dss/beta_diversity/RDA/")


ps_dss_only_flt <- prune_samples(sample_data(ps_flt_dss)$treatment == "dss", ps_flt_dss)
sample_data(ps_dss_only)$gg_group <- factor(sample_data(ps_dss_only)$gg_group, levels = c("54:50:dss", "54:500:dss", "final:50:dss", "final:500:dss"))
sample_data(ps_dss_only)$timepoint

sample_data(ps_dss_only_flt)$gg_group <- factor(sample_data(ps_dss_only_flt)$gg_group, levels = c("54:50:dss", "54:500:dss", "final:50:dss", "final:500:dss"))
sample_data(ps_dss_only_flt)$timepoint

# RDA method (constrained analysis), for both only treatment explanotory variables
betaDiversityTimepointsGroupedRDA(ps_dss_only, sample_id = "full_id", varToCompare = "gg_group", formula = "diet*timepoint",
                                  transform = "rel_ab", customColors = c("blue", "red", "blue4","red4"),
                                  font = "Arial", path = "../figures/thibault/dss_only/beta_diversity/RDA/")

# dbRDA method at last 2 timepoints and for DSS groups
betaDiversityTimepointsGroupedDbRDA(ps_dss_only_flt, sample_id = "full_id", varToCompare = "gg_group", formula = "diet*timepoint",
                                  transform = "none", customColors = c("blue", "red", "blue4","red4"),
                                  font = "Arial", path = "../figures/thibault/dss_only/beta_diversity/test/")



# Only final timepoint
ps_final_dss <- prune_samples(sample_data(ps_dss_only_flt)$timepoint == "final", ps_dss_only_flt)
betaDiversityTimepointsGroupedRDA(ps_final_dss, sample_id = "full_id", varToCompare = "gg_group", formula = "diet",
                                  transform = "rel_ab", customColors = c("blue", "red", "blue4","red4"),
                                  font = "Arial", path = "../figures/thibault/dss_only/beta_diversity/RDA_final/")




# dbRDA method at last timepoint and for DSS groups
ps_sub <- prune_samples(sample_data(ps_dss_only_flt)$timepoint == "final", ps_dss_only_flt)
betaDiversityTimepointsGroupedDbRDA(ps_sub, sample_id = "full_id", varToCompare = "gg_group", formula = "diet",
                                    transform = "none", customColors = c("blue4","red4"),
                                    font = "Arial", path = "../figures/thibault/dss_only/beta_diversity/test2/")







# dbRDA weighted unifrac at last timepoint and for all groups
ps_sub <- prune_samples(sample_data(ps_dss_flt)$timepoint == "final", ps_dss)
betaDiversityTimepointsGroupedDbRDA(ps_sub, sample_id = "full_id", varToCompare = "gg_group2", formula = "diet * treatment",
                                    transform = "none", customColors = c("blue", "red","blue4","red4"),
                                    font = "Arial", path = "../figures/thibault/dss_only/beta_diversity/test3/")

}







# Relative abundance analysis: finding differential abundant bugs at the species level, for DSS groups ONLY
#Path where to save graphs
existingDirCheck("~/Documents/CHUM_git/figures/thibault_new")
{
pathToSave <- "~/Documents/CHUM_git/figures/thibault_new/relative_abundance_by_timepoint/"
existingDirCheck(pathToSave)

#customColors for graph display
customColors = c("darkblue","darkred")

#Iterate through timepoints
for(timePoint in levels(sample_data(ps_dss_relab_flt)$timepoint)){
  
  #New path created for each week
  newPath <- paste(pathToSave, "timepoint_", timePoint, "/", sep = "")
  existingDirCheck(newPath)
  
  #Creating phyloseq objects for each timepoint
  ps_subset <- prune_samples(sample_data(ps_dss_relab_flt)$timepoint == timePoint & sample_data(ps_dss_relab_flt)$treatment == "dss", ps_dss_relab_flt)
  print(length(taxa_sums(ps_subset)))
  
  # Filtering
  # Function filtering out ASVs for which they were in total less than a threshold count
  ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
  # Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
  ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
  # print(sum(taxa_sums(ps_subset)))
  print(length(taxa_sums(ps_subset)))
  
  #Simple deseq object only accounting for the differences in diet
  deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ diet) 
  
  #Performing the deseq analysis
  deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
  
  print(resultsNames(deseq_subset))
  
  #For a given taxononical levels, creates graph for each timepoint, displaying which species were found to be differentially abundant
  relabSingleTimepoint(ps_subset, deseq_subset, measure = "log2fold", varToCompare = "diet",
                       timePoint = timePoint, taxa = "Species", threshold = 0.05, FDR = FALSE,
                       LDA = TRUE, customColors = customColors, path = newPath, displayPvalue = TRUE, additionnalAes = NULL)  
}

# At other taxonomic levels
taxonomicLevels <- c("Genus","Family","Order","Class","Phylum")

# Iterate through timepoints
for(timePoint in levels(sample_data(ps_dss_relab_flt)$timepoint)){
  
  # New path created for each week
  newPath <- paste(pathToSave, "timepoint_", timePoint, "/", sep = "")
  existingDirCheck(newPath)
  
  # Creating phyloseq objects for each timepoint
  ps_subset <- prune_samples(sample_data(ps_dss_relab_flt)$timepoint == timePoint & sample_data(ps_dss_relab_flt)$treatment == "dss", ps_dss_relab_flt)
  ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
  ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
  # Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)

  # print(length(taxa_sums(ps_subset)))
  # ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
  # print(length(taxa_sums(ps_subset)))
  # 
  # # Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
  # ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
  # print(length(taxa_sums(ps_subset)))
  
  for(txnLevel in taxonomicLevels){
    
    #Creates ps subset for taxonomical level of interest
    ps_taxa <- tax_glom(ps_subset, taxrank = txnLevel)
    colnames(sample_data(ps_taxa))[2] <- "sample_id"
    deseq_subset <- phyloseq_to_deseq2(ps_taxa, ~ diet)
    deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
    
    #For a given taxononical levels, creates graph for each timepoint, displaying which species were found to be differentially abundant
    relabSingleTimepoint(ps_taxa, deseq_subset, measure = "log2fold", "diet", timePoint = timePoint, taxa = txnLevel, threshold = 0.05, LDA = FALSE, customColors, newPath,
                         additionnalAes = theme(
                           axis.text.x = element_text(angle = 0, hjust = 0.5),
                           title = element_text(size = 12)),
                         dim = c(3.5,3.5), displayPvalue = TRUE) 
  }
  
  
  
}
}

# Relative abundance analysis: finding differential abundant bugs at the species level, all groups
#Path where to save graphs
{
  pathToSave <- "~/Documents/CHUM_git/figures/thibault_new/relative_abundance_by_timepoint_all_groups/"
  existingDirCheck(pathToSave)
  
  #customColors for graph display
  customColors = c("blue", "red", "darkblue","darkred")
  
  sample_data(ps_dss_relab_flt)$gg_group2 <- factor(sample_data(ps_dss_relab_flt)$gg_group2, levels = c("50:water", "500:water", "50:dss", "500:dss"))
  
  #Iterate through timepoints
  for(timePoint in levels(sample_data(ps_dss_relab_flt)$timepoint)){
    
    #New path created for each week
    newPath <- paste(pathToSave, "timepoint_", timePoint, "/", sep = "")
    existingDirCheck(newPath)
    
    #Creating phyloseq objects for each timepoint
    ps_subset <- prune_samples(sample_data(ps_dss_relab_flt)$timepoint == timePoint, ps_dss_relab_flt)
    print(length(taxa_sums(ps_subset)))
    
    # Filtering
    # Function filtering out ASVs for which they were in total less than a threshold count
    ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
    # Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
    ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
    # print(sum(taxa_sums(ps_subset)))
    print(length(taxa_sums(ps_subset)))
    
    #Simple deseq object only accounting for the differences in diet
    deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ treatment*diet) 
    
    #Performing the deseq analysis
    deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
    
    print(resultsNames(deseq_subset))
    
    #For a given taxononical levels, creates graph for each timepoint, displaying which species were found to be differentially abundant
    relabGroups(ps_subset, deseq_subset, measure = "log2fold", gg_group = "gg_group2", taxa = "Species", threshold = 0.01, FDR = FALSE,
                returnSigAsvs = FALSE, normalizeCounts = FALSE, customColors = customColors, 
                pairs = list(
                  list("50:water","500:water"), list("50:dss","500:dss"),
                  list("50:water","50:dss"), list("500:water","500:dss")),
                path = newPath, single_factor_design = FALSE,
                additionnalAes = NULL, dim = c(6,6), displayPvalue = FALSE)  
  }
  
  #customColors for graph display
  customColors = c("blue", "red", "darkblue","darkred")
  
  sample_data(ps_dss_relab_flt)$gg_group2 <- factor(sample_data(ps_dss_relab_flt)$gg_group2, levels = c("50:water", "500:water", "50:dss", "500:dss"))
  
  #At other taxonomic levels
  taxonomicLevels <- c("Genus","Family","Order","Class","Phylum")
  taxonomicLevels <- c("Phylum")
  
  #Iterate through timepoints
  for(timePoint in levels(sample_data(ps_dss_relab_flt)$timepoint)){
    
    #New path created for each week
    newPath <- paste(pathToSave, "timepoint_", timePoint, "/", sep = "")
    existingDirCheck(newPath)
    
    #Creating phyloseq objects for each timepoint
    ps_subset <- prune_samples(sample_data(ps_dss_relab_flt)$timepoint == timePoint, ps_dss_relab_flt)
    print(length(taxa_sums(ps_subset)))
    
    # # Filtering
    # # Function filtering out ASVs for which they were in total less than a threshold count
    # ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
    # # Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
    # ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
    # # print(sum(taxa_sums(ps_subset)))
    # print(length(taxa_sums(ps_subset)))
    
    for(txnLevel in taxonomicLevels){
      
      #Creates ps subset for taxonomical level of interest
      ps_taxa <- tax_glom(ps_subset, taxrank = txnLevel)
      colnames(sample_data(ps_taxa))[2] <- "sample_id"
      deseq_subset <- phyloseq_to_deseq2(ps_taxa, ~ treatment*diet)
      deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
      
      print(resultsNames(deseq_subset))
      
      #For a given taxononical levels, creates graph for each timepoint, displaying which species were found to be differentially abundant
      relabGroups(ps_taxa, deseq_subset, measure = "log2fold", gg_group = "gg_group2", taxa = txnLevel, threshold = 0.05, FDR = FALSE,
                  returnSigAsvs = FALSE, normalizeCounts = FALSE, customColors = customColors, 
                  pairs = list(
                    list("50:water","500:water"), list("50:dss","500:dss"),
                    list("50:water","50:dss"), list("500:water","500:dss")),
                  path = newPath, single_factor_design = FALSE,
                  additionnalAes = NULL, dim = c(6,6), displayPvalue = FALSE)  
    }
    
    
    
  }
}

# Relative abundance analysis: finding differential abundant bugs at the species level, for diet groups only
{
# Path where to save graphs
pathToSave <- "~/Documents/CHUM_git/figures/thibault_new/relative_abundance_diet/"
existingDirCheck(pathToSave)

#customColors for graph display
customColors = c("blue","red")

#Iterate through timepoints
for(timePoint in levels(sample_data(ps_flt_diet)$timepoint)){
  
  #New path created for each week
  newPath <- paste(pathToSave, "timepoint_", timePoint, "/", sep = "")
  existingDirCheck(newPath)
  
  #Creating phyloseq objects for each timepoint
  ps_subset <- prune_samples(sample_data(ps_flt_diet)$timepoint == timePoint, ps_flt_diet)
  # ps_sub <- prune_taxa(taxa_sums(ps_sub) > 10, ps_sub)
  colnames(sample_data(ps_subset))[2] <- "sample_id"
  # Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
  # ps_sub <- prune_taxa(colSums(otu_table(ps_sub) > 0) >= (0.05 * nsamples(ps_sub)), ps_sub)
  # print(length(taxa_sums(ps_subset)))
  # ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
  # print(length(taxa_sums(ps_subset)))
  # 
  # # Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
  # ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
  # print(length(taxa_sums(ps_subset)))
  
  #Simple deseq object only accounting for the differences in diet
  deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ diet) 
  
  #Performing the deseq analysis
  deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
  
  print(resultsNames(deseq_subset))
  
  #For a given taxononical levels, creates graph for each timepoint, displaying which species were found to be differentially abundant
  relabSingleTimepoint(ps_subset, deseq_subset, measure = "log2fold", "diet", timePoint = timePoint, taxa = "Species", threshold = 0.05, LDA = TRUE, FDR = TRUE, customColors = customColors, path = newPath)  
  
}

#At other taxonomic levels
taxonomicLevels <- c("Genus","Family","Order","Class","Phylum")

#Iterate through timepoints
for(timePoint in levels(sample_data(ps_flt_diet)$timepoint)){
  
  #New path created for each week
  newPath <- paste(pathToSave, "timepoint_", timePoint, "/", sep = "")
  existingDirCheck(newPath)
  
  #Creating phyloseq objects for each timepoint
  ps_subset <- prune_samples(sample_data(ps_flt_diet)$timepoint == timePoint & sample_data(ps_flt_diet)$treatment == "dss", ps_flt_diet)
  ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
  ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
  # Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
  
  # print(length(taxa_sums(ps_subset)))
  # ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
  # print(length(taxa_sums(ps_subset)))
  # 
  # # Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
  # ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
  # print(length(taxa_sums(ps_subset)))
  
  for(txnLevel in taxonomicLevels){
    
    #Creates ps subset for taxonomical level of interest
    ps_taxa <- tax_glom(ps_subset, taxrank = txnLevel)
    colnames(sample_data(ps_taxa))[2] <- "sample_id"
    deseq_subset <- phyloseq_to_deseq2(ps_taxa, ~ diet)
    deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
    
    #For a given taxononical levels, creates graph for each timepoint, displaying which species were found to be differentially abundant
    relabSingleTimepoint(ps_taxa, deseq_subset, measure = "log2fold", "diet", timePoint = timePoint, taxa = txnLevel, threshold = 0.05, customColors, newPath) 
  }
  
  
  
}
}



























# Relative abundance analysis, stackbar extended graphs
# For diet only timepoints
# First, timepoints and groups must be ordered properly and as factors
sample_data(ps_flt_diet)$diet 
sample_data(ps_flt_diet)$timepoint

# Define factor that is combination of diet and timepoint for graph visualization
sample_data(ps_flt_diet)$gg_group <- factor(paste(sample_data(ps_flt_diet)$diet, sample_data(ps_flt_diet)$timepoint, sep = ":"))
sample_data(ps_flt_diet)$gg_group

diet_phyla_fam <- plot_microbiota_timepoints(
  ps_object = ps_flt_diet,
  exp_group = "diet",
  timePoints = TRUE,
  time_variable = "timepoint",
  combined_group = 'gg_group',
  sample_name = 'full_id',
  hues = c("Purples", "Blues", "Greens", "Oranges"), # c("Purples", "Blues", "Reds", "Greens", "Oranges", "Greys", "BuPu")
  differential_analysis = T,
  sig_lab = T,
  n_row = 2,
  n_col = 3,
  threshold = 1,
  fdr_threshold = 0.05,
  main_level = "Phylum",
  sub_level = "Family",
  n_phy = 4, # number of taxa to show 
  mult_comp = F, # pairwise comparaisons for diff ab analysis
  selected_comparisons = list(c( "50:0",  "500:0"), c( "50:35",  "500:35"), c( "50:49",  "500:49")),
  showOnlySubLegend = TRUE
)

print(diet_phyla_fam$plot)
print(diet_phyla_fam$significant_table_main)
print(diet_phyla_fam$significant_table_sub)

library(ggh4x)

# Custom the plot
p <- diet_phyla_fam$plot + 
  facet_wrap2(~ gg_group, 
              scales  = "free_x", nrow = 2, ncol = 3,
              strip = strip_themed(background_x = elem_list_rect(fill = c("blue","blue","blue","red", "red", "red"))),
              labeller = as_labeller(c("50:0" = "50 ppm / 3w",
                                       "50:35" = "50 ppm / 8w",
                                       "50:49" = "50 ppm / 10w",
                                       "500:0" = "500 ppm / 3w",
                                       "500:35" = "500 ppm / 8w",
                                       "500:49"="500 ppm / 10w")))+
  theme(text = element_text(family = "Arial"),      # Global text settings
        strip.text = element_text(size = 14, face = "bold"),  # Facet titles
        plot.title = element_text(size = 20, face = "bold"),  # Main title
        axis.title = element_text(size = 15, face = "bold"),  # Axis titles
        axis.text = element_text(size = 12, face = "bold"),   # Axis text
        legend.title = element_text(face = "bold", size = 14)  # Legend title  # Legend text
  ) +
  scale_x_discrete(labels = function(x) substr(x, 1, 5))+
  labs(x = "Sample ID")
p

# Saving the plot and the associated stats
existingDirCheck("../figures/thibault/stackbar")
ggsave(plot = p, filename = "../figures/thibault/stackbar/diet_stackbar.png", width = 9, height = 7, dpi = 300)
writeStackbarExtendedSigTable(main_table = diet_phyla_fam$significant_table_main, includeSubTable = TRUE, sub_table = diet_phyla_fam$significant_table_sub, filepath = "../figures/thibault/stackbar/diet_stackbar_stats.xlsx")

# pvalues heatmap for the main lvl stats
pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/thibault/stackbar/diet_stackbar_stats.xlsx")),
            selected_comparisons = c("50:0_vs_500:0", "50:35_vs_500:35","50:49_vs_500:49"), displayChangeArrows = FALSE, displayPValues = TRUE,
            txn_lvl="Phylum", lvl = "main", taxons = diet_phyla_fam$main_names[!grepl("Others", x = diet_phyla_fam$main_names)], group = "gg_group", path)

# pvalues heatmap for the sub lvl stats
p = pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/thibault/stackbar/diet_stackbar_stats.xlsx")),
                selected_comparisons = c("50:0_vs_500:0", "50:35_vs_500:35","50:49_vs_500:49"),
                txn_lvl="Family", lvl = "sub", taxons =  diet_phyla_fam$sub_names, group = "gg_group", displayPValues = FALSE, displayChangeArrows = TRUE, path) # You can add [!grepl("Others", x = iron_exp_family$sub_names)] to remove "others"
p+scale_x_discrete(labels = c("50 vs 500 3w", "50 vs 500 8w", "50 vs 500 10w", "50 vs 500 14w"))+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 11))




# For diet_dss timepoints
# First, timepoints and groups must be ordered properly and as factors
sample_data(ps_dss_relab_flt)$diet 
sample_data(ps_dss_relab_flt)$treatment 
sample_data(ps_dss_relab_flt)$gg_group2
sample_data(ps_dss_relab_flt)$timepoint
ps_sub <- prune_samples(sample_data(ps_dss_relab_flt)$timepoint == "54", ps_dss_relab_flt)
print(length(taxa_sums(ps_sub)))
ps_sub <- prune_taxa(taxa_sums(ps_sub) > 10, ps_sub)
ps_sub <- prune_taxa(colSums(otu_table(ps_sub) > 0) >= (0.05 * nsamples(ps_sub)), ps_sub)
print(length(taxa_sums(ps_sub)))



# Selected comparisons should be a number of four and follow design as: "
diet_dss_phyla_fam <- plot_microbiota_2Fac(
  ps_object = ps_sub,
  exp_group = "gg_group2",
  twoFactor = TRUE,
  fac1 = "diet",
  refFac1 = "50",
  fac2 = "treatment",
  refFac2 = "water",
  sample_name = 'full_id',
  hues = c("Purples", "Blues", "Greens", "Oranges"), # c("Purples", "Blues", "Reds", "Greens", "Oranges", "Greys", "BuPu")
  differential_analysis = T,
  sig_lab = T,
  n_row = 2,
  n_col = 2,
  threshold = 1,
  fdr_threshold = 0.05,
  main_level = "Phylum",
  sub_level = "Family",
  n_phy = 4, # number of taxa to show 
  mult_comp = F, # pairwise comparaisons for diff ab analysis
  selected_comparisons = list(c("50:water", "50:dss"),
                              c("500:water", "500:dss"),
                              c("50:water", "500:water"),
                              c("50:dss", "500:dss")),
  showOnlySubLegend = FALSE
)

print(diet_dss_phyla_fam$plot)
print(diet_dss_phyla_fam$significant_table_main)
print(diet_dss_phyla_fam$significant_table_sub)

library(ggh4x)

# Custom the plot
p <- diet_dss_phyla_fam$plot + 
  facet_wrap2(~ gg_group2, 
              scales  = "free_x", nrow = 2, ncol = 2,
              strip = strip_themed(background_x = elem_list_rect(fill = c("blue","darkblue","red","darkred"))),
              labeller = as_labeller(c("50:water" = "50 ppm Ctrl",
                                       "50:dss" = "50 ppm DSS",
                                       "500:water" = "500 ppm Ctrl",
                                       "500:dss" = "500 ppm DSS")))+
  theme(text = element_text(family = "Arial"),      # Global text settings
        strip.text = element_text(size = 14, face = "bold", color = "white"),  # Facet titles
        plot.title = element_text(size = 20, face = "bold"),  # Main title
        axis.title = element_text(size = 15, face = "bold"),  # Axis titles
        axis.text = element_text(size = 12, face = "bold"),   # Axis text
        legend.title = element_text(face = "bold", size = 14)  # Legend title  # Legend text
  ) +
  scale_x_discrete(labels = function(x) substr(x, 1, 5))+
  labs(x = "Sample ID")
p

# Saving the plot and the associated stats
existingDirCheck("../figures/thibault/stackbar")
ggsave(plot = p, filename = "../figures/thibault/icm/diet_dss_stackbar.png", width = 7, height = 7, dpi = 300)
writeStackbarExtendedSigTable(main_table = diet_dss_phyla_fam$significant_table_main, includeSubTable = TRUE, sub_table = diet_dss_phyla_fam$significant_table_sub, filepath = "../figures/thibault/stackbar/diet_dss_stackbar_stats.xlsx")

# pvalues heatmap for the main lvl stats
pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/thibault/stackbar/diet_dss_stackbar_stats.xlsx")),
            selected_comparisons = c("50:water_vs_50:dss", "500:water_vs_500:dss","50:water_vs_500:water","50:dss_vs_500:dss"), displayChangeArrows = FALSE, displayPValues = TRUE,
            txn_lvl="Phylum", lvl = "main", taxons = diet_dss_phyla_fam$main_names[!grepl("Others", x = diet_dss_phyla_fam$main_names)], group = "gg_group2", path)

# pvalues heatmap for the sub lvl stats
p = pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/thibault/stackbar/diet_dss_stackbar_stats.xlsx")),
                selected_comparisons = c("50:water_vs_50:dss", "500:water_vs_500:dss","50:water_vs_500:water","50:dss_vs_500:dss"),
                txn_lvl="Family", lvl = "sub", taxons =  diet_dss_phyla_fam$sub_names, group = "gg_group2", displayPValues = FALSE, displayChangeArrows = TRUE, path) # You can add [!grepl("Others", x = iron_exp_family$sub_names)] to remove "others"
p+scale_x_discrete(labels = c("50 Ctrl vs 50 DSS", "500 Ctrl vs 500 DSS", "50 Ctrl vs 500 Ctrl", "50 DSS vs 500 DSS"))+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 11, face = "bold"))






sample_data(ps_dss_relab_flt)$diet 
sample_data(ps_dss_relab_flt)$treatment 
sample_data(ps_dss_relab_flt)$gg_group2
sample_data(ps_dss_relab_flt)$timepoint
ps_sub <- prune_samples(sample_data(ps_dss_relab_flt)$timepoint == "final", ps_dss_relab_flt)
sample_data(ps_sub)$diet 
sample_data(ps_sub)$treatment 
sample_data(ps_sub)$gg_group2
sample_data(ps_sub)$timepoint

length(taxa_sums(ps_flt))

# Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
# Function filtering out ASVs for which they were in total less than a threshold count
ps_sub <- prune_taxa(taxa_sums(ps_sub) > 10, ps_sub)
sum(taxa_sums(ps_flt))
length(taxa_sums(ps_flt))

# Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
ps_sub <- prune_taxa(colSums(otu_table(ps_sub) > 0) >= (0.05 * nsamples(ps_sub)), ps_sub)
sum(taxa_sums(ps_sub))
length(taxa_sums(ps_sub))

# Selected comparisons should be a number of four and follow design as: "
final_phyla_fam <- plot_microbiota_2Fac(
  ps_object = ps_sub,
  exp_group = "gg_group2",
  twoFactor = TRUE,
  fac1 = "diet",
  refFac1 = "50",
  fac2 = "treatment",
  refFac2 = "water",
  sample_name = 'full_id',
  hues = c("Purples", "Blues", "Greens", "Oranges"), # c("Purples", "Blues", "Reds", "Greens", "Oranges", "Greys", "BuPu")
  differential_analysis = T,
  sig_lab = T,
  n_row = 2,
  n_col = 2,
  threshold = 1,
  fdr_threshold = 0.05,
  main_level = "Phylum",
  sub_level = "Family",
  n_phy = 4, # number of taxa to show 
  mult_comp = F, # pairwise comparaisons for diff ab analysis
  selected_comparisons = list(c("50:water", "50:dss"),
                              c("500:water", "500:dss"),
                              c("50:water", "500:water"),
                              c("50:dss", "500:dss")),
  showOnlySubLegend = TRUE
)

print(final_phyla_fam$plot)
print(final_phyla_fam$significant_table_main)
print(final_phyla_fam$significant_table_sub)

library(ggh4x)

# Custom the plot
p <- final_phyla_fam$plot + 
  facet_wrap2(~ gg_group2, 
              scales  = "free_x", nrow = 2, ncol = 2,
              strip = strip_themed(background_x = elem_list_rect(fill = c("blue","blueviolet","red","darkorange"))),
              labeller = as_labeller(c("50:water" = "50 ppm Ctrl",
                                       "50:dss" = "50 ppm DSS",
                                       "500:water" = "500 ppm Ctrl",
                                       "500:dss" = "500 ppm DSS")))+
  theme(text = element_text(family = "Arial"),      # Global text settings
        strip.text = element_text(size = 14, face = "bold"),  # Facet titles
        plot.title = element_text(size = 20, face = "bold"),  # Main title
        axis.title = element_text(size = 15, face = "bold"),  # Axis titles
        axis.text = element_text(size = 12, face = "bold"),   # Axis text
        legend.title = element_text(face = "bold", size = 14)  # Legend title  # Legend text
  ) +
  scale_x_discrete(labels = function(x) substr(x, 1, 5))+
  labs(x = "Sample ID")
p

# Saving the plot and the associated stats
existingDirCheck("../figures/thibault/stackbar")
ggsave(plot = p, filename = "../figures/thibault/stackbar/final_stackbar.png", width = 6, height = 6, dpi = 300)
writeStackbarExtendedSigTable(main_table = final_phyla_fam$significant_table_main, includeSubTable = TRUE, sub_table = final_phyla_fam$significant_table_sub, filepath = "../figures/thibault/stackbar/final_stackbar_stats.xlsx")

# pvalues heatmap for the main lvl stats
pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/thibault/stackbar/final_stackbar_stats.xlsx")),
            selected_comparisons = c("50:water_vs_50:dss", "500:water_vs_500:dss","50:water_vs_500:water","50:dss_vs_500:dss"), displayChangeArrows = FALSE, displayPValues = TRUE,
            txn_lvl="Phylum", lvl = "main", taxons = final_phyla_fam$main_names[!grepl("Others", x = final_phyla_fam$main_names)], group = "gg_group2", path)

# pvalues heatmap for the sub lvl stats
p = pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/thibault/stackbar/final_stackbar_stats.xlsx")),
                selected_comparisons = c("50:water_vs_50:dss", "500:water_vs_500:dss","50:water_vs_500:water","50:dss_vs_500:dss"),
                txn_lvl="Family", lvl = "sub", taxons =  final_phyla_fam$sub_names, group = "gg_group2", displayPValues = FALSE, displayChangeArrows = TRUE, path) # You can add [!grepl("Others", x = iron_exp_family$sub_names)] to remove "others"
p+scale_x_discrete(labels = c("50 Ctrl vs 50 DSS", "500 Ctrl vs 500 DSS", "50 Ctrl vs 500 Ctrl", "50 DSS vs 500 DSS"))+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 11))



















