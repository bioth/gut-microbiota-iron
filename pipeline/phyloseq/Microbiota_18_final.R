# Loading required packages
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
  library(ggh4x)
  # library(ggpicrust2)
  
}

# Load custom functions for microbiota analysis
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/utilities.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/alpha_diversity_graphs_and_stats.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/beta_diversity_graphs_and_stats.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/correlation_graphs_and_stats.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/relab_analysis_graphs_and_stats.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/taxa_distrib_graphs_and_stats.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/plot_microbiota_extension.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/deseq2_log2fold_change_analysis.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/chronobiome.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/picrust2_graphs.R")

# For microbiota 18 (final data)
#set working directory
setwd("~/CHUM_git/Microbiota_18_final")
asv_table <- as.data.frame(fread("Microbiota18_final_data2/asv_table/asv_table.csv", sep = ";"))
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
  
  # Replace final by 
  metadata[metadata$timepoint=="final","timepoint"] <- "112"
  
  # Add week column 
  metadata$week <- as.character(round(as.numeric(metadata$timepoint)/7, 1))
  
  # Add age column 
  metadata$age <- as.character(round(as.numeric(metadata$timepoint)/7+3, 1))
  
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

# Load taxonomical assignments
taxa <- as.matrix(fread("taxonomy/taxa_annotation_final.csv", sep = ";"))
rownames(taxa) <- taxa[,1]  # Use the first column as row names
taxa <- taxa[,-1]  # Drop the first column

# Creating phyloseq object
ps <- phyloseq(otu_table(asv_table, taxa_are_rows = FALSE),
               tax_table(taxa), sample_data(metadata))

# Use short names for the asvs (eg ASV21) rather than full dna sequence name
# And add refseq entry to the phyloseq object
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
names(dna) <- taxa_names(ps)

# Creating phylogenetic tree
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
  write.tree(tree, file = "~/CHUM_git/Microbiota_18_final/taxonomy/phylogenetic_tree.newick")
  
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

sum(taxa_sums(ps)) # total number of reads
length(taxa_sums(ps)) # total number of ASVs
nrow(tax_table(ps))-sum(is.na(tax_table(ps)[,7])) # how many identified at species level
nrow(tax_table(ps))-sum(is.na(tax_table(ps)[,6])) # how many identified at genus level
nrow(tax_table(ps))-sum(is.na(tax_table(ps)[,5])) # how many identified at family level
nrow(tax_table(ps))-sum(is.na(tax_table(ps)[,4])) # how many identified at order level
nrow(tax_table(ps))-sum(is.na(tax_table(ps)[,3])) # how many identified at class level

# Put as factors variables that are going to be used
sample_data(ps)$gg_group2 <- factor(sample_data(ps)$gg_group2, levels = c("50:water", "500:water", "50:dss", "500:dss")) # Put gg_group2 as factor
sample_data(ps)$timepoint <- factor(sample_data(ps)$timepoint, levels = c("0","35","49","54","112")) # Put timepoint as factor
sample_data(ps)$week <- factor(sample_data(ps)$week, levels = c("0","5","7","7.7","16")) # Put week as factor
sample_data(ps)$treatment <- factor(sample_data(ps)$treatment, levels = c("water","dss")) # Put treatment as factor
sample_data(ps)$diet <- factor(sample_data(ps)$diet, levels = c("50","500")) # Put diet as factor

# Load phylogenetic tree if possible
tree <- read.tree("~/CHUM_git/Microbiota_18_final/taxonomy/phylogenetic_tree.newick")

# Add tree to phyloseq object
ps_tree <- merge_phyloseq(ps, phy_tree(tree))
phy_tree(ps_tree) <- midpoint(tree) # Root the tree

# Create single timepoint phyloseq objects and apply filter
ps_t0 <- prune_samples(sample_data(ps)$timepoint %in% c("0"), ps)
ps_t0_flt <- prune_taxa(taxa_sums(ps_t0) > 10, ps_t0)
ps_t0_flt <- prune_taxa(colSums(otu_table(ps_t0_flt) > 0) >= (0.5 * nsamples(ps_t0_flt)), ps_t0_flt)
length(taxa_sums(ps_t0_flt))

ps_t35 <- prune_samples(sample_data(ps)$timepoint %in% c("35"), ps)
ps_t35_flt <- prune_taxa(taxa_sums(ps_t35) > 10, ps_t35)
ps_t35_flt <- prune_taxa(colSums(otu_table(ps_t35_flt) > 0) >= (0.5 * nsamples(ps_t35_flt)), ps_t35_flt)
length(taxa_sums(ps_t35_flt))

ps_t49 <- prune_samples(sample_data(ps)$timepoint %in% c("49"), ps)
ps_t49_flt <- prune_taxa(taxa_sums(ps_t49) > 10, ps_t49)
ps_t49_flt <- prune_taxa(colSums(otu_table(ps_t49_flt) > 0) >= (0.5 * nsamples(ps_t49_flt)), ps_t49_flt)
length(taxa_sums(ps_t49_flt))

ps_t54 <- prune_samples(sample_data(ps)$timepoint %in% c("54"), ps)
ps_t54_flt <- prune_taxa(taxa_sums(ps_t54) > 10, ps_t54)
ps_t54_flt <- prune_taxa(colSums(otu_table(ps_t54_flt) > 0) >= (0.5 * nsamples(ps_t54_flt)), ps_t54_flt)
length(taxa_sums(ps_t54_flt))

ps_tfinal <- prune_samples(sample_data(ps)$timepoint %in% c("112"), ps)
ps_tfinal_flt <- prune_taxa(taxa_sums(ps_tfinal) > 10, ps_tfinal)
ps_tfinal_flt <- prune_taxa(colSums(otu_table(ps_tfinal_flt) > 0) >= (0.5 * nsamples(ps_tfinal_flt)), ps_tfinal_flt)
length(taxa_sums(ps_tfinal_flt))

# Create phyloseq obejcts that we need for the analysis
ps_diet <- merge_phyloseq(ps_t0, ps_t35, ps_t49)
ps_dss_alpha <- merge_phyloseq(ps_t49, ps_t54, ps_tfinal)
ps_dss_relab_flt <- merge_phyloseq(ps_t49_flt, ps_t54_flt, ps_tfinal_flt)
ps_dss <- merge_phyloseq(ps_t54, ps_tfinal)
ps_flt_diet <- merge_phyloseq(ps_t0_flt, ps_t35_flt, ps_t49_flt)
ps_flt_dss <- merge_phyloseq(ps_t54_flt, ps_tfinal_flt)
ps_flt_all <- merge_phyloseq(ps_t0_flt, ps_t35_flt, ps_t49_flt, ps_t54_flt, ps_tfinal_flt)

# Alpha diveristy
{
  existingDirCheck("../figures/Thibault_dss/diet/")
  graphs = alphaDiversityTimeSeries2(ps_diet, "../figures/Thibault_dss/diet/", time = "timepoint", group = "diet", writeData = TRUE)

  # Stats
  alpha_d <- read.xlsx("../figures/Thibault_dss/diet/alpha_diversity/alpha_diversity_data.xlsx")
  alpha_d$diet <- factor(alpha_d$diet, levels = c("50", "500"))
  alpha_d$timepoint <- factor(alpha_d$timepoint, levels = c("0", "35", "49"))
  
  # Chao1
  graphs[[1]]+
    scale_fill_manual(values = c("blue","red"),
                      labels = c("50 ppm","500 ppm"))+
    scale_x_discrete(labels = c("3 weeks\n(Weaning)","5 weeks", "7 weeks"))+
    labs(y = "Chao1 Index", x = "")+
    my_theme()+
    ylim(0,1800)+
    geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(0.8),           # left box in each timepoint
      xmax = c(1.2),
      annotations = "n.s.",
      y_position = c(1500), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = -0.1,
    )+geom_signif( # For second timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(1.8),           # left box in each timepoint
      xmax = c(2.2),
      annotations = "n.s.",
      y_position = c(1500), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+geom_signif( # For third timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(2.8),           # left box in each timepoint
      xmax = c(3.2),
      annotations = "p=0.08",
      y_position = c(1500), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = -0.1,
    )
  ggsave("../figures/Thibault_dss/diet/alpha_diversity/chao1.png",
         bg = "white",height = 4, width =5, dpi = 300)
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "0",], group = "diet", measure = "Chao1")
  wilcox.test(Chao1 ~ diet, data = alpha_d[alpha_d$timepoint == "0",])
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "35",], group = "diet", measure = "Chao1")
  t.test(Chao1 ~ diet, data = alpha_d[alpha_d$timepoint == "35",], var.equal = TRUE)
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "diet", measure = "Chao1")
  t.test(Chao1 ~ diet, data = alpha_d[alpha_d$timepoint == "49",], var.equal = TRUE)
  
  # Shannon
  graphs[[2]]+
    scale_fill_manual(values = c("blue","red"),
                      labels = c("50 ppm","500 ppm"))+
    scale_x_discrete(labels = c("3 weeks\n(Weaning)","5 weeks", "7 weeks"))+
    labs(y = "Shannon Index", x = "")+
    my_theme()+
  ylim(0,6)+
    geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(0.8),           # left box in each timepoint
      xmax = c(1.2),
      annotations = "n.s.",
      y_position = c(5.5), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = -0.1,
    )+geom_signif( # For second timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(1.8),           # left box in each timepoint
      xmax = c(2.2),
      annotations = "n.s.",
      y_position = c(5.5), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+geom_signif( # For third timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(2.8),           # left box in each timepoint
      xmax = c(3.2),
      annotations = "n.s.",
      y_position = c(5.2), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = -0.1,
    )
  ggsave("../figures/Thibault_dss/diet/alpha_diversity/shannon.png",
         bg = "white",height = 4, width =5, dpi = 300)
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "0",], group = "diet", measure = "Shannon")
  t.test(Shannon ~ diet, data = alpha_d[alpha_d$timepoint == "0",], var.equal = TRUE)
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "35",], group = "diet", measure = "Shannon")
  t.test(Shannon ~ diet, data = alpha_d[alpha_d$timepoint == "35",], var.equal = TRUE)
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "diet", measure = "Shannon")
  wilcox.test(Shannon ~ diet, data = alpha_d[alpha_d$timepoint == "49",])
  
  # InvSimpson
  graphs[[3]]+
    scale_fill_manual(values = c("blue","red"),
                      labels = c("50 ppm","500 ppm"))+
    scale_x_discrete(labels = c("3 weeks\n(Weaning)","5 weeks", "7 weeks"))+
    labs(y = "Inverse Simpson", x = "")+
    my_theme()+
    ylim(0,45)+
    geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(0.8),           # left box in each timepoint
      xmax = c(1.2),
      annotations = "n.s.",
      y_position = c(41), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = -0.1,
    )+geom_signif( # For second timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(1.8),           # left box in each timepoint
      xmax = c(2.2),
      annotations = "n.s.",
      y_position = c(32), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+geom_signif( # For third timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(2.8),           # left box in each timepoint
      xmax = c(3.2),
      annotations = "n.s.",
      y_position = c(20), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = -0.1,
    )
  ggsave("../figures/Thibault_dss/diet/alpha_diversity/invsimpson.png",
         bg = "white",height = 4, width =5, dpi = 300)
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "0",], group = "diet", measure = "InvSimpson")
  t.test(InvSimpson ~ diet, data = alpha_d[alpha_d$timepoint == "0",], var.equal = TRUE)
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "35",], group = "diet", measure = "InvSimpson")
  wilcox.test(InvSimpson ~ diet, data = alpha_d[alpha_d$timepoint == "35",])
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "diet", measure = "InvSimpson")
  wilcox.test(InvSimpson ~ diet, data = alpha_d[alpha_d$timepoint == "49",])
  
  # Alpha diveristy for dss_diet only
  existingDirCheck("../figures/Thibault_dss/diet_dss/")
  graphs = alphaDiversityTimeSeries2(ps_dss_alpha, "../figures/Thibault_dss/diet_dss/", time = "timepoint", group = "gg_group2", writeData = TRUE)
  
  # Load data for stats calculations
  alpha_d <- read.xlsx("../figures/Thibault_dss/diet_dss/alpha_diversity/alpha_diversity_data.xlsx")
  alpha_d$gg_group2 <- factor(alpha_d$gg_group2, levels = c("50:water","500:water","50:dss","500:dss"))
  alpha_d$timepoint <- factor(alpha_d$timepoint, levels = c("49", "54", "final"))
  
  # Chao1
  graphs[[1]]+
    scale_fill_manual(values = c("blue","red","darkblue","darkred"),
                      labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))+
    scale_pattern_manual(values = c("circle","stripe","circle","stripe"))+
    scale_x_discrete(labels = c("DSS day 0","DSS day 5", "End of\nrecovery"),
                     expand = c(0,0.5))+
    labs(y = "Species Richness", x = "")+
    guides(pattern = "none")+
    ylim(0,NA)+
    my_theme()
    geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(0.7),           # left box in each timepoint
      xmax = c(1.3),
      annotations = "n.s.",
      y_position = c(1500), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+
    geom_signif( # For last timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(2.7,3.1,2.7,2.9),           # left box in each timepoint
      xmax = c(2.9,3.3,3.1,3.3),
      annotations = c("n.s.","n.s.","p=0.08","p=0.08"),
      y_position = c(1500,1100,1750,2000), #
      tip_length = c(0.02,0.02,0.02,0.02),
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+
    geom_signif( # For second timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(1.7,2.1,1.7,1.9), # left position
      xmax = c(1.9,2.3,2.1,2.3), # left position
      annotations = c("n.s.","n.s.","***","***"),
      y_position = c(1350,1200,1600,1800), #
      tip_length = c(0.02,0.02,0.02,0.02),
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )
  ggsave("../figures/Thibault_dss/diet_dss/alpha_diversity/chao1.png",
         bg = "white",height = 5, width =7, dpi = 300)
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "gg_group2", measure = "Chao1")
  TukeyHSD(aov(Chao1 ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "49",]))
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "54",], group = "gg_group2", measure = "Chao1")
  TukeyHSD(aov(Chao1 ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "54",]))
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "final",], group = "gg_group2", measure = "Chao1")
  TukeyHSD(aov(Chao1 ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "final",]))
  
  # Shannon
  graphs[[2]]+
    scale_fill_manual(values = c("blue","red","darkblue","darkred"),
                      labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))+
    scale_pattern_manual(values = c("circle","stripe","circle","stripe"))+
    scale_x_discrete(labels = c("DSS day 0","DSS day 5", "End of\nrecovery"),
                     expand = c(0,0.5))+
    labs(y = "Shannon Index", x = "", title = "Shannon index")+
    guides(pattern = "none")+
    ylim(0,NA)+
    my_theme()
    geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(0.7),           # left box in each timepoint
      xmax = c(1.3),
      annotations = "n.s.",
      y_position = c(4), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+
    geom_signif( # For last timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(2.7,3.1,2.7,2.9),           # left box in each timepoint
      xmax = c(2.9,3.3,3.1,3.3),
      annotations = c("n.s.","p=0.09","n.s.","p=0.09"),
      y_position = c(4,4,4.5,5), #
      tip_length = c(0.02,0.02,0.02,0.02),
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+
    geom_signif( # For second timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(1.7,2.1,1.7,1.9), # left position
      xmax = c(1.9,2.3,2.1,2.3), # left position
      annotations = c("n.s.","n.s.","*","p=0.08"),
      y_position = c(4.2,4.2,4.5,4.9), #
      tip_length = c(0.02,0.02,0.02,0.02),
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )
  ggsave("../figures/Thibault_dss/diet_dss/alpha_diversity/shannon.png",
         bg = "white",height = 5, width =7, dpi = 300)
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "gg_group2", measure = "Shannon")
  pairwise.wilcox.test(alpha_d[alpha_d$timepoint == "49",]$Shannon, alpha_d[alpha_d$timepoint == "49",]$gg_group2, p.adjust.method = "BH")
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "54",], group = "gg_group2", measure = "Shannon")
  TukeyHSD(aov(Shannon ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "54",]))
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "final",], group = "gg_group2", measure = "Shannon")
  TukeyHSD(aov(Shannon ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "final",]))
  pairwise.wilcox.test(alpha_d[alpha_d$timepoint == "final",]$Shannon, alpha_d[alpha_d$timepoint == "final",]$gg_group2, p.adjust.method = "BH")
  
  
  # InvSimpson
  graphs[[3]]+
    scale_fill_manual(values = c("blue","red","darkblue","darkred"),
                      labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))+
    scale_pattern_manual(values = c("circle","stripe","circle","stripe"))+
    scale_x_discrete(labels = c("DSS day 0","DSS day 5", "End of\nrecovery"),
                     expand = c(0,0.5))+
    labs(y = "Inverse Simpson Index", x = "", title = "Inverse Simpson index")+
    guides(pattern = "none")+
    ylim(0,NA)+
    my_theme()
    geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(0.7),           # left box in each timepoint
      xmax = c(1.3),
      annotations = "n.s.",
      y_position = c(14), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+    
    geom_signif( # For second timepoint
    # comparisons = list(c(groups[1],groups[4])),
    xmin = c(1.7,2.1,1.7,1.9), # left position
    xmax = c(1.9,2.3,2.1,2.3), # left position
    annotations = c("n.s.","n.s.","***","**"),
    y_position = c(12,18.5,20.5,23), #
    tip_length = c(0.02,0.02,0.02,0.02),
    color = "black",
    size = 0.5,
    textsize = 4,
    margin_top = 0.1, # Moves the top according to this value
    vjust = 0,
  )+
    geom_signif( # For last timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(2.7,3.1,2.7,2.9),           # left box in each timepoint
      xmax = c(2.9,3.3,3.1,3.3),
      annotations = c("n.s.","p=0.063","n.s.","p=0.063"),
      y_position = c(10.5,16,18.5,20.5), #
      tip_length = c(0.02,0.02,0.02,0.02),
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )
  ggsave("../figures/Thibault_dss/diet_dss/alpha_diversity/InvSimpson.png",
         bg = "white",height = 5, width =7, dpi = 300)
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "gg_group2", measure = "InvSimpson")
  pairwise.wilcox.test(alpha_d[alpha_d$timepoint == "49",]$InvSimpson, alpha_d[alpha_d$timepoint == "49",]$gg_group2, p.adjust.method = "BH")
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "54",], group = "gg_group2", measure = "InvSimpson")
  pairwise.wilcox.test(alpha_d[alpha_d$timepoint == "54",]$InvSimpson, alpha_d[alpha_d$timepoint == "54",]$gg_group2, p.adjust.method = "BH")
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "final",], group = "gg_group2", measure = "InvSimpson")
  pairwise.wilcox.test(alpha_d[alpha_d$timepoint == "final",]$InvSimpson, alpha_d[alpha_d$timepoint == "final",]$gg_group2, p.adjust.method = "BH")
  
  
}

# Alpha diveristy - but timeline
{
  existingDirCheck("../figures/Thibault_dss/alpha_diversity_timeline/")
  
  # Week as numeric
  sample_data(ps)$timepoint  <- as.numeric(as.character(sample_data(ps)$timepoint))
  sample_data(ps)$gg_group2 
  sample_data(ps)$gg_group2 <- factor(sample_data(ps)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
  
  sample_data(ps)$gg_group3 <- ifelse(
    sample_data(ps)$week %in% c(3, 8),
    as.character(sample_data(ps)$diet),
    as.character(sample_data(ps)$gg_group2)
  )
  sample_data(ps)$gg_group3 <- factor(sample_data(ps)$gg_group3, levels = c("50","500","50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
  custom_colors <- c("blue","red","blue","red","darkblue", "darkred")
  
  #Estinate richness measures for dataset
  richness_data <- estimate_richness(ps, measures = c("Chao1", "Shannon", "InvSimpson", "Observed"))
  alpha_d <- cbind(as.data.frame(sample_data(ps)), richness_data)
  
  custom_colors <- c("blue","red","darkblue", "darkred")
  
  
  graphs <- alphaDiversityTimeline(ps, time = "timepoint", group = "gg_group2", custom_colors, semRibbons = TRUE)
  graphs <- alphaDiversityTimeline(ps, time = "week", group = "gg_group3", custom_colors)
  
  # Observed
  p1 <- graphs[[4]]+
    scale_x_continuous(
      breaks = seq(min(sample_data(ps)$timepoint), max(sample_data(ps)$timepoint), by = 7), # breaks every 7 days
      labels = function(x) paste0("T", x)                    # "T" before each label
    ) +
    labs(y = "Species observed", x = "Timepoint", color = "Group", fill = "", title = 'Species observed')+
    guides(fill = "none")+
    theme(axis.ticks.x = element_line(),
          legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5),
          axis.text.x = element_text(size = 6)
          )+
    ylim(0,NA)
  
  ggsave("../figures/Thibault_dss/alpha_diversity_timeline/observed.png",
         bg = "white",height = 3, width =7, dpi = 300)
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "gg_group2", measure = "Chao1")
  TukeyHSD(aov(Chao1 ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "49",]))
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "54",], group = "gg_group2", measure = "Chao1")
  TukeyHSD(aov(Chao1 ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "54",]))
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "final",], group = "gg_group2", measure = "Chao1")
  TukeyHSD(aov(Chao1 ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "final",]))
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "54",], group = "gg_group2", measure = "Observed")
  TukeyHSD(aov(Observed ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "54",]))
  
  # Shannon
  p2 <- graphs[[2]]+
    scale_x_continuous(
      breaks = seq(min(sample_data(ps)$timepoint), max(sample_data(ps)$timepoint), by = 7), # breaks every 7 days
      labels = function(x) paste0("T", x)                    # "T" before each label
    ) +
    labs(y = "Shannon index", x = "Timepoint", color = "Group", fill = "", title = 'Shannon')+
    guides(fill = "none")+
    theme(axis.ticks.x = element_line(),
          legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5),
          axis.text.x = element_text(size = 6))+
    ylim(0,NA)
  
  ggsave("../figures/Thibault_dss/alpha_diversity_timeline/shannon.png",
         bg = "white",height = 3, width =7, dpi = 300)
  
  # Inverse Simpson
  p3 <- graphs[[3]]+
    scale_x_continuous(
      breaks = seq(min(sample_data(ps)$timepoint), max(sample_data(ps)$timepoint), by = 7), # breaks every 7 days
      labels = function(x) paste0("T", x)                    # "T" before each label
    ) +
    labs(y = "Inverse Simpson", x = "Timepoint", color = "Group", fill = "", title = 'Inverse Simpson')+
    guides(fill = "none")+
    theme(axis.ticks.x = element_line(),
          legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5),
          axis.text.x = element_text(size = 6))+
    ylim(0,NA)
  
  ggsave("../figures/Thibault_dss/alpha_diversity_timeline/invsimpson.png",
         bg = "white",height = 3, width =7, dpi = 300)
  
  final_plot <- p1 + p2 +p3 +plot_layout(guides = "collect")
  final_plot
  ggsave("../figures/Thibault_dss/alpha_diversity_timeline/figure.png",
         bg = "white",height = 5, width =15, dpi = 300)
}

# Beta diversity
{
  set.seed(100)
  
  phy_tree(ps_flt_diet) <- midpoint(tree) # Add rooted tree to phyloseq object
  
  # For diet timepoints
  sample_data(ps_flt_diet)$diet <- factor(sample_data(ps_flt_diet)$diet, labels = c("50 ppm","500 ppm"))
  
  # Weighted unifrac
  betaDiversityTimepoint2Factors(ps_flt_diet, sample_id = "sample_id", timeVariable = "timepoint",
                                 varToCompare =  "diet", distMethod ="wunifrac",
                                 transform = "rel_ab", customColors = c("blue","red"),
                                 font = "Arial", path = "../figures/Thibault_dss/beta_diversity/diet/", 
                                 additionnalAes = my_theme()+theme(plot.title = element_text(size = 16)), dim = c(4,12), displayPValue = TRUE,
                                 combineGraphs = TRUE)

  # For diet + treatment at t54 and tfinal
  phy_tree(ps_dss_relab_flt) <- midpoint(tree) # Add rooted tree to phyloseq object
  sample_data(ps_dss_relab_flt)$gg_group2 <- factor(sample_data(ps_dss_relab_flt)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl", "50 ppm DSS", "500 ppm DSS"))
  # Weighted unifrac => not very informative
  betaDiversityTimepoint2Factors(ps_dss_relab_flt, sample_id = "sample_id", timeVariable = "timepoint",
                                 varToCompare =  "gg_group2", distMethod ="wunifrac",
                                 customColors = c("blue","red","darkblue","darkred"),
                                 font = "Arial", path = "../figures/Thibault_dss/beta_diversity/diet_dss/",
                                 additionnalAes = my_theme()+theme(plot.title = element_text(size = 16)), dim = c(4,12), combineGraphs = TRUE)
  
  # Weighted unifrac for dss groups - t54 and last timepoint
  ps_sub <- prune_samples(sample_data(ps_flt_dss)$treatment == "dss", ps_flt_dss)
  sample_data(ps_sub)$diet <- factor(sample_data(ps_sub)$diet, labels = c("50 ppm DSS","500 ppm DSS"))
  phy_tree(ps_sub) <- midpoint(tree) # Add rooted tree to phyloseq object
  betaDiversityTimepoint2Factors(ps_sub , sample_id = "sample_id", timeVariable = "timepoint",
                                 varToCompare =  "diet", distMethod ="wunifrac",
                                 customColors = c("darkblue","darkred"),
                                 font = "Arial", path = "../figures/Thibault_dss/beta_diversity/dss_only/",
                                 additionnalAes = my_theme()+theme(plot.title = element_blank()), dim = c(2.5,4))

  # dbRDA method at t54 and for DSS groups
  betaDiversityTimepointsGroupedDbRDA(ps_sub, sample_id = "sample_id", varToCompare = "diet", formula = "diet",
                                      transform = "none", distMethod = "bray", customColors = c("darkblue","darkred"),
                                      font = "Arial", path = "../figures/Thibault_dss/dss_only/beta_diversity/dbRDA_t54/", dim = c(3,4),
                                      additionnalAes = my_theme(), displayPValue = TRUE)
  
  # dbRDA method at last timepoint and for DSS groups
  ps_sub <- prune_samples(sample_data(ps_flt_dss)$timepoint == "final" & sample_data(ps_flt_dss)$treatment == "dss", ps_flt_dss)
  sample_data(ps_sub)$diet <- factor(sample_data(ps_sub)$diet, labels = c("50 ppm DSS","500 ppm DSS"))
  betaDiversityTimepointsGroupedDbRDA(ps_sub, sample_id = "sample_id", varToCompare = "diet", formula = "diet",
                                      transform = "none", distMethod = "bray", customColors = c("darkblue","darkred"),
                                      font = "Arial", path = "../figures/Thibault_dss/dss_only/beta_diversity/dbRDA_tfinal/", dim = c(3,4),
                                      additionnalAes = 
                                        my_theme(), displayPValue = TRUE)
  
  # dbRDA method at last timepoint and for all groups
  ps_sub <- prune_samples(sample_data(ps_flt_dss)$timepoint == "final", ps_flt_dss)
  sample_data(ps_sub)$gg_group2
  betaDiversityTimepointsGroupedDbRDA(ps_sub, sample_id = "sample_id", varToCompare = "gg_group2", formula = "diet+treatment",
                                      transform = "none", distMethod = "bray", customColors = c("blue","red","darkblue","darkred"),
                                      font = "Arial", path = "../figures/Thibault_dss/dss_only/beta_diversity/dbRDA-all-groups/",
                                      additionnalAes = my_theme(), dim = c(4,5))
  
  # dbRDA method at t54 and last timepoint and for DSS groups
  ps_sub <- prune_samples(sample_data(ps_flt_dss)$treatment == "dss", ps_flt_dss)
  sample_data(ps_sub)$diet
  sample_data(ps_sub)$gg_group <- factor(sample_data(ps_sub)$gg_group, levels = c("54:50:dss","54:500:dss","final:50:dss" ,"final:500:dss"))
  betaDiversityTimepointsGroupedDbRDA(ps_sub, sample_id = "sample_id", varToCompare = "gg_group", formula = "diet*timepoint",
                                      transform = "none", distMethod = "bray", customColors = c("blue","red","darkblue","darkred"),
                                      font = "Arial", path = "../figures/Thibault_dss/dss_only/beta_diversity/dbRDA_2Last/",
                                      additionnalAes = my_theme(), dim = c(4,5))
}

# Relative abundance analysis: finding differential abundant bugs at the species level, for diet groups only
{
  # Path where to save graphs
  pathToSave <- "~/CHUM_git/figures/Thibault_dss/newTaxAnnotation/relative_abundance_diet/"
  existingDirCheck(pathToSave)
  
  #customColors for graph display
  customColors = c("blue","red")
  customPhylaColors = c("#e6550d","#31a354", "#583093")
  
  #Iterate through timepoints
  for(timePoint in levels(sample_data(ps_flt_diet)$timepoint)){
    
    #New path created for each week
    newPath <- paste(pathToSave, "timepoint_", timePoint, "/", sep = "")
    existingDirCheck(newPath)
    
    #Creating phyloseq objects for each timepoint
    ps_subset <- prune_samples(sample_data(ps_flt_diet)$timepoint == timePoint, ps_flt_diet)
    ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)
    print(length(taxa_sums(ps_subset)))
    
    #Simple deseq object only accounting for the differences in diet
    deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ diet) 
    deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric") #Performing the deseq analysis
    print(resultsNames(deseq_subset))
    
    #For a given taxononical levels, creates graph for each timepoint, displaying which species were found to be differentially abundant
    relabSingleTimepoint(ps_subset, deseq_subset, measure = "log2fold", "diet", timePoint = timePoint, taxa = "Species", threshold = 0.05, LDA = FALSE, FDR = TRUE, customColors = customColors, path = newPath)
    
    # log2fold change graph based on deseq2 results
    log2foldChangeGraphSingleTimepoint(ps_subset, deseq_subset, timePoint = timePoint, taxa = "Species", threshold = 0.05, customColors = customColors, customPhylaColors = customPhylaColors, path = newPath, dim =c(6,10))
    
  }
  
  #At other taxonomic levels
  taxonomicLevels <- c("Genus","Family","Order","Class","Phylum")
  
  #Iterate through timepoints
  for(timePoint in levels(sample_data(ps_flt_diet)$timepoint)){
    
    #New path created for each week
    newPath <- paste(pathToSave, "timepoint_", timePoint, "/", sep = "")
    existingDirCheck(newPath)
    
    #Creating phyloseq objects for each timepoint
    ps_subset <- prune_samples(sample_data(ps_flt_diet)$timepoint == timePoint, ps_flt_diet)
    ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)
    print(length(taxa_sums(ps_subset)))
    
    for(txnLevel in taxonomicLevels){
      
      #Creates ps subset for taxonomical level of interest
      ps_taxa <- tax_glom(ps_subset, taxrank = txnLevel)
      colnames(sample_data(ps_taxa))[2] <- "sample_id"
      deseq_subset <- phyloseq_to_deseq2(ps_taxa, ~ diet)
      deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
      
      #For a given taxononical levels, creates graph for each timepoint, displaying which species were found to be differentially abundant
      relabSingleTimepoint(ps_taxa, deseq_subset, measure = "log2fold", "diet", timePoint = timePoint, taxa = txnLevel, threshold = 0.05, LDA = FALSE, FDR = TRUE, customColors = customColors, path = newPath) 
    }
  }
}

# Relative abundance analysis: finding differential abundant bugs at the species level, all groups, for t49 t54 and tfinal timepoints
# Path where to save graphs
{
  pathToSave <- "~/CHUM_git/figures/Thibault_dss/newTaxAnnotation/relative_abundance_dss_diet_all_groups/"
  existingDirCheck(pathToSave)
  
  #customColors for graph display
  customColors = c("blue", "red", "darkblue","darkred")
  
  #Iterate through timepoints
  for(timePoint in levels(sample_data(ps_dss_relab_flt)$timepoint)){
    
    #New path created for each week
    newPath <- paste(pathToSave, "timepoint_", timePoint, "/", sep = "")
    existingDirCheck(newPath)
    
    #Creating phyloseq objects for each timepoint
    ps_subset <- prune_samples(sample_data(ps_dss_relab_flt)$timepoint == timePoint, ps_dss_relab_flt)
    ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)
    print(length(taxa_sums(ps_subset)))
    
    # DESeq analysis for full model + correcting for cage effects
    deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ treatment+diet+diet:treatment) 
    deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric") #Performing the deseq analysis
    # res <- results(deseq_subset)
    # View(res[is.na(res$padj), ])
    print(resultsNames(deseq_subset))
    
    #For a given taxononical levels, creates graph for each timepoint, displaying which species were found to be differentially abundant
    relabGroups(ps_subset, deseq_subset, measure = "log2fold", gg_group = "gg_group2", taxa = "Species", threshold = 0.05, FDR = TRUE,
                returnSigAsvs = FALSE, normalizeCounts = FALSE, customColors = customColors,
                pairs = list(
                  list("50:water","500:water"), list("50:dss","500:dss"),
                  list("50:water","50:dss"), list("500:water","500:dss")),
                path = newPath, single_factor_design = FALSE,
                dim = c(4,4.5), displayPvalue = FALSE, displaySignificance = TRUE, additionnalAes =
                  list(scale_x_discrete(labels = c("50 ppm\nCtrl","500 ppm\nCtrl","50 ppm\nDSS","500 ppm\nDSS")),
                       my_theme(),
                       labs(color = "", x=""))) # Include axis lines  # Include axis bar)
  }
  
  #customColors for graph display
  customColors = c("blue", "red", "darkblue","darkred")
  
  #At other taxonomic levels
  taxonomicLevels <- c("Genus","Family","Order","Class","Phylum")
  
  #Iterate through timepoints
  for(timePoint in levels(sample_data(ps_dss_relab_flt)$timepoint)){
    
    #New path created for each week
    newPath <- paste(pathToSave, "timepoint_", timePoint, "/", sep = "")
    existingDirCheck(newPath)
    
    #Creating phyloseq objects for each timepoint
    ps_subset <- prune_samples(sample_data(ps_dss_relab_flt)$timepoint == timePoint, ps_dss_relab_flt)
    ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)
    print(length(taxa_sums(ps_subset)))
    
    for(txnLevel in taxonomicLevels){
      
      #Creates ps subset for taxonomical level of interest
      ps_taxa <- tax_glom(ps_subset, taxrank = txnLevel)
      colnames(sample_data(ps_taxa))[2] <- "sample_id"
      deseq_subset <- phyloseq_to_deseq2(ps_taxa, ~ treatment+diet+diet:treatment)
      deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
      
      print(resultsNames(deseq_subset))
      
      #For a given taxononical levels, creates graph for each timepoint, displaying which species were found to be differentially abundant
      relabGroups(ps_taxa, deseq_subset, measure = "log2fold", gg_group = "gg_group2", taxa = txnLevel, threshold = 0.05, FDR = TRUE,
                  returnSigAsvs = FALSE, normalizeCounts = FALSE, customColors = customColors, 
                  pairs = list(
                    list("50:water","500:water"), list("50:dss","500:dss"),
                    list("50:water","50:dss"), list("500:water","500:dss")),
                  path = newPath, single_factor_design = FALSE,
                  dim = c(5,5), displayPvalue = FALSE, displaySignificance = TRUE, additionnalAes =
                    list(scale_x_discrete(labels = c("50 ppm\nCtrl","500 ppm\nCtrl","50 ppm\nDSS","500 ppm\nDSS")),
                         my_theme(),
                         labs(color = "", x="")))  
    }
  }
}

# Relative abundance analysis: finding differential abundant bugs at different taxonomical levels - only DSS groups
{
  pathToSave <- "~/CHUM_git/figures/Thibault_dss/newTaxAnnotation/relative_abundance_dss/"
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
    ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)
    print(length(taxa_sums(ps_subset)))
    
    #Simple deseq object only accounting for the differences in diet
    deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ diet) 
    deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric") #Performing the deseq analysis
    print(resultsNames(deseq_subset))
    
    #For a given taxononical levels, creates graph for each timepoint, displaying which species were found to be differentially abundant
    relabSingleTimepoint(ps_subset, deseq_subset, measure = "log2fold", varToCompare = "diet",
                         timePoint = timePoint, taxa = "Species", threshold = 0.05, FDR = TRUE, blockFactor = FALSE,
                         LDA = TRUE, customColors = customColors, path = newPath, displayPvalue = FALSE, additionnalAes = my_theme(), displaySampleID = FALSE)  
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
    ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)
    print(length(taxa_sums(ps_subset)))
    
    for(txnLevel in taxonomicLevels){
      
      #Creates ps subset for taxonomical level of interest
      ps_taxa <- tax_glom(ps_subset, taxrank = txnLevel)
      colnames(sample_data(ps_taxa))[2] <- "sample_id"
      deseq_subset <- phyloseq_to_deseq2(ps_taxa, ~ diet)
      deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
      print(resultsNames(deseq_subset))
      
      #For a given taxononical levels, creates graph for each timepoint, displaying which species were found to be differentially abundant
      relabSingleTimepoint(ps_taxa, deseq_subset, measure = "log2fold", "diet", timePoint = timePoint, taxa = txnLevel, threshold = 0.05, FDR = TRUE, LDA = FALSE, customColors, newPath,
                           additionnalAes = my_theme(), blockFactor = FALSE,
                           dim = c(3.5,3.5), displayPvalue = TRUE) 
    }
  }
}

# Relative abundance analysis, stackbar extended graphs
{
# For diet only timepoints
# Define factor that is combination of diet and timepoint for graph visualization
sample_data(ps_flt_diet)$diet 
sample_data(ps_flt_diet)$diet <- factor(sample_data(ps_flt_diet)$diet, labels = c("50","500"))
sample_data(ps_flt_diet)$gg_group <- factor(paste(sample_data(ps_flt_diet)$diet, sample_data(ps_flt_diet)$timepoint, sep = ":"))
sample_data(ps_flt_diet)$gg_group

diet_phyla_fam <- plot_microbiota_timepoints(
  ps_object = ps_flt_diet,
  exp_group = "diet",
  timePoints = TRUE,
  time_variable = "timepoint",
  combined_group = 'gg_group',
  sample_name = 'sample_id',
  hues = c("Blues","Greens","Purples","Oranges"), # c("Purples", "Blues", "Reds", "Greens", "Oranges", "Greys", "BuPu")
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
  showOnlySubLegend = FALSE
)

print(diet_phyla_fam$plot)
print(diet_phyla_fam$significant_table_main)
print(diet_phyla_fam$significant_table_sub)

# Define custom x-axis labels for each facet
id3 = unique(substring(sample_data(ps_flt_diet)[sample_data(ps_flt_diet)$week == "3" & sample_data(ps_flt_diet)$diet == "50",]$sample_id, 1,5))
id8 = unique(substring(text = sample_data(ps_flt_diet)[sample_data(ps_flt_diet)$week == "8" & sample_data(ps_flt_diet)$diet == "50",]$sample_id, 1,5))
list1 = which(id8 %in% id3)

id3 = unique(substring(sample_data(ps_flt_diet)[sample_data(ps_flt_diet)$week == "3" & sample_data(ps_flt_diet)$diet == "500",]$sample_id, 1,5))
id8 = unique(substring(text = sample_data(ps_flt_diet)[sample_data(ps_flt_diet)$week == "8" & sample_data(ps_flt_diet)$diet == "500",]$sample_id, 1,5))
list2 = which(id8 %in% id3)+24

facet_scales <- list(
  scale_x_discrete(labels = as.character(list1)),
  scale_x_discrete(labels = as.character(1:24)),
  scale_x_discrete(labels = as.character(1:24)),
  scale_x_discrete(labels = as.character(list2)),
  scale_x_discrete(labels = as.character(25:48)),
  scale_x_discrete(labels = as.character(25:48))
)

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
        strip.text = element_text(size = 14, face = "bold", color = "white"),  # Facet titles
        plot.title = element_text(size = 20, face = "bold"),  # Main title
        axis.title = element_text(size = 15, face = "bold"),  # Axis titles
        axis.text = element_text(size = 12, face = "bold"),   # Axis text
        legend.title = element_text(face = "bold", size = 14)  # Legend title  # Legend text
  ) +
  facetted_pos_scales(x = facet_scales)+
  scale_x_discrete(labels = function(x) substr(x, 1, 5))+
  labs(x = "Mouse ID")
p

# Saving the plot and the associated stats
existingDirCheck("../figures/Thibault_dss/stackbar")
ggsave(plot = p, filename = "../figures/Thibault_dss/stackbar/diet_stackbar.png", width = 9, height = 7, dpi = 300)
writeStackbarExtendedSigTable(main_table = diet_phyla_fam$significant_table_main, includeSubTable = TRUE, sub_table = diet_phyla_fam$significant_table_sub, filepath = "../figures/Thibault_dss/new_filtering/stackbar/diet_stackbar_stats.xlsx")

# pvalues heatmap for the main lvl stats
pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/thibault/stackbar/diet_stackbar_stats.xlsx")),
            selected_comparisons = c("50:0_vs_500:0", "50:35_vs_500:35","50:49_vs_500:49"), displayChangeArrows = FALSE, displayPValues = TRUE,
            txn_lvl="Phylum", lvl = "main", taxons = diet_phyla_fam$main_names[!grepl("Others", x = diet_phyla_fam$main_names)], group = "gg_group", path)

# pvalues heatmap for the sub lvl stats
p = pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_dss/new_filtering/stackbar/diet_stackbar_stats.xlsx")),
                selected_comparisons = c("50:0_vs_500:0", "50:35_vs_500:35","50:49_vs_500:49"),
                txn_lvl="Family", lvl = "sub", taxons =  diet_phyla_fam$sub_names, group = "gg_group", displayPValues = FALSE, displayChangeArrows = TRUE, path) # You can add [!grepl("Others", x = iron_exp_family$sub_names)] to remove "others"
p+scale_x_discrete(labels = c("50 vs 500 3w", "50 vs 500 8w", "50 vs 500 10w", "50 vs 500 14w"))+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 11))




# For diet_dss timepoints
# First, timepoints and groups must be ordered properly and as factors
sample_data(ps_t54_flt)$diet 
sample_data(ps_t54_flt)$treatment 
sample_data(ps_t54_flt)$gg_group2 <- factor(sample_data(ps_t54_flt)$gg_group2, levels = c("50:water", "50:dss", "500:water","500:dss"))
sample_data(ps_t54_flt)$timepoint

# Selected comparisons should be a number of four and follow design as: "
diet_dss_phyla_fam <- plot_microbiota_2Fac(
  ps_object = ps_t54_flt,
  exp_group = "gg_group2",
  twoFactor = TRUE,
  fac1 = "diet",
  refFac1 = "50",
  fac2 = "treatment",
  refFac2 = "water",
  sample_name = 'sample_id',
  hues = c("Blues","Greens","Purples","Oranges"), # c("Purples", "Blues", "Reds", "Greens", "Oranges", "Greys", "BuPu")
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
  labs(x = "Mouse ID")
p

# Saving the plot and the associated stats
existingDirCheck("../figures/Thibault_dss/stackbar")
ggsave(plot = p, filename = "../figures/Thibault_dss/stackbar/diet_dss_stackbar.png", width = 7, height = 7, dpi = 300)
writeStackbarExtendedSigTable(main_table = diet_dss_phyla_fam$significant_table_main, includeSubTable = TRUE, sub_table = diet_dss_phyla_fam$significant_table_sub, filepath = "../figures/Thibault_dss/new_filtering/stackbar/diet_dss_stackbar_stats.xlsx")

# pvalues heatmap for the main lvl stats
pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_dss/new_filtering/stackbar/diet_dss_stackbar_stats.xlsx")),
            selected_comparisons = c("50:water_vs_50:dss", "500:water_vs_500:dss","50:water_vs_500:water","50:dss_vs_500:dss"), displayChangeArrows = TRUE, displayPValues = FALSE,
            txn_lvl="Phylum", lvl = "main", taxons = diet_dss_phyla_fam$main_names[!grepl("Others", x = diet_dss_phyla_fam$main_names)], group = "gg_group2", path)

# pvalues heatmap for the sub lvl stats
p = pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_dss/new_filtering/stackbar/diet_dss_stackbar_stats.xlsx")),
                selected_comparisons = c("50:water_vs_50:dss", "500:water_vs_500:dss","50:water_vs_500:water","50:dss_vs_500:dss"),
                txn_lvl="Family", lvl = "sub", taxons =  diet_dss_phyla_fam$sub_names, group = "gg_group2", displayPValues = FALSE, displayChangeArrows = TRUE, path) # You can add [!grepl("Others", x = iron_exp_family$sub_names)] to remove "others"
p+scale_x_discrete(labels = c("50 Ctrl vs 50 DSS", "500 Ctrl vs 500 DSS", "50 Ctrl vs 500 Ctrl", "50 DSS vs 500 DSS"))+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 11, face = "bold"))






sample_data(ps_tfinal_flt)$diet 
sample_data(ps_tfinal_flt)$treatment 
sample_data(ps_tfinal_flt)$gg_group2 <- factor(sample_data(ps_tfinal_flt)$gg_group2, levels = c("50:water", "50:dss", "500:water", "500:dss"))
sample_data(ps_tfinal_flt)$timepoint
length(taxa_sums(ps_tfinal_flt))

# Selected comparisons should be a number of four and follow design as: "
final_phyla_fam <- plot_microbiota_2Fac(
  ps_object = ps_tfinal_flt,
  exp_group = "gg_group2",
  twoFactor = TRUE,
  fac1 = "treatment",
  refFac1 = "water",
  fac2 = "diet",
  refFac2 = "50",
  sample_name = 'sample_id',
  hues = c("Blues","Greens","Purples","Oranges"), # c("Purples", "Blues", "Reds", "Greens", "Oranges", "Greys", "BuPu")
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
  selected_comparisons = list(c("50:water", "500:water"),
                              c("50:dss", "500:dss"),
                              c("50:water", "50:dss"),
                              c("500:water", "500:dss")),
  showOnlySubLegend = FALSE
)

print(final_phyla_fam$plot)
print(final_phyla_fam$significant_table_main)
print(final_phyla_fam$significant_table_sub)

library(ggh4x)

# Custom the plot
p <- final_phyla_fam$plot + 
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
  labs(x = "Mouse ID")
p

# Saving the plot and the associated stats
existingDirCheck("../figures/Thibault_dss/stackbar")
ggsave(plot = p, filename = "../figures/Thibault_dss/stackbar/final_stackbar.png", width = 6, height = 6, dpi = 300)
writeStackbarExtendedSigTable(main_table = final_phyla_fam$significant_table_main, includeSubTable = TRUE, sub_table = final_phyla_fam$significant_table_sub, filepath = "../figures/Thibault_dss/new_filtering/stackbar/final_stackbar_stats.xlsx")

# pvalues heatmap for the main lvl stats
pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_dss/new_filtering/stackbar/final_stackbar_stats.xlsx")),
            selected_comparisons = c("50:water_vs_50:dss", "500:water_vs_500:dss","50:water_vs_500:water","50:dss_vs_500:dss"), displayChangeArrows = FALSE, displayPValues = TRUE,
            txn_lvl="Phylum", lvl = "main", taxons = final_phyla_fam$main_names[!grepl("Others", x = final_phyla_fam$main_names)], group = "gg_group2", path)

# pvalues heatmap for the sub lvl stats
p = pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_dss/new_filtering/stackbar/final_stackbar_stats.xlsx")),
                selected_comparisons = c("50:water_vs_500:water", "50:dss_vs_500:dss","50:water_vs_50:dss","500:water_vs_500:dss"),
                txn_lvl="Family", lvl = "sub", taxons =  final_phyla_fam$sub_names, group = "gg_group2", displayPValues = FALSE, displayChangeArrows = TRUE, path) # You can add [!grepl("Others", x = iron_exp_family$sub_names)] to remove "others"
p+scale_x_discrete(labels = c("50 Ctrl vs 500 Ctrl", "50 DSS vs 500 DSS", "50 Ctrl vs 50 DSS", "500 Ctrl vs 500 DSS"))+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 11))
}

# Firmicutes to bacteroidites ratio
{
  ps_phyla <- tax_glom(ps, taxrank = "Phylum") # Agglom ps at phylum level
  ps_phyla <- transformCounts(ps_phyla, transformation = "rel_ab") # Produce relative abundance table and calculate F/B ratio
  rel_ab <- as.data.frame(otu_table(ps_phyla))
  rel_ab <- rel_ab[rownames(tax_table(ps_phyla)[tax_table(ps_phyla)[,"Phylum"] %in% c("Firmicutes","Bacteroidota"),]),]
  tax_table(ps_phyla)
  rownames(rel_ab) <- c("Firmicutes","Bacteroidota")
  rel_ab <- rel_ab %>%
    tibble::rownames_to_column("phyla") %>%
    pivot_longer(
      cols = -phyla,
      names_to = "sample_id",
      values_to = "abundance"
    ) %>%
    pivot_wider(
      names_from = phyla,
      values_from = abundance,
      names_glue = "{phyla}_abundance"
    ) %>%
    mutate(
      fb_ratio = Firmicutes_abundance / Bacteroidota_abundance
    )
  
  fb_ratio_df <- merge(rel_ab, metadata, by = "sample_id") # Bind metadata information
  fb_ratio_df_diet <- fb_ratio_df[fb_ratio_df$timepoint %in% c("0","35","49"),]
  fb_ratio_df_diet$diet <- factor(fb_ratio_df_diet$diet, levels = c("50","500"), labels = c("50 ppm","500 ppm"))
  fb_ratio_df_diet$week <- factor(fb_ratio_df_diet$week, levels = c("3","8", "10"), labels = c("3 weeks","8 weeks","10 weeks"))
  fbRatioGraphTimeSeries(fb_ratio_df_diet, group = "diet",time = "week", measure = "fb_ratio", custom_colors = c("blue","red"), custom_theme = my_theme())+
    ylim(0,10)+
    geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(0.8),           # left box in each timepoint
      xmax = c(1.2),
      annotations = "n.s.",
      y_position = c(3), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = -0.1,
    )+geom_signif( # For second timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(1.8),           # left box in each timepoint
      xmax = c(2.2),
      annotations = "n.s.",
      y_position = c(8), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+geom_signif( # For third timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(2.8),           # left box in each timepoint
      xmax = c(3.2),
      annotations = "p=0.07",
      y_position = c(9), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = -0.1,
    )
  ggsave("../figures/Thibault_dss/fb_ratio/fb_ratio_diet.png",
         bg = "white",height = 4, width =5, dpi = 300)
  # Stats
  verifyStatsAssumptions(fb_ratio_df_diet[fb_ratio_df_diet$timepoint == "0",], group = "diet",measure =  "fb_ratio")
  wilcox.test(fb_ratio ~ diet, data = fb_ratio_df_diet[fb_ratio_df_diet$timepoint == "0",])
  
  verifyStatsAssumptions(fb_ratio_df_diet[fb_ratio_df_diet$timepoint == "35",], group = "diet",measure =  "fb_ratio")
  t.test(fb_ratio ~ diet, data = fb_ratio_df_diet[fb_ratio_df_diet$timepoint == "35",], var.equal = TRUE)
  
  verifyStatsAssumptions(fb_ratio_df_diet[fb_ratio_df_diet$timepoint == "49",], group = "diet",measure =  "fb_ratio")
  t.test(fb_ratio ~ diet, data = fb_ratio_df_diet[fb_ratio_df_diet$timepoint == "49",], var.equal = FALSE)
  
  # For dss_diet groups
  fb_ratio_df_gg_group2 <- fb_ratio_df[fb_ratio_df$timepoint %in% c("49","54","final"),]
  fb_ratio_df_gg_group2$gg_group2 <- factor(fb_ratio_df_gg_group2$gg_group2, levels = c("50:water","500:water", "50:dss","500:dss"), labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
  fb_ratio_df_gg_group2$week <- factor(fb_ratio_df_gg_group2$week, levels = c("10","10.7", "18"), labels = c("DSS Day 0","DSS Day 5","End of recovery"))
  fb_ratio_df_gg_group2 <- fb_ratio_df_gg_group2 %>%
    mutate(log_FB_ratio = log10(fb_ratio))
  
  fbRatioGraphTimeSeries(fb_ratio_df_gg_group2, group = "gg_group2",time = "week", measure = "fb_ratio", custom_colors = c("blue","red","darkblue", "darkred"), custom_theme = my_theme())+
    ylim(0,9)+
    labs(y = "F/B ratio", title = "F/B ratio overtime\naccording to DSS exposure ")+
    geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(0.7),           # left box in each timepoint
      xmax = c(1.3),
      annotations = "n.s.",
      y_position = c(9), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+
    geom_signif( # For second timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(1.7,2.1,1.7,1.9), # left position
      xmax = c(1.9,2.3,2.1,2.3), # left position
      annotations = c("n.s.","n.s.","***","***"),
      y_position = c(6.9,1.8,7.5,5.6), #
      tip_length = c(0.02,0.02,0.02,0.02),
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+
    geom_signif( # For last timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(2.7,3.1,2.7,2.9),           # left box in each timepoint
      xmax = c(2.9,3.3,3.1,3.3),
      annotations = c("n.s.","n.s.","*","n.s."),
      y_position = c(5.4,6,6.4,6.9), #
      tip_length = c(0.02,0.02,0.02,0.02),
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )
  ggsave("../figures/Thibault_dss/fb_ratio/fb_ratio_dss.png",
         bg = "white",height = 5, width =7, dpi = 300)
  
  # Stats
  verifyStatsAssumptions(fb_ratio_df_gg_group2[fb_ratio_df_gg_group2$timepoint == "49",], group = "gg_group2",measure =  "fb_ratio")
  TukeyHSD(aov(log_FB_ratio ~ gg_group2 , data = fb_ratio_df_gg_group2[fb_ratio_df_gg_group2$timepoint == "49",]))
  
  verifyStatsAssumptions(fb_ratio_df_gg_group2[fb_ratio_df_gg_group2$timepoint == "54",], group = "gg_group2",measure =  "fb_ratio")
  pairwise.wilcox.test(fb_ratio_df_gg_group2[fb_ratio_df_gg_group2$timepoint == "54",]$fb_ratio, fb_ratio_df_gg_group2[fb_ratio_df_gg_group2$timepoint == "54",]$gg_group2, p.adjust.method = "BH")
  
  verifyStatsAssumptions(fb_ratio_df_gg_group2[fb_ratio_df_gg_group2$timepoint == "final",], group = "gg_group2",measure =  "fb_ratio")
  pairwise.wilcox.test(fb_ratio_df_gg_group2[fb_ratio_df_gg_group2$timepoint == "final",]$fb_ratio, fb_ratio_df_gg_group2[fb_ratio_df_gg_group2$timepoint == "final",]$gg_group2, p.adjust.method = "BH")
  
}

# Testing recovery index
{
  library(microbiome)
  
  # For 50 ppm
  # Define baseline centroids
  baseline <- ps %>%
    subset_samples(timepoint == "49" & diet == "50") %>%
    microbiome::transform("compositional")
  
  # Calculate centroid vector
  centroid <- colMeans(t(microbiome::abundances(baseline)))
  
  # Define ps_subset which comprise only 50 ppm samples
  ps_subset <- ps %>%
    subset_samples(diet == "50" & timepoint %in% c("49","54","final"))
  
  # Calculate Bray–Curtis distance to baseline
  ps_comp <- ps_subset %>% microbiome::transform("compositional")
  ab <- microbiome::abundances(ps_comp)
  
  # Pairwise Bray–Curtis from each sample to baseline centroid
  dist_to_base <- apply(ab, 2, function(samp) vegan::vegdist(rbind(centroid, samp), method = "bray")[1])
  
  # Define recovery index = 1 - (BC sample (t49)/max(BC))
  rec_idx <- 1 - dist_to_base / max(dist_to_base)
  sample_data(ps_comp)$RecoveryIndex <- rec_idx
  df_1 <- as.data.frame(sample_data(ps_comp))
  
  
  # For 500 ppm
  # Define baseline centroids
  baseline <- ps %>%
    subset_samples(timepoint == "49" & diet == "500") %>%
    microbiome::transform("compositional")
  
  # Calculate centroid vector
  centroid <- colMeans(t(microbiome::abundances(baseline)))
  
  # Define ps_subset which comprise only 50 ppm samples
  ps_subset <- ps %>%
    subset_samples(diet == "500" & timepoint %in% c("49","54","final"))
  
  # Calculate Bray–Curtis distance to baseline
  ps_comp <- ps_subset %>% microbiome::transform("compositional")
  ab <- microbiome::abundances(ps_comp)
  
  # Pairwise Bray–Curtis from each sample to baseline centroid
  dist_to_base <- apply(ab, 2, function(samp) vegan::vegdist(rbind(centroid, samp), method = "bray")[1])
  
  # Define recovery index = 1 - (BC sample (t49)/max(BC))
  rec_idx <- 1 - dist_to_base / max(dist_to_base)
  sample_data(ps_comp)$RecoveryIndex <- rec_idx
  df_2 <- as.data.frame(sample_data(ps_comp))
  
  # Plot recovery index over time
  df <- rbind(df_1, df_2) # Combine 50/500 data
  df_plot <- df
  df_plot$gg_group2 <- factor(df_plot$gg_group2, levels = c("50:water","500:water","50:dss","500:dss"), labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
  df_plot$timepoint <- factor(df_plot$timepoint, labels = c("DSS day 0","DSS day 5","End of recovery"))
  
  ggplot(df_plot, aes(x = timepoint, y = RecoveryIndex, fill = gg_group2)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.6,
                 color = "black") +
    scale_fill_manual(values = c("blue","red","darkblue", "darkred"))+
    labs(title = "Bray-Curtis recovery index over time", y = "Recovery Index", fill = "Group", x = "Time")+ # , pattern = NA
    my_theme()+
    ylim(0,1)+
    geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(0.7),           # left box in each timepoint
      xmax = c(1.3),
      annotations = "n.s.",
      y_position = c(0.8), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+
    geom_signif( # For second timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(1.7,2.1,1.7,1.9), # left position
      xmax = c(1.9,2.3,2.1,2.3), # left position
      annotations = c("n.s.","n.s.","***","***"),
      y_position = c(0.8,0.45,0.85,0.9), #
      tip_length = c(0.02,0.02,0.02,0.02),
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+
    geom_signif( # For last timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(2.7,3.1,2.7,2.9),           # left box in each timepoint
      xmax = c(2.9,3.3,3.1,3.3),
      annotations = c("n.s.","*","*","n.s."),
      y_position = c(0.8,0.6,0.85,0.9), #
      tip_length = c(0.02,0.02,0.02,0.02),
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )
  ggsave("../figures/Thibault_dss/recovery_index/recovery_index.png",
         bg = "white",height = 5, width =7, dpi = 300)
  
  # Stats
  verifyStatsAssumptions(df[df$timepoint == "49",], group = "gg_group2",measure =  "RecoveryIndex")
  TukeyHSD(aov(RecoveryIndex ~ gg_group2 , data = df[df$timepoint == "49",]))
  
  verifyStatsAssumptions(df[df$timepoint == "54",], group = "gg_group2",measure =  "RecoveryIndex")  
  TukeyHSD(aov(RecoveryIndex ~ gg_group2 , data = df[df$timepoint == "54",]))
  
  verifyStatsAssumptions(df[df$timepoint == "final",], group = "gg_group2",measure =  "RecoveryIndex")  
  TukeyHSD(aov(RecoveryIndex ~ gg_group2 , data = df[df$timepoint == "final",]))
}

# Correlation with iron at t35
{
  library(readxl)
  
  #Preparing dataframe for correlation
  # Load iron measurements data
  df <- as.data.frame(read_xlsx("~/CHUM_git/gut-microbiota-iron/experiments/finished exp/young-DSS-exp3/young48_dss_ferrozine_t35.xlsx"))
  df <- df[,c(2,14)]
  colnames(df) <- df[2,]
  colnames(df)[1:2] <- c("sample_id","iron_concentration")
  df <- df[-c(1,2,27),]
  df$sample_id <- gsub(" ", "_", df$sample_id)
  rownames(df) <- df$sample_id
  df <- df[,-1,drop = FALSE]
  df$iron_concentration <- as.numeric(df$iron_concentration)

  # ps_subset <- tax_glom(ps_subset, taxrank = "Family")
  deseq_subset <- phyloseq_to_deseq2(ps_t35_flt, ~ diet) 
  deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
  res <- results(deseq_subset, name = resultsNames(deseq_subset)[2]) #50 vs 500
  print(resultsNames(deseq_subset))
  
  #For species level
  #One heatmap for all groups
  p <- correlation2Var(ps_t35_flt, deseq_subset, measure = "log2fold", "diet", taxa = "Species", displayPvalue = FALSE, threshold = 0.05, FDR = TRUE, "~/CHUM_git/figures/Thibault_dss/newTaxAnnotation/correlation_heatmaps/", df = df, global = FALSE, singleVariable = TRUE, showIndivCor = FALSE, transformation = "CLR", displayOnlySig = FALSE, returnMainFig = TRUE, displaySpeciesASVNumber = TRUE)
  p <- p+
    labs(title = "Correlation between iron in stools at t35\nand differentially abundant species abundance", y = "", x = "")+
    coord_fixed() + # Makes thing squared
    theme(
      text = element_text(family = "Arial"),      # Global text settings
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),  # Axis titles
      axis.text.y = element_blank(),  # Axis titles
      axis.text.x = element_text(size = 12, face = "bold",hjust = 0.5, vjust = 0.75),   # Axis text
      legend.title = element_text(face = "bold", size = 12),  # Legend title  # Legend text
      axis.ticks = element_blank()
    ) #axis.text.x = element_text(angle = -45, hjust = 1)
  p
  ggsave(filename = "~/CHUM_git/figures/Thibault_dss/correlation_heatmaps/species_iront35_hmap.png",
         plot = p, bg = "white", height = 6, width = 6, dpi = 300)
  
}

# Correlation with DSS-associated parameters
{
  library(readxl)
  
  # Load DSS related results
  variables <- read.csv("~/CHUM_git/gut-microbiota-iron/experiments/finished exp/young-DSS-exp3/young_dss_followup_day5.csv", sep = ",")
  variables <- variables[-22,]
  variables$id <- gsub(x = variables$id, replacement = "", pattern = "[A-Z]")
  variables$id <- paste0(variables$id, "_T54")
  rownames(variables) <- variables$id
  variables <- variables[,c(8,9,13)]
  
  ps_subset <- prune_samples(sample_data(ps_t54_flt)$treatment == "dss", ps_t54_flt)
  ps_subset <- tax_glom(ps_subset, taxrank = "Genus")
  deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ diet) 
  deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
  res <- results(deseq_subset, name = resultsNames(deseq_subset)[2]) #50 vs 500
  print(resultsNames(deseq_subset))
  
  #For species level
  #One heatmap for all groups
  p <- correlation2Var(ps_subset, deseq_subset, measure = "log2fold", "diet", taxa = "Genus", displayPvalue = FALSE, threshold = 0.05, FDR = TRUE, "~/CHUM_git/figures/Thibault_dss/newTaxAnnotation/correlation_heatmaps/", df = variables, global = TRUE, singleVariable = FALSE, showIndivCor = FALSE, transformation = "CLR", displayOnlySig = FALSE, returnMainFig = TRUE, displaySpeciesASVNumber = FALSE)
  p <- p+
    labs(title = "Correlation between iron in stools at t35\nand differentially abundant species abundance", y = "", x = "")+
    coord_fixed() + # Makes thing squared
    theme(
      text = element_text(family = "Arial"),      # Global text settings
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),  # Axis titles
      axis.text.y = element_blank(),  # Axis titles
      axis.text.x = element_text(size = 12, face = "bold",hjust = 0.5, vjust = 0.75),   # Axis text
      legend.title = element_text(face = "bold", size = 12),  # Legend title  # Legend text
      axis.ticks = element_blank()
    ) #axis.text.x = element_text(angle = -45, hjust = 1)
  p
  ggsave(filename = "~/CHUM_git/figures/Thibault_dss/correlation_heatmaps/species_iront35_hmap.png",
         plot = p, bg = "white", height = 6, width = 6, dpi = 300)
}

# Chronobiome 
{
  theme_chronobiome <- function() {
    theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        panel.spacing = unit(0, "lines"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", color = "white", size = 12)
      )
  }
  
  sample_data(ps_flt_diet)$diet <- factor(sample_data(ps_flt_diet)$diet, labels = c("50 ppm","500 ppm"))
  
  p <- plot_timeline_2_groups(
    ps_object = ps_flt_diet,
    exp_group =  "diet", # must be as factor
    time_group = "week", # must be as factor
    sample_name = "sample_id",
    main_level = 'Phylum',
    sub_level = 'Family',
    average_relab_per_group = TRUE,
    smoothing = FALSE,
    hues = c("Blues", "Greens", "Purples", "Oranges"),
    color_bias = 2,
    n_phy = 4,
    custom_theme = theme_chronobiome()
  )
  p+
    facet_wrap2(~ diet, 
                scales  = "free_x", nrow = 2, ncol = 1,
                strip = strip_themed(background_x = elem_list_rect(fill = c("blue", "red"))))+
    labs(x = "Time (weeks)")+
    
    
  existingDirCheck("../figures/Thibault_dss/chronobiome")
  ggsave("../figures/Thibault_dss/chronobiome/diet_chronobiome.png", width = 7, height = 8, dpi = 500, bg = "white")
  
  sample_data(ps_flt_all)$gg_group2 <- factor(sample_data(ps_flt_all)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
  
  p <- plot_timeline_2_groups(
    ps_object = ps_flt_all,
    exp_group =  "gg_group2", # must be as factor
    time_group = "week", # must be as factor
    sample_name = "sample_id",
    main_level = 'Phylum',
    sub_level = 'Family',
    average_relab_per_group = TRUE,
    smoothing = FALSE,
    n_phy = 4,
    hues = c("Blues", "Greens", "Purples", "Oranges"),
    color_bias = 2,
    custom_theme = theme_chronobiome()
  )
  
  p+
    facet_wrap2(~ gg_group2, 
                scales  = "free_x", nrow = 2, ncol = 2,
                strip = strip_themed(background_x = elem_list_rect(fill = c("blue", "red","darkblue","darkred"))))+
    scale_x_continuous(n.breaks = 12)+
    labs(x = "Time (weeks)")
  
  
  ggsave("../figures/Thibault_dss/chronobiome/diet_dss_chronobiome.png", width = 10, height = 8, dpi = 800, bg = "white")
  
  sample_data(ps_flt_all)$gg_group2 <- factor(sample_data(ps_flt_all)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
  
  plot_timeline_taxa(ps_object = ps_flt_all,
                     exp_group =  "gg_group2", # must be as factor
                     time_group = "week", # must be as factor
                     sample_name = "sample_id",
                     main_level = 'Phylum',
                     sub_level = 'Family',
                     average_relab_per_group = TRUE,
                     smoothing = FALSE,
                     threshold = 1,
                     n_phy = 4,
                     hues = c("Blues", "Greens", "Purples", "Oranges"),
                     color_bias = 2,
                     path = "~/CHUM_git/figures/Thibault_dss/chronobiome/",
                     dim = c(9,7),
                     custom_theme = theme_chronobiome(),
                     additionnalAes = list(facet_wrap2(~ gg_group2, 
                                                       scales  = "free_x", nrow = 2, ncol = 2,
                                                       strip = strip_themed(background_x = elem_list_rect(fill = c("blue", "red","darkblue","darkred")))),
                                           labs(x = "Time (weeks)"), scale_x_continuous(n.breaks = 12)))
}

# Produce picrust2 input for whole dataset
producePicrust2Inputs(ps_flt_all, "~/CHUM_git/Microbiota_18_final/")

# Picrust2 - 50 vs 50 dss and 500 vs 500 dss, at t54
{
  meta <- metadata[metadata$diet == "50" & metadata$timepoint == "54",]
  meta <- metadata[metadata$diet == "500" & metadata$timepoint == "54",]
  ids <- meta$id[duplicated(meta$id)]
  ids <- paste0(ids, "_T54")
  meta <- meta[!meta$sample_id %in% ids,]
  meta <- meta[meta$id != "10994",] # For 50 ppm
  meta <- meta[meta$id != "33105",] # For 500 ppm
  
  ko_ab <- read.table("~/CHUM_git/Microbiota_18_final/picrust2/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE) # Load KO annotations
  colnames(ko_ab)[2:ncol(ko_ab)] <- substring(colnames(ko_ab)[2:ncol(ko_ab)], 2)
  pattern <- paste(meta$sample_id, collapse = "|")
  indexes <-  grep(pattern, colnames(ko_ab)) 
  ko_ab <- ko_ab[,c(1,indexes)] # Keep only samples for 10 weeks
  ko_ab <- ko_ab[rowSums(ko_ab[,-1])!=0,]
  kegg_ab <- ko2kegg_abundance(data = ko_ab) # KO to kegg pathways
  
  # Perform differential abundance analysis
  kegg_daa_results_df <- pathway_daa(
    abundance = kegg_ab,
    metadata = meta,
    group = "treatment",
    daa_method = "DESeq2"
  )
  
  # Filter features with p < 0.05
  feature_with_p_0.05 <- kegg_daa_results_df %>%
    filter(p_adjust < 0.00001)
  
  # Retrieve kegg brite hierarchies information
  features <- feature_with_p_0.05$feature
  brite_mapping <- getBriteFromKeggPathID(features)
  
  meta <- meta[,-1] # get rid of id col in metadata
  meta$treatment <- factor(meta$treatment, levels = c("water","dss"), labels = c("50 ppm Ctrl", "50 ppm DSS"))
  meta$treatment <- factor(meta$treatment, levels = c("water","dss"), labels = c("500 ppm Ctrl", "500 ppm DSS"))
  
  
  # custom_col_cat <- terrain.colors(11)
  # custom_col_cat <- heat.colors(11)
  custom_col_cat <- brewer.pal(11, "Set3")
  custom_col_cat <- alpha(custom_col_cat, 0.3)
  
  KeggPathwayHmap(kegg_ab = kegg_ab, brite_mapping = brite_mapping, metadata = meta, group = "treatment",custom_colors_group = c("blue","darkblue"), custom_col_cat, hierarchy = "2")
  KeggPathwayHmap(kegg_ab = kegg_ab, brite_mapping = brite_mapping, metadata = meta, group = "treatment",custom_colors_group = c("red","darkred"), custom_col_cat, hierarchy = "2")
  
  existingDirCheck("~/CHUM_git/figures/Thibault_dss/picrust2")
  ggsave("~/CHUM_git/figures/Thibault_dss/picrust2/kegg_t54_50ppm_hmap.png", bg = "white", height = 13, width = 13, dpi = 300)
  ggsave("~/CHUM_git/figures/Thibault_dss/picrust2/kegg_t54_500ppm_hmap.png", bg = "white", height = 13, width = 13, dpi = 300)
  
}

# Picrust2 - 50 vs 50 dss and 500 vs 500 dss, at tFinal
{
  meta <- metadata[metadata$diet == "50" & metadata$timepoint == "final",]
  meta <- metadata[metadata$diet == "500" & metadata$timepoint == "final",]
  
  ko_ab <- read.table("~/CHUM_git/Microbiota_18_final/picrust2/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE) # Load KO annotations
  colnames(ko_ab)[2:ncol(ko_ab)] <- substring(colnames(ko_ab)[2:ncol(ko_ab)], 2)
  pattern <- paste(meta$sample_id, collapse = "|")
  indexes <-  grep(pattern, colnames(ko_ab)) 
  ko_ab <- ko_ab[,c(1,indexes)] # Keep only samples for 10 weeks
  ko_ab <- ko_ab[rowSums(ko_ab[,-1])!=0,]
  kegg_ab <- ko2kegg_abundance(data = ko_ab) # KO to kegg pathways
  
  # Perform differential abundance analysis
  kegg_daa_results_df <- pathway_daa(
    abundance = kegg_ab,
    metadata = meta,
    group = "treatment",
    daa_method = "DESeq2"
  )
  
  # Filter features with p < 0.05
  feature_with_p_0.05 <- kegg_daa_results_df %>%
    filter(p_adjust < 0.05)
  
  # Retrieve kegg brite hierarchies information
  features <- feature_with_p_0.05$feature
  brite_mapping <- getBriteFromKeggPathID(features)
  
  meta <- meta[,-1] # get rid of id col in metadata
  meta$treatment <- factor(meta$treatment, levels = c("water","dss"), labels = c("50 ppm Ctrl", "50 ppm DSS"))
  meta$treatment <- factor(meta$treatment, levels = c("water","dss"), labels = c("500 ppm Ctrl", "500 ppm DSS"))
  
  
  # custom_col_cat <- terrain.colors(11)
  # custom_col_cat <- heat.colors(11)
  custom_col_cat <- brewer.pal(11, "Set3")
  custom_col_cat <- alpha(custom_col_cat, 0.3)
  
  KeggPathwayHmap(kegg_ab = kegg_ab, brite_mapping = brite_mapping, metadata = meta, group = "treatment",custom_colors_group = c("blue","darkblue"), custom_col_cat, hierarchy = "2")
  KeggPathwayHmap(kegg_ab = kegg_ab, brite_mapping = brite_mapping, metadata = meta, group = "treatment",custom_colors_group = c("red","darkred"), custom_col_cat, hierarchy = "2")
  
  existingDirCheck("~/CHUM_git/figures/Thibault_dss/picrust2")
  ggsave("~/CHUM_git/figures/Thibault_dss/picrust2/kegg_tfinal_50ppm_hmap.png", bg = "white", height = 10, width = 13, dpi = 300)
  ggsave("~/CHUM_git/figures/Thibault_dss/picrust2/kegg_tfinal_500ppm_hmap.png", bg = "white", height = 10, width = 13, dpi = 300)
  
}

# Picrust2 - 50 dss vs 500 dss at t54 / at tfinal
{
  meta <- metadata[metadata$treatment == "dss" & metadata$timepoint == "54",]
  ids <- meta$id[duplicated(meta$id)]
  ids <- paste0(ids, "_T54")
  meta <- meta[!meta$sample_id %in% ids,]
  meta <- meta[!meta$id %in% c("10994","33105"),] # Samples not included
  
  meta <- metadata[metadata$treatment == "dss" & metadata$timepoint == "final",]
  
  ko_ab <- read.table("~/CHUM_git/Microbiota_18_final/picrust2/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE) # Load KO annotations
  colnames(ko_ab)[2:ncol(ko_ab)] <- substring(colnames(ko_ab)[2:ncol(ko_ab)], 2)
  pattern <- paste(meta$sample_id, collapse = "|")
  indexes <-  grep(pattern, colnames(ko_ab)) 
  ko_ab <- ko_ab[,c(1,indexes)] # Keep only samples for 10 weeks
  ko_ab <- ko_ab[rowSums(ko_ab[,-1])!=0,]
  kegg_ab <- ko2kegg_abundance(data = ko_ab) # KO to kegg pathways
  
  # Perform differential abundance analysis
  kegg_daa_results_df <- pathway_daa(
    abundance = kegg_ab,
    metadata = meta,
    group = "diet",
    daa_method = "DESeq2"
  )
  
  # Filter features with p < 0.05 (!!!NON ADJUSTED!!!)
  feature_with_p_0.05 <- kegg_daa_results_df %>%
    filter(p_values < 0.05)
  
  # Retrieve kegg brite hierarchies information
  features <- feature_with_p_0.05$feature
  brite_mapping <- getBriteFromKeggPathID(features)
  
  meta <- meta[,-1] # get rid of id col in metadata
  meta$treatment <- factor(meta$diet, levels = c("50","500"), labels = c("50 ppm DSS", "500 ppm DSS"))
  
  # custom_col_cat <- terrain.colors(11)
  # custom_col_cat <- heat.colors(11)
  custom_col_cat <- brewer.pal(11, "Set3")
  custom_col_cat <- alpha(custom_col_cat, 0.3)
  
  KeggPathwayHmap(kegg_ab = kegg_ab, brite_mapping = brite_mapping, metadata = meta, group = "diet",custom_colors_group = c("darkblue","darkred"), custom_col_cat, hierarchy = "2")
  
  existingDirCheck("~/CHUM_git/figures/Thibault_dss/picrust2")
  ggsave("~/CHUM_git/figures/Thibault_dss/picrust2/kegg_t54_50dss_vs_500dss_hmap.png", bg = "white", height = 4, width = 9, dpi = 300)
  ggsave("~/CHUM_git/figures/Thibault_dss/picrust2/kegg_tfinal_50dss_vs_500dss_hmap.png", bg = "white", height = 4, width = 9, dpi = 300)
  
}

# Creating a cauliflower phylogenetic tree
{
  library(ggtree)
  library(ape)
  library(metacoder)
  library(RColorBrewer)
  
  tax_table(ps_tree)[is.na(tax_table(ps_tree))] <- "Unknown"
  
  tax_data <- parse_phyloseq(ps_tree)
  all_ids  <- tax_data$taxon_ids()
  tip_df   <- tax_data$data$tax_data
  
  taxa_df <- data.frame(taxon_id = all_ids, Phylum = NA_character_, stringsAsFactors = FALSE)
  taxa_df$Phylum[match(tip_df$taxon_id, all_ids)] <- as.character(tip_df$Phylum)
  tax_data$data$taxa <- taxa_df
  
  el <- tax_data$edge_list
  
  # Initialize propagated Phylum
  prop_phylum <- setNames(taxa_df$Phylum, taxa_df$taxon_id)
  
  # Fill internal nodes by inheriting from their descendants
  repeat {
    # Find internals missing Phylum
    no_phylum_ids <- taxa_df$taxon_id[is.na(prop_phylum)]
    updated <- FALSE
    
    for (nid in no_phylum_ids) {
      kids <- el$to[el$from == nid]
      ph_values <- unique(na.omit(prop_phylum[kids]))
      if (length(ph_values) > 0) {
        prop_phylum[nid] <- ph_values[1]
        updated <- TRUE
      }
    }
    if (!updated) break
  }
  
  tax_data$data$taxa$Phylum <- prop_phylum
  
  unique_phyla <- sort(unique(prop_phylum))
  # palette_vec <- colorRampPalette(brewer.pal(min(8, length(unique_phyla)), "Set3"))(length(unique_phyla))
  palette_vec <- c("#483C46","#70AE6E","#FFC300","#900C3F","#C70039","#3C6E71","#BEEE62","#FF5733")
  palette_vec <- c("#c7522a","#e5c185","#f0daa5","#fbf2c4","#b8cdab","#74a892","#008585","#004343")
  palette_vec <- c("#003f5c","#58508d","#8a508f","#bc5090","#de5a79","#ff6361","#ff8531","#ffa600")
  palette_vec <- c()
  palette_vec <- c("#715660","#846470","#566071","#727f96","#c39f72","#d5bf95","#667762","#879e82","gray")
  
  phylum_pal  <- setNames(palette_vec, unique_phyla)
  # Shuffle (randomize) the values, keeping the same names
  phylum_pal[] <- sample(phylum_pal)
  tax_data$data$taxa$color <- phylum_pal[prop_phylum]
  
  heat_tree(
    tax_data,
    node_color       = color,
    edge_color       = color,
    node_label       = "",
    node_size_trans = "area",
    node_size = n_obs,
    node_color_range = phylum_pal,
    edge_color_range = phylum_pal,
    make_node_legend = FALSE,
    node_legend_title = "Phylum",
    layout           = "davidson-harel",
    initial_layout = "reingold-tilford"
    # initial_layout = "fruchterman-reingold"
  )
  
  
   phylum_palibrary(patchwork)  # For combining plots
  
  # Data frame for legend
  legend_df <- data.frame(
    Phylum = names(phylum_pal),
    Color  = phylum_pal
  )
  
  # Make dummy points for legend
  p_legend <- ggplot(legend_df, aes(x = 1, y = Phylum, color = Phylum)) +
    geom_point(size = 3) +
    scale_color_manual(values = phylum_pal) +
    theme_void() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 10, hjust = 0),
      plot.margin = margin(2, 2, 2, 2)
    ) +
    labs(title = "Phylum")
  
  final_plot <- p_tree + p_legend + plot_layout(widths = c(5, 1))
  final_plot
  
  ggsave(filename = "~/CHUM_git/figures/Thibault_iron/cauliflower phytree/phy_tree.png", plot = final_plot, bg = "white", height = 7, width = 11, dpi = 300)
}




# Prepare picrust2 input for last timepoint
ps_subset <- prune_samples(sample_data(ps_dss_relab_flt)$timepoint == "final", ps_dss_relab_flt)

# Filtering
# Function filtering out ASVs for which they were in total less than a threshold count
ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
# Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
producePicrust2Inputs(ps_subset, "~/CHUM_git/Microbiota_18/")




# Prepare picrust2 input for t35
ps_subset <- prune_samples(sample_data(ps_flt_diet)$timepoint == "35", ps_flt_diet)

# Filtering
# Function filtering out ASVs for which they were in total less than a threshold count
ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
# Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
existingDirCheck("~/CHUM_git/Microbiota_18/t35")
producePicrust2Inputs(ps_subset, "~/CHUM_git/Microbiota_18/t35/")




# Prepare picrust2 input for t49
ps_subset <- prune_samples(sample_data(ps_dss_relab_flt)$timepoint == "49", ps_dss_relab_flt)

# Filtering
# Function filtering out ASVs for which they were in total less than a threshold count
ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
# Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
existingDirCheck("~/CHUM_git/Microbiota_18/t49")
producePicrust2Inputs(ps_subset, "~/CHUM_git/Microbiota_18/t49/")





# Prepare picrust2 input for t49
ps_subset <- prune_samples(sample_data(ps_dss_relab_flt)$timepoint == "54", ps_dss_relab_flt)

# Filtering
# Function filtering out ASVs for which they were in total less than a threshold count
ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
# Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
existingDirCheck("~/CHUM_git/Microbiota_18/t54")
producePicrust2Inputs(ps_subset, "~/CHUM_git/Microbiota_18/t54/")


























































































# New tests for beta diversity

# Beta diversity (for only diets timepoints)
# Bray curtis filtered (good results)
betaDiversityTimepoint2Factors(ps_flt_diet, sample_id = "sample_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="bray",
                               transform = "rel_ab", customColors = c("blue","red"),
                               font = "Arial", path = "../figures/Thibault_dss/diet/beta_diversity/filtered/")

# Beta diversity (for diet + treatment)
# Bray curtis filtered (good for t54, less interpretable for tfinal)
betaDiversityTimepoint2Factors(ps_flt_dss, sample_id = "sample_id", timeVariable = "timepoint",
                               varToCompare =  "gg_group2", distMethod ="bray",
                               customColors = c("blue","red","blue4","red4"),
                               font = "Arial", path = "../figures/Thibault_dss/diet_dss/beta_diversity/filtered/")

# Beta diversity, for control only at tfinal
ps_subset <- prune_samples(sample_data(ps_flt_dss)$timepoint == "final" & sample_data(ps_flt_dss)$treatment == "water", ps_flt_dss) # Create subset
# Bray curtis filtered 
betaDiversityTimepoint2Factors(ps_subset, sample_id = "sample_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="bray",
                               transform = "rel_ab", customColors = c("blue","red"),
                               font = "Arial", path = "../figures/Thibault_dss/diet/beta_diversity/filtered/")

# Beta diversity, for dss only at t54 and tfinal
ps_subset <- prune_samples(sample_data(ps_flt_dss)$treatment == "dss", ps_flt_dss) # Create subset
# Bray curtis filtered (results show that they are more different at the end of the recovery than at last timepoint of dss)
betaDiversityTimepoint2Factors(ps_subset, sample_id = "sample_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="bray",
                               transform = "rel_ab", customColors = c("darkblue","darkred"),
                               font = "Arial", path = "../figures/Thibault_dss/dss/beta_diversity/filtered/")

# Beta diversity, for dss only at t54 and tfinal, test with Jaccard index 
ps_subset <- prune_samples(sample_data(ps_flt_dss)$treatment == "dss", ps_flt_dss) # Create subset
length(taxa_sums(ps))
ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)
length(taxa_sums(ps))
betaDiversityTimepoint2Factors(ps_subset, sample_id = "sample_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="jaccard",
                               transform = "none", customColors = c("darkblue","darkred"),
                               font = "Arial", path = "../figures/Thibault_dss/dss/beta_diversity/jaccard/")

# Beta diversity, for 50 only at t54 and tfinal (showing that they remain different)
ps_subset <- prune_samples(sample_data(ps_dss)$diet == "50", ps_dss) # Create subset
# Bray curtis filtered (results show that they are more different at the end of the recovery than at last timepoint of dss)
betaDiversityTimepoint2Factors(ps_subset, sample_id = "sample_id", timeVariable = "timepoint",
                               varToCompare =  "treatment", distMethod ="bray",
                               transform = "rel_ab", customColors = c("blue","darkblue"),
                               font = "Arial", path = "../figures/Thibault_dss/diet_50/beta_diversity/filtered/")

# Beta diversity, for 500 only at t54 and tfinal (showing that they remain different)
ps_subset <- prune_samples(sample_data(ps_flt_dss)$diet == "500", ps_flt_dss) # Create subset
# Bray curtis filtered (results show that they are more different at the end of the recovery than at last timepoint of dss)
betaDiversityTimepoint2Factors(ps_subset, sample_id = "sample_id", timeVariable = "timepoint",
                               varToCompare =  "treatment", distMethod ="bray",
                               transform = "rel_ab", customColors = c("red","darkred"),
                               font = "Arial", path = "../figures/Thibault_dss/diet_500/beta_diversity/filtered/")


# Now we apply dbRDA for more complex designs with all groups for t54 and tfinal

ps_subset <- prune_samples(sample_data(ps_flt_dss)$timepoint == "54", ps_flt_dss)
betaDiversityTimepointsGroupedDbRDA(ps = ps_subset, sample_id = "sample_id", varToCompare = "gg_group2",distMethod = "hellinger",
                                    customColors = c("blue","red","blue4","red4"), formula = "treatment * diet + cage", transform = "rel_ab",
                                    path = "../figures/Thibault_dss/diet_dss/dbRDA_t54/", font = "Arial", additionnalAes = NULL)

ps_subset <- prune_samples(sample_data(ps_flt_dss)$timepoint == "final", ps_flt_dss)
betaDiversityTimepointsGroupedDbRDA(ps = ps_subset, sample_id = "sample_id", varToCompare = "gg_group2",distMethod = "hellinger",
                                    customColors = c("blue","red","blue4","red4"), formula = "treatment * diet + cage", transform = "rel_ab",
                                    path = "../figures/Thibault_dss/diet_dss/dbRDA_final/", font = "Arial", additionnalAes = NULL)


















# Testing a presence/absence approach on the tfinal for the diets
ps_subset <- prune_samples(sample_data(ps)$timepoint %in% c("final") & sample_data(ps)$treatment == "dss", ps)
ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.2 * nsamples(ps_subset)), ps_subset)
length(taxa_sums(ps_subset))
ps_subset <- subset_taxa(ps_subset, !is.na(tax_table(ps_subset)[, "Species"])) # Keep ASVs that were identified at the species level
print(length(taxa_sums(ps_subset)))
# jaccard_dist <- phyloseq::distance(ps_subset, method = "jaccard", binary = TRUE)
otu_table <- as.data.frame(t(otu_table(ps_subset)))
presence_absence <- as.data.frame(otu_table > 0) # Transform into matrix of presence absence
presence_absence[] <- lapply(presence_absence, as.numeric) # Ensure data is numeric and not boolean
presence_absence <- as.matrix(presence_absence)

# Subset of mice from 50:dss group
dss50 <- presence_absence[,sample_names(ps_subset)[sample_data(ps_subset)$diet == "50"]]
dss500 <- presence_absence[,sample_names(ps_subset)[sample_data(ps_subset)$diet == "500"]]


asvList50 <- row.names(dss50[!rowSums(dss50) > 0,])
asvList500 <- row.names(dss500[!rowSums(dss500) > 0,])
# Identify ASVs that are only present for dss50 and dss500
# # Remove columns (taxa) with all zeros (or constant values)
# presence_absence <- presence_absence[, apply(presence_absence, 2, function(x) length(unique(x)) > 1), drop = FALSE]


if (any(is.na(rownames(presence_absence))) || any(duplicated(rownames(presence_absence)))) {
  stop("There are missing or duplicated sample names in your presence_absence matrix.")
}

# Check for duplicate or missing column names (taxa)
if (any(is.na(colnames(presence_absence))) || any(duplicated(colnames(presence_absence)))) {
  stop("There are missing or duplicated taxa names in your presence_absence matrix.")
}

library(vegan)


# Compute Jaccard dissimilarity matrix
jaccard_dist <- vegdist(presence_absence, method = "jaccard", binary = TRUE)

# Perform NMDS
nmds <- metaMDS(jaccard_dist, k = 2, trymax = 100)

# Plot NMDS
plot(nmds, type = "t")
points(nmds, col = as.factor(sample_data$DietGroup), pch = 16)































































#### Tests for timeline graphs ####
source("~/CHUM_git/gut-microbiota-iron/pipelinelinux/microbiota_analysis/chronobiome.R")

plot_timeline_2_groups(
  ps_object = ps_flt_diet,
  exp_group =  "diet", # must be as factor
  time_group = "week", # must be as factor
  sample_name = "sample_id",
  main_level = 'Phylum',
  sub_level = 'Family',
  average_relab_per_group = TRUE,
  n_phy = 4,
  differential_analysis = FALSE,
  test = c("Wald", "LRT")[1],
  fdr_threshold = 0.05,
  sig_lab = FALSE,
  fitType = c("parametric", "local", "mean", "glmGamPoi")[1],
  sfType = c("ratio", "poscounts", "iterate")[1],
  betaPrior = FALSE,
  reduced = FALSE,
  quiet = TRUE,
  minReplicatesForReplace = 7,
  modelMatrixType = c("standard", "expanded")[1],
  useT = FALSE,
  minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5,
  parallel = FALSE
  
)



# Run statistical tests, for now, independent from the chronobiome functions
# LRT to detect time effects within each group
# Subset for 50 ppm, and for family level taxa
ps_subset <- prune_samples(sample_data(ps_flt_diet)$diet == "50", ps_flt_diet)
ps_taxa <- tax_glom(ps_subset, taxrank = "Phylum")
deseq_subset <- phyloseq_to_deseq2(ps_taxa, ~ week)
deseq_subset <- DESeq(deseq_subset, test="LRT", reduced = ~1)


# deseq_subset <- phyloseq_to_deseq2(ps_taxa, ~ diet+week+diet:week)
# deseq_subset <- DESeq(deseq_subset, test="LRT", reduced = ~diet+week)
resultsNames(deseq_subset)
print(results(deseq_subset))
res <- results(deseq_subset)
print(res)
res <- results(deseq_subset, contrast = c("week", "3", "8"), test = "Wald")
print(res)
res <- results(deseq_subset, contrast = c("week", "8", "10"), test = "Wald")
print(res)

sigtab <- NULL
res_sub <- results(deseq_subset, contrast = c("week", , cmp[[2]]))
res_sub <- subset(res, padj < 0.05) # Keep only significant fematures
res_sub <- cbind(as(res_sub, "data.frame"), as(tax_table(ps_taxa)[rownames(res_sub), ], "matrix"))
res_sub$comparaison <- paste0(cmp[[1]],"_vs_",cmp[[2]])
sigtab <- bind_rows(sigtab, res_sub) # Append res_sub to final complete sigtab 


selected_comparisons = list(c("3", "8"),
                            c("8", "10"))

sigtab <- NULL
for(cmp in selected_comparisons){
  res_sub <- results(deseq_subset, contrast = c("week", cmp[[1]], cmp[[2]]))
  res_sub <- subset(res_sub, padj < 0.05) # Keep only significant features
  res_sub <- cbind(as(res_sub, "data.frame"), as(tax_table(ps_taxa)[rownames(res_sub), ], "matrix"))
  res_sub$comparaison <- paste0(cmp[[1]],"_vs_",cmp[[2]])
  sigtab <- bind_rows(sigtab, res_sub) # Append res_sub to final complete sigtab 
}
sigtab$ASV <- gsub("^(ASV\\d{1,3}).*", "\\1", rownames(sigtab))
View(sigtab)

if (nrow(significant_features) == 0) {
  message("No significant feature detected at")
} else {
  significant_features <- as.data.frame(cbind(significant_features,
                                              tax_table(fam_glom)[rownames(tax_table(fam_glom)) %in%
                                                                    rownames(significant_features), ]))
  significant_features_sub <- significant_features
  df_long$differential_abundance <- FALSE
  df_long$differential_abundance[df_long$plot_taxa %in% significant_features[, sub_level]] <- TRUE
  significant_features$stars <- ""
  if (sig_lab == TRUE) {
    significant_features$stars <- symnum(significant_features$padj,
                                         symbols = c("***", "**", "*", ""),
                                         cutpoints = c(0, .001, .01, .05, 1),
                                         corr = FALSE)
    star_vec <- significant_features$stars[match(df_long$plot_taxa, significant_features[, sub_level])]
    star_vec[is.na(star_vec)] <- ""
    df_long$plot_taxa <- paste0(df_long$plot_taxa, " ", star_vec)
  }
  df_long$legend_label <- ifelse(df_long$differential_abundance,
                                 paste0("<b>", df_long$plot_taxa, "</b>"),
                                 as.character(df_long$plot_taxa))
  df_long$plot_taxa <- df_long$legend_label
  df_long$plot_taxa <- factor(df_long$plot_taxa, levels = unique(df_long$plot_taxa))
}







# LRT analysis
# Subset for 50 ppm
relabTimepoints(ps_flt_diet, deseq_subset, timeVariable = "week", varToCompare = "diet", taxa = "Species", threshold == 0.05, customColors = c("blue","red"), path = "~/CHUM_git/figures/Thibault_dss/lrt_diet/")
# Subset for 50 ppm, and for family level taxa
ps_subset <- prune_samples(sample_data(ps_flt_diet)$diet == "50", ps_flt_diet)
ps_taxa <- tax_glom(ps_subset, taxrank = "Phylum")
deseq_subset <- phyloseq_to_deseq2(ps_taxa, ~ week)
deseq_subset <- DESeq(deseq_subset, test="LRT", reduced = ~1)