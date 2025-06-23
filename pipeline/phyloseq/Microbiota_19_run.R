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
  
}

# Load custom functions for microbiota analysis
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/utilities.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/alpha_diversity_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/beta_diversity_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/correlation_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/relab_analysis_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/taxa_distrib_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/plot_microbiota_extension.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/deseq2_log2fold_change_analysis.R")

# For microbiota 19
#set working directory
setwd("~/Documents/CHUM_git/Microbiota_19/")
asv_table <- as.data.frame(fread("asv_table/asv_table.csv", sep = ";"))
rownames(asv_table) <- asv_table[,1]  # Use the first column as row names
asv_table <- asv_table[,-1]  # Drop the first column
rownames(asv_table) <- gsub("_16S", "", rownames(asv_table))

# Metadata handling
{
  #loading metadata of interest
  metadata <- read.xlsx("metadata/dissection.xlsx")
  metadata <- metadata[,-c(5:8)] # Remove the non-metadata stuff (liver measures and stuff)
  metadata$ID <- substring(metadata$ID, 1, 5) # Remove the letter at the end of id
  rownames(metadata) <- metadata$id #adding id col as rownames too
  metadata <- metadata[-17,] # Remove dead mouse
  colnames(metadata)[2] <- "id"
  
  # Extract 16S reads sample ids
  samples <- read.xlsx("metadata/Microbiota_19_NextSeqReadSet_2025-06-22.xlsx")
  samples <- as.data.frame(samples$Nom)
  colnames(samples) <- "sample_id"
  samples$id <- substring(samples$sample_id, 1, 5)
  samples$timepoint <- substring(samples$sample_id, 8, nchar(samples$sample_id)) 
  samples <- samples[-129,] # Remove pcrBlank sample
  
  # Bind both metadata df to link timepoints with their metadata (diet and treatment)
  metadata <- merge(samples, metadata, by = "id")
  metadata$week <- ifelse(metadata$timepoint == "final", "18", as.character(round(as.numeric(metadata$timepoint)/7, 1)+3)) # Add week column 
  metadata$gg_group <- # Adding gg_group variable (combination of time, diet and treatment)
    paste(metadata$timepoint, 
          metadata$diet,
          metadata$treatment, 
          sep = ":")
  metadata$gg_group2 <- # Another gg_group variable (diet and treatment only)
    paste(metadata$diet,
          metadata$treatment, 
          sep = ":")
  rownames(metadata) <- metadata$sample_id # Put full_id as rownames
}

# Load taxonomical assignments
taxa <- as.matrix(fread("taxonomy/taxa_annotation.csv", sep = ";"))
rownames(taxa) <- taxa[,1]  # Use the first column as row names
taxa <- taxa[,-1]  # Drop the first column

# Load phylogenetic tree if possible
tree <- read.tree("~/Documents/CHUM_git/Microbiota_18/taxonomy/phylogenetic_tree.newick")

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

sum(taxa_sums(ps)) # total number of reads
length(taxa_sums(ps)) # total number of ASVs
nrow(tax_table(ps))-sum(is.na(tax_table(ps)[,7])) # how many identified at species level
nrow(tax_table(ps))-sum(is.na(tax_table(ps)[,6])) # how many identified at genus level
nrow(tax_table(ps))-sum(is.na(tax_table(ps)[,5])) # how many identified at family level
nrow(tax_table(ps))-sum(is.na(tax_table(ps)[,4])) # how many identified at order level
nrow(tax_table(ps))-sum(is.na(tax_table(ps)[,3])) # how many identified at class level

# Put as factors variables that are going to be used
sample_data(ps)$gg_group2 <- factor(sample_data(ps)$gg_group2, levels = c("50:water", "500:water", "50:abx", "500:abx")) # Put gg_group2 as factor
sample_data(ps)$timepoint <- factor(sample_data(ps)$timepoint, levels = c("0","35","49","56","final")) # Put timepoint as factor
sample_data(ps)$week <- factor(sample_data(ps)$week, levels = c("3","8","10","11","18")) # Put week as factor
sample_data(ps)$treatment <- factor(sample_data(ps)$treatment, levels = c("water","abx")) # Put treatment as factor
sample_data(ps)$diet <- factor(sample_data(ps)$diet, levels = c("50","500")) # Put diet as factor

# Create single timepoint phyloseq objects and apply filter
ps_t0 <- prune_samples(sample_data(ps)$timepoint %in% c("0"), ps)
ps_t0_flt <- prune_taxa(taxa_sums(ps_t0) > 10, ps_t0)
ps_t0_flt <- prune_taxa(colSums(otu_table(ps_t0_flt) > 0) >= (0.3 * nsamples(ps_t0_flt)), ps_t0_flt)
length(taxa_sums(ps_t0_flt))

ps_t35 <- prune_samples(sample_data(ps)$timepoint %in% c("35"), ps)
ps_t35_flt <- prune_taxa(taxa_sums(ps_t35) > 10, ps_t35)
ps_t35_flt <- prune_taxa(colSums(otu_table(ps_t35_flt) > 0) >= (0.3 * nsamples(ps_t35_flt)), ps_t35_flt)
length(taxa_sums(ps_t35_flt))

ps_t49 <- prune_samples(sample_data(ps)$timepoint %in% c("49"), ps)
ps_t49_flt <- prune_taxa(taxa_sums(ps_t49) > 10, ps_t49)
ps_t49_flt <- prune_taxa(colSums(otu_table(ps_t49_flt) > 0) >= (0.3 * nsamples(ps_t49_flt)), ps_t49_flt)
length(taxa_sums(ps_t49_flt))

ps_t56 <- prune_samples(sample_data(ps)$timepoint %in% c("56"), ps)
ps_t56_flt <- prune_taxa(taxa_sums(ps_t56) > 10, ps_t56)
ps_t56_flt <- prune_taxa(colSums(otu_table(ps_t56_flt) > 0) >= (0.3 * nsamples(ps_t56_flt)), ps_t56_flt)
length(taxa_sums(ps_t56_flt))

ps_tfinal <- prune_samples(sample_data(ps)$timepoint %in% c("final"), ps)
ps_tfinal_flt <- prune_taxa(taxa_sums(ps_tfinal) > 10, ps_tfinal)
ps_tfinal_flt <- prune_taxa(colSums(otu_table(ps_tfinal_flt) > 0) >= (0.3 * nsamples(ps_tfinal_flt)), ps_tfinal_flt)
length(taxa_sums(ps_tfinal_flt))

# Create phyloseq obejcts that we need for the analysis
ps_diet <- merge_phyloseq(ps_t0, ps_t35, ps_t49)
ps_abx_alpha <- merge_phyloseq(ps_t49, ps_t56, ps_tfinal)
ps_abx_relab_flt <- merge_phyloseq(ps_t49_flt, ps_t56_flt, ps_tfinal_flt)
ps_abx <- merge_phyloseq(ps_t56, ps_tfinal)
ps_flt_diet <- merge_phyloseq(ps_t0_flt, ps_t35_flt, ps_t49_flt)
ps_flt_abx <- merge_phyloseq(ps_t56_flt, ps_tfinal_flt)
ps_flt_all <- merge_phyloseq(ps_t0_flt, ps_t35_flt, ps_t49_flt, ps_t56_flt, ps_tfinal_flt)

# Alpha diveristy
{
  existingDirCheck("../figures/Thibault_abx/diet/")
  graphs = alphaDiversityTimeSeries2(ps_diet, "../figures/Thibault_abx/diet/", time = "timepoint", group = "diet", writeData = TRUE)
  
  # Stats
  alpha_d <- read.xlsx("../figures/Thibault_abx/diet/alpha_diversity/alpha_diversity_data.xlsx")
  alpha_d$diet <- factor(alpha_d$diet, levels = c("50", "500"))
  alpha_d$timepoint <- factor(alpha_d$timepoint, levels = c("0", "35", "49"))
  
  # Chao1
  graphs[[1]]+
    scale_fill_manual(values = c("blue","red"),
                      labels = c("50 ppm","500 ppm"))+
    scale_x_discrete(labels = c("3 weeks\n(Weaning)","5 weeks", "7 weeks"))+
    labs(y = "Chao1 Index", x = "")+
    my_theme()+
    ylim(0,1200)+
    geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(0.8),           # left box in each timepoint
      xmax = c(1.2),
      annotations = "n.s.",
      y_position = c(1150), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(1.8),           # left box in each timepoint
      xmax = c(2.2),
      annotations = "n.s.",
      y_position = c(1150), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(2.8),           # left box in each timepoint
      xmax = c(3.2),
      annotations = "**",
      y_position = c(1150), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )
  ggsave("../figures/Thibault_abx/diet/alpha_diversity/chao1.png",
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
    ylim(0,5.5)+
    geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(0.8),           # left box in each timepoint
      xmax = c(1.2),
      annotations = "n.s.",
      y_position = c(5), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(1.8),           # left box in each timepoint
      xmax = c(2.2),
      annotations = "p=0.07",
      y_position = c(5), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 2.5,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(2.8),           # left box in each timepoint
      xmax = c(3.2),
      annotations = "n.s.",
      y_position = c(5), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )
  ggsave("../figures/Thibault_abx/diet/alpha_diversity/shannon.png",
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
    ylim(0,62)+
    geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(0.8),           # left box in each timepoint
      xmax = c(1.2),
      annotations = "n.s.",
      y_position = c(60), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(1.8),           # left box in each timepoint
      xmax = c(2.2),
      annotations = "p=0.06",
      y_position = c(36), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 2.5,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(2.8),           # left box in each timepoint
      xmax = c(3.2),
      annotations = "n.s.",
      y_position = c(36), #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )
  ggsave("../figures/Thibault_abx/diet/alpha_diversity/invsimpson.png",
         bg = "white",height = 4, width =5, dpi = 300)
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "0",], group = "diet", measure = "InvSimpson")
  wilcox.test(InvSimpson ~ diet, data = alpha_d[alpha_d$timepoint == "0",])
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "35",], group = "diet", measure = "InvSimpson")
  wilcox.test(InvSimpson ~ diet, data = alpha_d[alpha_d$timepoint == "35",])
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "diet", measure = "InvSimpson")
  wilcox.test(InvSimpson ~ diet, data = alpha_d[alpha_d$timepoint == "49",])
  
  # Alpha diveristy for dss_diet only
  existingDirCheck("../figures/Thibault_abx/diet_dss/")
  graphs = alphaDiversityTimeSeries2(ps_abx_alpha, "../figures/Thibault_abx/diet_dss/", time = "timepoint", group = "gg_group2", writeData = TRUE)
  
  # Load data for stats calculations
  alpha_d <- read.xlsx("../figures/Thibault_abx/diet_abx/alpha_diversity/alpha_diversity_data.xlsx")
  alpha_d$gg_group2 <- factor(alpha_d$gg_group2, levels = c("50:water","500:water","50:abx","500:abx"))
  alpha_d$timepoint <- factor(alpha_d$timepoint, levels = c("49", "56", "final"))
  
  # Chao1
  graphs[[1]]+
    scale_fill_manual(values = c("blue","red","deepskyblue", "brown1"),
                      labels = c("50 ppm control","500 ppm control","50 ppm Abx","500 ppm Abx"))+
    scale_pattern_manual(values = c("circle","stripe","circle","stripe"))+
    scale_x_discrete(labels = c("Abx day 0","Abx day 7", "End of\nrecovery"),
                     expand = c(0,0.5))+
    labs(y = "Species Richness", x = "")+
    guides(pattern = "none")+
    ylim(0,1350)+
    my_theme()+
    geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(0.7),           # left box in each timepoint
      xmax = c(1.3),
      annotations = "n.s.",
      y_position = c(1120), #
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
      annotations = c("n.s.","*","***","***"),
      y_position = c(1100,600,1200,1300), #
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
        annotations = c("n.s.","n.s.","**","**"),
        y_position = c(1120,800,1200,1300), #
        tip_length = c(0.02,0.02,0.02,0.02),
        color = "black",
        size = 0.5,
        textsize = 4,
        margin_top = 0.1, # Moves the top according to this value
        vjust = 0,
      )
    
  ggsave("../figures/Thibault_abx/diet_abx/alpha_diversity/chao1.png",
         bg = "white",height = 5, width =7, dpi = 300)
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "gg_group2", measure = "Chao1")
  TukeyHSD(aov(Chao1 ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "49",]))
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "56",], group = "gg_group2", measure = "Chao1")
  TukeyHSD(aov(Chao1 ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "56",]))
  pairwise.wilcox.test(alpha_d[alpha_d$timepoint == "56",]$Chao1, alpha_d[alpha_d$timepoint == "56",]$gg_group2, p.adjust.method = "BH")
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "final",], group = "gg_group2", measure = "Chao1")
  TukeyHSD(aov(Chao1 ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "final",]))

  
  # Shannon
  graphs[[2]]+
    scale_fill_manual(values = c("blue","red","deepskyblue", "brown1"),
                      labels = c("50 ppm control","500 ppm control","50 ppm Abx","500 ppm Abx"))+
    scale_pattern_manual(values = c("circle","stripe","circle","stripe"))+
    scale_x_discrete(labels = c("Abx day 0","Abx day 7", "End of\nrecovery"),
                     expand = c(0,0.5))+
    labs(y = "Shannon Index", x = "", title = "Shannon index")+
    guides(pattern = "none")+
    ylim(0,5)+
    my_theme()+
    geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(0.7),           # left box in each timepoint
      xmax = c(1.3),
      annotations = "n.s.",
      y_position = c(4.2), #
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
      annotations = c("n.s.","*","***","***"),
      y_position = c(4.2,3.5,4.5,4.9), #
      tip_length = c(0.02,0.02,0.02,0.02),
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+
      geom_signif( # For last timepoint
        # comparisons = list(c(groups[1],groups[4])),
        xmin = c(2.7),           # left box in each timepoint
        xmax = c(3.3),
        annotations = "n.s.",
        y_position = c(4.2), #
        tip_length = 0,
        color = "black",
        size = 0.5,
        textsize = 4,
        margin_top = 0.1, # Moves the top according to this value
        vjust = 0,
      )
  ggsave("../figures/Thibault_abx/diet_abx/alpha_diversity/shannon.png",
         bg = "white",height = 5, width =7, dpi = 300)
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "gg_group2", measure = "Shannon")
  TukeyHSD(aov(Shannon ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "49",]))
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "56",], group = "gg_group2", measure = "Shannon")
  pairwise.wilcox.test(alpha_d[alpha_d$timepoint == "56",]$Shannon, alpha_d[alpha_d$timepoint == "56",]$gg_group2, p.adjust.method = "BH")
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "final",], group = "gg_group2", measure = "Shannon")
  TukeyHSD(aov(Shannon ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "final",]))
  
  # InvSimpson
  graphs[[3]]+
    scale_fill_manual(values = c("blue","red","deepskyblue", "brown1"),
                      labels = c("50 ppm control","500 ppm control","50 ppm Abx","500 ppm Abx"))+
    scale_pattern_manual(values = c("circle","stripe","circle","stripe"))+
    scale_x_discrete(labels = c("Abx day 0","Abx day 7", "End of\nrecovery"),
                     expand = c(0,0.5))+
    labs(y = "Inverse Simpson Index", x = "", title = "Inverse Simpson index")+
    guides(pattern = "none")+
    ylim(0,22)+
    my_theme()+
    geom_signif( # For first timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(0.7),           # left box in each timepoint
      xmax = c(1.3),
      annotations = "n.s.",
      y_position = c(13), #
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
      annotations = c("n.s.","*","***","**"),
      y_position = c(16.5,11,18,19.5), #
      tip_length = c(0.02,0.02,0.02,0.02),
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )+
    geom_signif( # For last timepoint
      # comparisons = list(c(groups[1],groups[4])),
      xmin = c(2.7),           # left box in each timepoint
      xmax = c(3.3),
      annotations = "n.s.",
      y_position = 13, #
      tip_length = 0,
      color = "black",
      size = 0.5,
      textsize = 4,
      margin_top = 0.1, # Moves the top according to this value
      vjust = 0,
    )
  ggsave("../figures/Thibault_abx/diet_abx/alpha_diversity/InvSimpson.png",
         bg = "white",height = 5, width =7, dpi = 300)
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "gg_group2", measure = "InvSimpson")
  pairwise.wilcox.test(alpha_d[alpha_d$timepoint == "49",]$InvSimpson, alpha_d[alpha_d$timepoint == "49",]$gg_group2, p.adjust.method = "BH")
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "56",], group = "gg_group2", measure = "InvSimpson")
  pairwise.wilcox.test(alpha_d[alpha_d$timepoint == "56",]$InvSimpson, alpha_d[alpha_d$timepoint == "56",]$gg_group2, p.adjust.method = "BH")
  
  verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "final",], group = "gg_group2", measure = "InvSimpson")
  pairwise.wilcox.test(alpha_d[alpha_d$timepoint == "final",]$InvSimpson, alpha_d[alpha_d$timepoint == "final",]$gg_group2, p.adjust.method = "BH")
  
  
}

# Beta diversity
{
  # For diet timepoints
  sample_data(ps_flt_diet)$diet <- factor(sample_data(ps_flt_diet)$diet, labels = c("50 ppm","500 ppm"))
  # Bray curtis filtered
  betaDiversityTimepoint2Factors(ps_flt_diet, sample_id = "sample_id", timeVariable = "week",
                                 varToCompare =  "diet", distMethod ="bray",
                                 transform = "rel_ab", customColors = c("blue","red"),
                                 font = "Arial", path = "../figures/Thibault_abx/diet/beta_diversity/filtered/", 
                                 additionnalAes = my_theme(), dim = c(4,5), displayPValue = TRUE)
  
  # For diet + treatment at t56 and tfinal
  sample_data(ps_flt_abx)$gg_group2 <- factor(sample_data(ps_flt_abx)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl", "50 ppm Abx", "500 ppm Abx"))
  # Bray curtis filtered
  betaDiversityTimepoint2Factors(ps_flt_abx, sample_id = "sample_id", timeVariable = "timepoint",
                                 varToCompare =  "gg_group2", distMethod ="bray",
                                 customColors = c("blue","red","deepskyblue", "brown1"),
                                 font = "Arial", path = "../figures/Thibault_abx/diet_abx/beta_diversity/filtered/",
                                 additionnalAes = my_theme())
  
  # For abx groups only, t49, t56 and tfinal
  ps_sub <- prune_samples(sample_data(ps_abx_relab_flt)$treatment == "abx", ps_abx_relab_flt)
  # Bray curtis filtered
  betaDiversityTimepoint2Factors(ps_sub, sample_id = "sample_id", timeVariable = "timepoint",
                                 varToCompare =  "diet", distMethod ="bray",
                                 customColors = c("deepskyblue", "brown1"),
                                 font = "Arial", path = "../figures/Thibault_abx/diet_abx/beta_diversity/abxOnly/",
                                 additionnalAes = my_theme(), dim = c(4,5), displayPValue = TRUE)
  
  
  # For abx groups only, comparison between t56 and tfinal
  ps_sub <- prune_samples(sample_data(ps_flt_abx)$treatment == "abx", ps_flt_abx)
  sample_data(ps_sub)$gg_group <- factor(sample_data(ps_sub)$gg_group, levels = c("56:50:abx","56:500:abx","final:50:abx","final:500:abx"))
  # Bray curtis filtered
  betaDiversityTimepoint2Factors(ps_sub, sample_id = "sample_id", timeVariable = "treatment",
                                 varToCompare =  "gg_group", distMethod ="bray",
                                 customColors = c("deepskyblue", "brown1","deepskyblue4", "brown4"),
                                 font = "Arial", path = "../figures/Thibault_abx/diet_abx/beta_diversity/abxOnly/",
                                 additionnalAes = my_theme(), dim = c(4,5), displayPValue = FALSE)
  

}

# Relative abundance analysis: finding differential abundant bugs at the species level, for diet groups only
{
  # Path where to save graphs
  pathToSave <- "~/Documents/CHUM_git/figures/Thibault_abx/relative_abundance_diet/"
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
    relabSingleTimepoint(ps_subset, deseq_subset, measure = "log2fold", "diet", timePoint = timePoint, taxa = "Species", threshold = 0.05, LDA = TRUE, FDR = TRUE, customColors = customColors, path = newPath)
    
    # log2fold change graph based on deseq2 results
    #log2foldChangeGraphSingleTimepoint(ps_subset, deseq_subset, timePoint = timePoint, taxa = "Species", threshold = 0.05, customColors = customColors, customPhylaColors = customPhylaColors, path = newPath, dim =c(6,10))
    
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
  pathToSave <- "~/Documents/CHUM_git/figures/Thibault_abx/relative_abundance_abx_diet_all_groups/"
  existingDirCheck(pathToSave)
  
  #customColors for graph display
  customColors = c("blue", "red", "deepskyblue", "brown1")
  
  #Iterate through timepoints
  for(timePoint in levels(sample_data(ps_abx_relab_flt)$timepoint)){
    
    #New path created for each week
    newPath <- paste(pathToSave, "timepoint_", timePoint, "/", sep = "")
    existingDirCheck(newPath)
    
    #Creating phyloseq objects for each timepoint
    ps_subset <- prune_samples(sample_data(ps_abx_relab_flt)$timepoint == timePoint, ps_abx_relab_flt)
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
                  list("50:water","500:water"), list("50:abx","500:abx"),
                  list("50:water","50:abx"), list("500:water","500:abx")),
                path = newPath, single_factor_design = FALSE,
                dim = c(4,4.5), displayPvalue = FALSE, displaySignificance = TRUE, additionnalAes =
                  list(scale_x_discrete(labels = c("50 ppm\ncontrol","500 ppm\ncontrol","50 ppm\nabx","500 ppm\nabx")),
                       my_theme(),
                       labs(color = "", x=""))) # Include axis lines  # Include axis bar)
  }
  
  #customColors for graph display
  customColors = c("blue", "red", "deepskyblue", "brown1")
  
  #At other taxonomic levels
  taxonomicLevels <- c("Genus","Family","Order","Class","Phylum")
  
  #Iterate through timepoints
  for(timePoint in levels(sample_data(ps_abx_relab_flt)$timepoint)){
    
    #New path created for each week
    newPath <- paste(pathToSave, "timepoint_", timePoint, "/", sep = "")
    existingDirCheck(newPath)
    
    #Creating phyloseq objects for each timepoint
    ps_subset <- prune_samples(sample_data(ps_abx_relab_flt)$timepoint == timePoint, ps_abx_relab_flt)
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
                    list("50:water","500:water"), list("50:abx","500:abx"),
                    list("50:water","50:abx"), list("500:water","500:abx")),
                  path = newPath, single_factor_design = FALSE,
                  dim = c(5,5), displayPvalue = FALSE, displaySignificance = TRUE, additionnalAes =
                    list(scale_x_discrete(labels = c("50 ppm\ncontrol","500 ppm\ncontrol","50 ppm\nabx","500 ppm\nabx")),
                         my_theme(),
                         labs(color = "", x="")))  
    }
  }
}

# Relative abundance analysis: finding differential abundant bugs at different taxonomical levels - only abx groups
{
  pathToSave <- "~/Documents/CHUM_git/figures/Thibault_abx/relative_abundance_abx/"
  existingDirCheck(pathToSave)
  
  #customColors for graph display
  customColors = c("deepskyblue", "brown1")
  customPhylaColors = c("#e6550d","#31a354", "#583093")
  
  #Iterate through timepoints
  for(timePoint in levels(sample_data(ps_abx_relab_flt)$timepoint)[3]){
    
    #New path created for each week
    newPath <- paste(pathToSave, "timepoint_", timePoint, "/", sep = "")
    existingDirCheck(newPath)
    
    #Creating phyloseq objects for each timepoint
    ps_subset <- prune_samples(sample_data(ps_abx_relab_flt)$timepoint == timePoint & sample_data(ps_abx_relab_flt)$treatment == "abx", ps_abx_relab_flt)
    ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)
    print(length(taxa_sums(ps_subset)))
    
    #Simple deseq object only accounting for the differences in diet
    deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ diet) 
    deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric") #Performing the deseq analysis
    print(resultsNames(deseq_subset))
    
    # For a given taxononical levels, creates graph for each timepoint, displaying which species were found to be differentially abundant
    relabSingleTimepoint(ps_subset, deseq_subset, measure = "log2fold", varToCompare = "diet",
                         timePoint = timePoint, taxa = "Species", threshold = 0.05, FDR = TRUE, blockFactor = FALSE,
                         LDA = TRUE, customColors = customColors, path = newPath, displayPvalue = TRUE, additionnalAes = my_theme(), displaySampleID = FALSE)
    

    log2foldChangeGraphSingleTimepoint(ps_subset, deseq_subset, timePoint = timePoint, taxa = "Species", threshold = 0.05, customColors = customColors, customPhylaColors = customPhylaColors, path = newPath, dim =c(6,10))
  }
  
  # At other taxonomic levels
  taxonomicLevels <- c("Genus","Family","Order","Class","Phylum")
  
  # Iterate through timepoints
  for(timePoint in levels(sample_data(ps_abx_relab_flt)$timepoint)[3]){
    
    # New path created for each week
    newPath <- paste(pathToSave, "timepoint_", timePoint, "/", sep = "")
    existingDirCheck(newPath)
    
    # Creating phyloseq objects for each timepoint
    ps_subset <- prune_samples(sample_data(ps_abx_relab_flt)$timepoint == timePoint & sample_data(ps_abx_relab_flt)$treatment == "abx", ps_abx_relab_flt)
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
# For diet only timepoints
{
# Define factor that is combination of diet and timepoint for graph visualization
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
existingDirCheck("../figures/Thibault_abx/stackbar")
ggsave(plot = p, filename = "../figures/Thibault_abx/stackbar/diet_stackbar.png", width = 9, height = 7, dpi = 300)
writeStackbarExtendedSigTable(main_table = diet_phyla_fam$significant_table_main, includeSubTable = TRUE, sub_table = diet_phyla_fam$significant_table_sub, filepath = "../figures/Thibault_abx/stackbar/diet_stackbar_stats.xlsx")

# pvalues heatmap for the main lvl stats
p = pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_abx/stackbar/diet_stackbar_stats.xlsx")),
            selected_comparisons = c("50:0_vs_500:0", "50:35_vs_500:35","50:49_vs_500:49"), displayChangeArrows = FALSE, displayPValues = FALSE,
            txn_lvl="Phylum", lvl = "main", taxons = diet_phyla_fam$main_names[!grepl("Others", x = diet_phyla_fam$main_names)], group = "gg_group", path)
p <- p+scale_x_discrete(labels = c("50 vs 500 3w", "50 vs 500 8w", "50 vs 500 10w", "50 vs 500 14w"))+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 11))
ggsave(plot = p, filename = "../figures/Thibault_abx/stackbar/diet_stackbar_main_table.png", width = 3, height = 5, dpi = 300, bg = "white")


# pvalues heatmap for the sub lvl stats
p = pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_abx/stackbar/diet_stackbar_stats.xlsx")),
                selected_comparisons = c("50:0_vs_500:0", "50:35_vs_500:35","50:49_vs_500:49"),
                txn_lvl="Family", lvl = "sub", taxons =  diet_phyla_fam$sub_names, group = "gg_group", displayPValues = FALSE, displayChangeArrows = TRUE, path) # You can add [!grepl("Others", x = iron_exp_family$sub_names)] to remove "others"
p <- p+scale_x_discrete(labels = c("50 vs 500 3w", "50 vs 500 8w", "50 vs 500 10w", "50 vs 500 14w"))+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 11))
ggsave(plot = p, filename = "../figures/Thibault_abx/stackbar/diet_stackbar_sub_table.png", width = 3.5, height = 8, dpi = 300, bg = "white")
}

# For diet_abx timepoints
{
# First, timepoints and groups must be ordered properly and as factors
sample_data(ps_t56_flt)$diet 
sample_data(ps_t56_flt)$treatment 
sample_data(ps_t56_flt)$gg_group2 <- factor(sample_data(ps_t56_flt)$gg_group2, levels = c("50:water", "50:abx", "500:water","500:abx"))
sample_data(ps_t56_flt)$timepoint

# Selected comparisons should be a number of four and follow design as: "
diet_abx_phyla_fam <- plot_microbiota_2Fac(
  ps_object = ps_t56_flt,
  exp_group = "gg_group2",
  twoFactor = TRUE,
  fac1 = "diet",
  refFac1 = "50",
  fac2 = "treatment",
  refFac2 = "water",
  sample_name = 'sample_id',
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
  selected_comparisons = list(c("50:water", "50:abx"),
                              c("500:water", "500:abx"),
                              c("50:water", "500:water"),
                              c("50:abx", "500:abx")),
  showOnlySubLegend = FALSE
)

print(diet_abx_phyla_fam$plot)
print(diet_abx_phyla_fam$significant_table_main)
print(diet_abx_phyla_fam$significant_table_sub)

# Custom the plot
p <- diet_abx_phyla_fam$plot + 
  facet_wrap2(~ gg_group2, 
              scales  = "free_x", nrow = 2, ncol = 2,
              strip = strip_themed(background_x = elem_list_rect(fill = c("blue","deepskyblue","red","brown1"))),
              labeller = as_labeller(c("50:water" = "50 ppm Ctrl",
                                       "50:abx" = "50 ppm Abx",
                                       "500:water" = "500 ppm Ctrl",
                                       "500:abx" = "500 ppm Abx")))+
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
existingDirCheck("../figures/Thibault_abx/stackbar")
ggsave(plot = p, filename = "../figures/Thibault_abx/stackbar/diet_abx_stackbar.png", width = 7, height = 7, dpi = 300)
writeStackbarExtendedSigTable(main_table = diet_abx_phyla_fam$significant_table_main, includeSubTable = TRUE, sub_table = diet_abx_phyla_fam$significant_table_sub, filepath = "../figures/Thibault_abx/stackbar/diet_abx_stackbar_stats.xlsx")

# pvalues heatmap for the main lvl stats
p = pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_abx/stackbar/diet_abx_stackbar_stats.xlsx")),
            selected_comparisons = c("50:water_vs_50:abx", "500:water_vs_500:abx","50:water_vs_500:water","50:abx_vs_500:abx"), displayChangeArrows = TRUE, displayPValues = FALSE,
            txn_lvl="Phylum", lvl = "main", taxons = diet_abx_phyla_fam$main_names[!grepl("Others", x = diet_abx_phyla_fam$main_names)], group = "gg_group2", path)
p <- p+scale_x_discrete(labels = c("50 Ctrl vs 50 Abx", "500 Ctrl vs 500 Abx", "50 Ctrl vs 500 Ctrl", "50 Abx vs 500 Abx"))+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 11, face = "bold"))
ggsave(plot = p, filename = "../figures/Thibault_abx/stackbar/diet_abx_stackbar_main_table.png", width = 3.5, height = 5, dpi = 300, bg = "white")


# pvalues heatmap for the sub lvl stats
p = pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_abx/stackbar/diet_abx_stackbar_stats.xlsx")),
                selected_comparisons = c("50:water_vs_50:abx", "500:water_vs_500:abx","50:water_vs_500:water","50:abx_vs_500:abx"),
                txn_lvl="Family", lvl = "sub", taxons =  diet_abx_phyla_fam$sub_names, group = "gg_group2", displayPValues = FALSE, displayChangeArrows = TRUE, path) # You can add [!grepl("Others", x = iron_exp_family$sub_names)] to remove "others"
p <- p+scale_x_discrete(labels = c("50 Ctrl vs 50 abx", "500 Ctrl vs 500 Abx", "50 Ctrl vs 500 Ctrl", "50 abx vs 500 Abx"))+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 11, face = "bold"))
ggsave(plot = p, filename = "../figures/Thibault_abx/stackbar/diet_abx_stackbar_sub_table.png", width = 3.5, height = 8, dpi = 300, bg = "white")
}

# For diet_abx last timepoint
{
sample_data(ps_tfinal_flt)$diet 
sample_data(ps_tfinal_flt)$treatment 
sample_data(ps_tfinal_flt)$gg_group2 <- factor(sample_data(ps_tfinal_flt)$gg_group2, levels = c("50:water", "50:abx", "500:water", "500:abx"))
sample_data(ps_tfinal_flt)$timepoint

# Selected comparisons should be a number of four and follow design as: "
final_phyla_fam <- plot_microbiota_2Fac(
  ps_object = ps_tfinal_flt,
  exp_group = "gg_group2",
  twoFactor = TRUE,
  fac1 = "diet",
  refFac1 = "50",
  fac2 = "treatment",
  refFac2 = "water",
  sample_name = 'sample_id',
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
  selected_comparisons = list(c("50:water", "50:abx"),
                              c("500:water", "500:abx"),
                              c("50:water", "500:water"),
                              c("50:abx", "500:abx")),
  showOnlySubLegend = FALSE
)

print(final_phyla_fam$plot)
print(final_phyla_fam$significant_table_main)
print(final_phyla_fam$significant_table_sub)

# Custom the plot
p <- final_phyla_fam$plot + 
  facet_wrap2(~ gg_group2, 
              scales  = "free_x", nrow = 2, ncol = 2,
              strip = strip_themed(background_x = elem_list_rect(fill = c("blue","deepskyblue","red","brown1"))),
              labeller = as_labeller(c("50:water" = "50 ppm Ctrl",
                                       "50:abx" = "50 ppm Abx",
                                       "500:water" = "500 ppm Ctrl",
                                       "500:abx" = "500 ppm Abx")))+
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
existingDirCheck("../figures/Thibault_abx/stackbar")
ggsave(plot = p, filename = "../figures/Thibault_abx/stackbar/final_stackbar.png", width = 7, height = 7, dpi = 300)
writeStackbarExtendedSigTable(main_table = final_phyla_fam$significant_table_main, includeSubTable = TRUE, sub_table = final_phyla_fam$significant_table_sub, filepath = "../figures/Thibault_abx/stackbar/final_stackbar_stats.xlsx")

# pvalues heatmap for the main lvl stats
p = pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_abx/stackbar/final_stackbar_stats.xlsx")),
            selected_comparisons = c("50:water_vs_50:abx", "500:water_vs_500:abx","50:water_vs_500:water","50:abx_vs_500:abx"), displayChangeArrows = TRUE, displayPValues = FALSE,
            txn_lvl="Phylum", lvl = "main", taxons = final_phyla_fam$main_names[!grepl("Others", x = final_phyla_fam$main_names)], group = "gg_group2", path)
p <- p+scale_x_discrete(labels = c("50 Ctrl vs 50 Abx","500 Ctrl vs 500 Abx","50 Ctrl vs 500 Ctrl","50 Abx vs 500 Abx"))+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 11))
ggsave(plot = p, filename = "../figures/Thibault_abx/stackbar/final_stackbar_main_table.png", width = 3.5, height = 5, dpi = 300, bg = "white")


# pvalues heatmap for the sub lvl stats
p = pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_abx/stackbar/final_stackbar_stats.xlsx")),
                selected_comparisons = c("50:water_vs_50:abx", "500:water_vs_500:abx","50:water_vs_500:water","50:abx_vs_500:abx"),
                txn_lvl="Family", lvl = "sub", taxons =  final_phyla_fam$sub_names, group = "gg_group2", displayPValues = FALSE, displayChangeArrows = TRUE, path) # You can add [!grepl("Others", x = iron_exp_family$sub_names)] to remove "others"
p <- p+scale_x_discrete(labels = c("50 Ctrl vs 50 Abx","500 Ctrl vs 500 Abx","50 Ctrl vs 500 Ctrl","50 Abx vs 500 Abx"))+
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 11))
ggsave(plot = p, filename = "../figures/Thibault_abx/stackbar/final_stackbar_sub_table.png", width = 3.5, height = 8, dpi = 300, bg = "white")

}




# Prepare picrust2 input for last timepoint
ps_subset <- prune_samples(sample_data(ps_dss_relab_flt)$timepoint == "final", ps_dss_relab_flt)

# Filtering
# Function filtering out ASVs for which they were in total less than a threshold count
ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
# Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
producePicrust2Inputs(ps_subset, "~/Documents/CHUM_git/Microbiota_18/")




# Prepare picrust2 input for t35
ps_subset <- prune_samples(sample_data(ps_flt_diet)$timepoint == "35", ps_flt_diet)

# Filtering
# Function filtering out ASVs for which they were in total less than a threshold count
ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
# Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
existingDirCheck("~/Documents/CHUM_git/Microbiota_18/t35")
producePicrust2Inputs(ps_subset, "~/Documents/CHUM_git/Microbiota_18/t35/")




# Prepare picrust2 input for t49
ps_subset <- prune_samples(sample_data(ps_dss_relab_flt)$timepoint == "49", ps_dss_relab_flt)

# Filtering
# Function filtering out ASVs for which they were in total less than a threshold count
ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
# Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
existingDirCheck("~/Documents/CHUM_git/Microbiota_18/t49")
producePicrust2Inputs(ps_subset, "~/Documents/CHUM_git/Microbiota_18/t49/")





# Prepare picrust2 input for t49
ps_subset <- prune_samples(sample_data(ps_dss_relab_flt)$timepoint == "54", ps_dss_relab_flt)

# Filtering
# Function filtering out ASVs for which they were in total less than a threshold count
ps_subset <- prune_taxa(taxa_sums(ps_subset) > 10, ps_subset)
# Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
ps_subset <- prune_taxa(colSums(otu_table(ps_subset) > 0) >= (0.05 * nsamples(ps_subset)), ps_subset)
existingDirCheck("~/Documents/CHUM_git/Microbiota_18/t54")
producePicrust2Inputs(ps_subset, "~/Documents/CHUM_git/Microbiota_18/t54/")


























































































# New tests for beta diversity

# Beta diversity (for only diets timepoints)
# Bray curtis filtered (good results)
betaDiversityTimepoint2Factors(ps_flt_diet, sample_id = "sample_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="bray",
                               transform = "rel_ab", customColors = c("blue","red"),
                               font = "Arial", path = "../figures/Thibault_abx/diet/beta_diversity/filtered/")

# Beta diversity (for diet + treatment)
# Bray curtis filtered (good for t54, less interpretable for tfinal)
betaDiversityTimepoint2Factors(ps_flt_dss, sample_id = "sample_id", timeVariable = "timepoint",
                               varToCompare =  "gg_group2", distMethod ="bray",
                               customColors = c("blue","red","blue4","red4"),
                               font = "Arial", path = "../figures/Thibault_abx/diet_dss/beta_diversity/filtered/")

# Beta diversity, for control only at tfinal
ps_subset <- prune_samples(sample_data(ps_flt_dss)$timepoint == "final" & sample_data(ps_flt_dss)$treatment == "water", ps_flt_dss) # Create subset
# Bray curtis filtered 
betaDiversityTimepoint2Factors(ps_subset, sample_id = "sample_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="bray",
                               transform = "rel_ab", customColors = c("blue","red"),
                               font = "Arial", path = "../figures/Thibault_abx/diet/beta_diversity/filtered/")

# Beta diversity, for dss only at t54 and tfinal
ps_subset <- prune_samples(sample_data(ps_flt_dss)$treatment == "dss", ps_flt_dss) # Create subset
# Bray curtis filtered (results show that they are more different at the end of the recovery than at last timepoint of dss)
betaDiversityTimepoint2Factors(ps_subset, sample_id = "sample_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="bray",
                               transform = "rel_ab", customColors = c("darkblue","darkred"),
                               font = "Arial", path = "../figures/Thibault_abx/dss/beta_diversity/filtered/")

# Beta diversity, for dss only at t54 and tfinal, test with Jaccard index 
ps_subset <- prune_samples(sample_data(ps_flt_dss)$treatment == "dss", ps_flt_dss) # Create subset
length(taxa_sums(ps))
ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)
length(taxa_sums(ps))
betaDiversityTimepoint2Factors(ps_subset, sample_id = "sample_id", timeVariable = "timepoint",
                               varToCompare =  "diet", distMethod ="jaccard",
                               transform = "none", customColors = c("darkblue","darkred"),
                               font = "Arial", path = "../figures/Thibault_abx/dss/beta_diversity/jaccard/")

# Beta diversity, for 50 only at t54 and tfinal (showing that they remain different)
ps_subset <- prune_samples(sample_data(ps_dss)$diet == "50", ps_dss) # Create subset
# Bray curtis filtered (results show that they are more different at the end of the recovery than at last timepoint of dss)
betaDiversityTimepoint2Factors(ps_subset, sample_id = "sample_id", timeVariable = "timepoint",
                               varToCompare =  "treatment", distMethod ="bray",
                               transform = "rel_ab", customColors = c("blue","darkblue"),
                               font = "Arial", path = "../figures/Thibault_abx/diet_50/beta_diversity/filtered/")

# Beta diversity, for 500 only at t54 and tfinal (showing that they remain different)
ps_subset <- prune_samples(sample_data(ps_flt_dss)$diet == "500", ps_flt_dss) # Create subset
# Bray curtis filtered (results show that they are more different at the end of the recovery than at last timepoint of dss)
betaDiversityTimepoint2Factors(ps_subset, sample_id = "sample_id", timeVariable = "timepoint",
                               varToCompare =  "treatment", distMethod ="bray",
                               transform = "rel_ab", customColors = c("red","darkred"),
                               font = "Arial", path = "../figures/Thibault_abx/diet_500/beta_diversity/filtered/")


# Now we apply dbRDA for more complex designs with all groups for t54 and tfinal

ps_subset <- prune_samples(sample_data(ps_flt_dss)$timepoint == "54", ps_flt_dss)
betaDiversityTimepointsGroupedDbRDA(ps = ps_subset, sample_id = "sample_id", varToCompare = "gg_group2",distMethod = "hellinger",
                                    customColors = c("blue","red","blue4","red4"), formula = "treatment * diet + cage", transform = "rel_ab",
                                    path = "../figures/Thibault_abx/diet_dss/dbRDA_t54/", font = "Arial", additionnalAes = NULL)

ps_subset <- prune_samples(sample_data(ps_flt_dss)$timepoint == "final", ps_flt_dss)
betaDiversityTimepointsGroupedDbRDA(ps = ps_subset, sample_id = "sample_id", varToCompare = "gg_group2",distMethod = "hellinger",
                                    customColors = c("blue","red","blue4","red4"), formula = "treatment * diet + cage", transform = "rel_ab",
                                    path = "../figures/Thibault_abx/diet_dss/dbRDA_final/", font = "Arial", additionnalAes = NULL)


















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
source("~/Documents/CHUM_git/gut-microbiota-iron/pipelinelinux/microbiota_analysis/chronobiome.R")

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
relabTimepoints(ps_flt_diet, deseq_subset, timeVariable = "week", varToCompare = "diet", taxa = "Species", threshold == 0.05, customColors = c("blue","red"), path = "~/Documents/CHUM_git/figures/Thibault_abx/lrt_diet/")
# Subset for 50 ppm, and for family level taxa
ps_subset <- prune_samples(sample_data(ps_flt_diet)$diet == "50", ps_flt_diet)
ps_taxa <- tax_glom(ps_subset, taxrank = "Phylum")
deseq_subset <- phyloseq_to_deseq2(ps_taxa, ~ week)
deseq_subset <- DESeq(deseq_subset, test="LRT", reduced = ~1)