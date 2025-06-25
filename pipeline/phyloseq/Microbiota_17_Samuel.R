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
  library(caret)
  
}


#Load custom functions for microbiota analysis
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/utilities.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/alpha_diversity_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/beta_diversity_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/correlation_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/relab_analysis_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/taxa_distrib_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/plot_microbiota_extension.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/plot_microbiota_ext.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/deseq2_log2fold_change_analysis.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/other scripts/dataManipFunctions.R")

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

# Alpha diversity
{
  richness_data <- estimate_richness(ps_samuel, measures = c("Shannon", "InvSimpson","Chao1"))
  richness_data <- cbind(as.data.frame(sample_data(ps_samuel)), richness_data)
  richness_data$treatment <- factor(richness_data$treatment, levels = c("Vehicle","Putrescine"))
  verifyStatsAssumptions(richness_data, "treatment", measure = "Chao1")
  t.test(Chao1 ~ treatment, data = richness_data, var.equal = TRUE)
  ironBoxplot(richness_data, "Chao1", group = "treatment", title = "Chao1",y_axis_title = "Chao1", custom_colors = c('black','#A22004'), stats = TRUE, test_results = "**",text_sizes = 5)
  
  verifyStatsAssumptions(richness_data, "treatment", measure = "Shannon")
  t.test(Shannon ~ treatment, data = richness_data, var.equal = TRUE)
  ironBoxplot(richness_data, "Shannon", group = "treatment", title = "Shannon index",y_axis_title = "Shannon index", custom_colors = c('black','#A22004'), stats = TRUE, test_results = "*",text_sizes = 5)
  
  verifyStatsAssumptions(richness_data, "treatment", measure = "InvSimpson")
  t.test(InvSimpson ~ treatment, data = richness_data, var.equal = TRUE)
  ironBoxplot(richness_data, "InvSimpson", group = "treatment", title = "Inverse simpson index",y_axis_title = "Inverse simpson index", custom_colors = c('black','#A22004'), stats = TRUE, test_results = "*",text_sizes = 5)
}

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
  betaDiversityPairwise(ps_samuel, "treatment", pairs = list(list("Vehicle","Putrescine")), "bray", customColors, font = "Times New Roman", displayPValue = TRUE, dim = c(3.75,3.75), transform = "none", path = "~/Documents/CHUM_git/figures/Samuel_final_wt/beta_diversity/")
  betaDiversityPairwise(ps_samuel, "treatment", pairs = list(list("Vehicle","Putrescine")), "wunifrac", customColors, font = "Times New Roman", displayPValue = FALSE, dim = c(3.75,3.75), transform = "none", path = "~/Documents/CHUM_git/figures/Samuel_final_wt/beta_diversity/", title = FALSE)
}

# Differential abundance analysis
{
  deseq_samuel <- phyloseq_to_deseq2(ps_samuel, ~ treatment)
  deseq_samuel <- DESeq(deseq_samuel, test="Wald", fitType = "parametric")
  resultsNames(deseq_samuel)
  customColors = c('black','#A22004')
  relabSingleTimepoint(ps_samuel, deseq_samuel, measure = "log2fold", "treatment", timePoint = "diet", taxa = "Species", displayPvalue = TRUE, LDA = TRUE, threshold = 0.05, customColors = customColors, path = "~/Documents/CHUM_git/figures/Samuel_final_wt/differential_abundance/", dim = c(5,5))
  
  #At other taxonomic levels
  taxonomicLevels <- c("Genus","Family","Order","Class","Phylum")
  
  for(txnLevel in taxonomicLevels){
    
    #Creates ps subset for taxonomical level of interest
    ps_subset <- tax_glom(ps_samuel, taxrank = txnLevel)
    deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ treatment) 
    deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
    relabSingleTimepoint(ps_subset, deseq_subset, measure = "log2fold", "treatment", timePoint = "diet", taxa = txnLevel, displayPvalue = TRUE, LDA = FALSE, threshold = 0.05, customColors = customColors, path = "~/Documents/CHUM_git/figures/Samuel_final_wt/differential_abundance/", dim = c(5,5))
  }
  
}

# Correlation heatmaps
{
  deseq_samuel <- phyloseq_to_deseq2(ps_samuel, ~ treatment)
  deseq_samuel <- DESeq(deseq_samuel, test="Wald", fitType = "parametric")
  
  #Preparing dataframe for correlation
  variables <- read.xlsx("~/Documents/CHUM_git/figures/old/samuel/correlation/Samuel's Combine Data_Correlation.xlsx")
  rownames(variables) <- variables$X3
  variables <- variables[,c(6:11)]
  colnames(variables) <- c("DAI_end_point","EcNC101_colonisation","lcn-2","il-6","tnf-α","colon_length")
  
  customColors = list('black','#A22004')
  p <- correlation2Var(ps_samuel, deseq_samuel, measure = "log2fold", "treatment", taxa = "Species",
                  displayPvalue = FALSE, threshold = 0.05, path = "~/Documents/CHUM_git/figures/Samuel_final_wt/correlation/",
                  df = variables, global = TRUE, showIndivCor = TRUE, transformation = "rel_ab",
                  displayOnlySig = FALSE, returnMainFig = TRUE, displaySpeciesASVNumber = FALSE, colorsHmap = c("#002bff","#ff002b"))
  
  plot <- p+
    coord_fixed() + # Makes thing squared
    scale_x_discrete(labels = pretty_string, position = "top", expand = c(0,0)) +  # Replace space with newline
    scale_y_discrete(limits =  c("Escherichia-Shigella albertii","Lactobacillus johnsonii","Bifidobacterium pseudolongum","Adlercreutzia equolifaciens"),
                     labels =  c("Escherichia-Shigella albertii","Lactobacillus johnsonii","Bifidobacterium pseudolongum","Adlercreutzia equolifaciens"),
                     expand = c(0,0))+
    geom_tile(color = "black", lwd = 0.75, linetype = 1) +
    labs(y = "")+
    geom_text(aes(label = significance), color = "black", size = 6) +
    theme(
      text = element_text(family = "Times New Roman"),      # Global text settings
      axis.title.x = element_blank(),  # Axis titles
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),   # Axis text
      legend.title = element_text(face = "bold", size = 14, margin = margin(b = 20)),  # Legend title  # Legend text
      legend.text = element_text(size = 12),
      axis.text.x = element_text(angle = -45, hjust = 1, size = 14),
      axis.text.y = element_text(face = "bold.italic", size = 14),
      legend.key.height = unit(1.25, "cm"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank())
  plot
  ggsave(filename = "~/Documents/CHUM_git/figures/Samuel_final_wt/correlation/species_hmap.png", plot = plot, bg = "white", height = 8, width = 8, dpi = 300)
}

sample_data(ps_samuel)$treatment <- factor(sample_data(ps_samuel)$treatment, levels = c("Vehicle","Putrescine"))
# Stackbar extended 
il22_exp_family <- plot_microbiota(
  ps_object = ps_samuel,
  exp_group = 'treatment',
  sample_name = 'sample_id',
  hues = c("Blues","Greens", "Purples", "Oranges"),
  differential_analysis = T,
  sig_lab = T,
  fdr_threshold = 0.05,
  main_level = "Phylum",
  sub_level = "Family",
  n_phy = 4, # number of taxa to show
  )
plot <- il22_exp_family$plot
il22_exp_family$plot

# Replace the sample names by some number as an id
facet_levels <- unique(plot$data$treatment)
num_samples <- 13
facet_scales <- lapply(seq_along(facet_levels), function(i) {
  start_label <- 1+(i-1)*num_samples
  end_label <- num_samples*i
  scale_x_discrete(labels = as.character(start_label:end_label))
})


p <- plot + 
  facetted_pos_scales(x = facet_scales) +
  theme(
    text = element_text(family = "Times New Roman"),      # Global text settings
    strip.text.x = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),  # Main title
    axis.title = element_text(size = 15, face = "bold", margin = margin(t = 30)),  # Axis titles
    axis.text = element_text(size = 12, face = "bold"),   # Axis text
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0, margin = margin(t = -17)),
    axis.title.y = element_text(margin = margin(r = -15, unit = "pt")),
    legend.title = element_text(face = "bold", size = 14),  # Legend title  # Legend text
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.01, "lines")      # Size of item labels
  ) +
  labs(x = "Mouse ID")
p
ggsave(filename = "~/Documents/CHUM_git/figures/Samuel_final_wt/stackbar/stackbar.png", plot = p, bg = "white", height = 6, width = 8, dpi = 300)


# For microbiota 17 - specific to Samuel's data, and il22 ko only
{
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

ps_samuel <- subset_samples(ps, student == "Samuel" & genotype == "IL-22ra1-/-")
sum(taxa_sums(ps_samuel))
length(taxa_sums(ps_samuel))

# Alpha diversity
{
  richness_data <- estimate_richness(ps_samuel, measures = c("Shannon", "InvSimpson","Chao1"))
  richness_data <- cbind(as.data.frame(sample_data(ps_samuel)), richness_data)
  richness_data$treatment <- factor(richness_data$treatment, levels = c("Vehicle","Putrescine"))
  verifyStatsAssumptions(richness_data, "treatment", measure = "Chao1")
  t.test(Chao1 ~ treatment, data = richness_data, var.equal = TRUE)
  ironBoxplot(richness_data, "Chao1", group = "treatment", title = "Chao1",y_axis_title = "Chao1", custom_colors = c('black','#A22004'), stats = TRUE, test_results = "n.s.",text_sizes = 5)
  
  verifyStatsAssumptions(richness_data, "treatment", measure = "Shannon")
  t.test(Shannon ~ treatment, data = richness_data, var.equal = TRUE)
  ironBoxplot(richness_data, "Shannon", group = "treatment", title = "Shannon index",y_axis_title = "Shannon index", custom_colors = c('black','#A22004'), stats = TRUE, test_results = "*",text_sizes = 5)
  
  verifyStatsAssumptions(richness_data, "treatment", measure = "InvSimpson")
  t.test(InvSimpson ~ treatment, data = richness_data, var.equal = TRUE)
  ironBoxplot(richness_data, "InvSimpson", group = "treatment", title = "Inverse simpson index",y_axis_title = "Inverse simpson index", custom_colors = c('black','#A22004'), stats = TRUE, test_results = "**",text_sizes = 5)
}

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
  betaDiversityPairwise(ps_samuel, "treatment", pairs = list(list("Vehicle","Putrescine")), "bray", customColors, font = "Times New Roman", displayPValue = TRUE, dim = c(3.75,3.75), transform = "none", path = "~/Documents/CHUM_git/figures/Samuel_final_il22/beta_diversity/")
  betaDiversityPairwise(ps_samuel, "treatment", pairs = list(list("Vehicle","Putrescine")), "wunifrac", customColors, font = "Times New Roman", displayPValue = TRUE, dim = c(3.75,3.75), transform = "none", path = "~/Documents/CHUM_git/figures/Samuel_final_il22/beta_diversity/")
}

# Differential abundance analysis
{
  deseq_samuel <- phyloseq_to_deseq2(ps_samuel, ~ treatment)
  deseq_samuel <- DESeq(deseq_samuel, test="Wald", fitType = "parametric")
  resultsNames(deseq_samuel)
  customColors = c('black','#A22004')
  relabSingleTimepoint(ps_samuel, deseq_samuel, measure = "log2fold", "treatment", timePoint = "diet", taxa = "Species", displayPvalue = TRUE, LDA = TRUE, threshold = 0.05, customColors = customColors, path = "~/Documents/CHUM_git/figures/Samuel_final_il22/differential_abundance/", dim = c(5,5))
  
  #At other taxonomic levels
  taxonomicLevels <- c("Genus","Family","Order","Class","Phylum")
  
  for(txnLevel in taxonomicLevels){
    
    #Creates ps subset for taxonomical level of interest
    ps_subset <- tax_glom(ps_samuel, taxrank = txnLevel)
    deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ treatment) 
    deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
    relabSingleTimepoint(ps_subset, deseq_subset, measure = "log2fold", "treatment", timePoint = "diet", taxa = txnLevel, displayPvalue = TRUE, LDA = FALSE, threshold = 0.05, customColors = customColors, path = "~/Documents/CHUM_git/figures/Samuel_final_il22/differential_abundance/", dim = c(5,5))
  }
  
}

# Correlation heatmaps
{
  deseq_samuel <- phyloseq_to_deseq2(ps_samuel, ~ treatment)
  deseq_samuel <- DESeq(deseq_samuel, test="Wald", fitType = "parametric")
  
  #Preparing dataframe for correlation
  variables <- read.xlsx("~/Documents/CHUM_git/figures/old/samuel/correlation/Samuel's Combine Data_Correlation.xlsx")
  rownames(variables) <- variables$X3
  variables <- variables[,c(6:11)]
  colnames(variables) <- c("DAI_end_point","EcNC101_Colonisation_end_point","lcn-2","il-6","tnf-alpha","colon_length")
  
  customColors = list('black','#A22004')
  p <- correlation2Var(ps_samuel, deseq_samuel, measure = "log2fold", "treatment", taxa = "Species",
                       displayPvalue = FALSE, threshold = 0.05, path = "~/Documents/CHUM_git/figures/Samuel_final_il22/correlation/",
                       df = variables, global = TRUE, showIndivCor = TRUE, transformation = "rel_ab",
                       displayOnlySig = FALSE, returnMainFig = TRUE, displaySpeciesASVNumber = FALSE, colorsHmap = c("#051671","#820909"))
  
  p <- p+
    coord_fixed() + # Makes thing squared
    scale_x_discrete(labels = pretty_string, position = "top") +  # Replace space with newline
    # scale_x_discrete(labels = pretty_string, position = "top")+
    scale_y_discrete(limits =  c("Gordonibacter urolithinfaciens","Escherichia-Shigella albertii"),
                     labels =  c("Gordonibacter urolithinfaciens","Escherichia-Shigella albertii"))+
    geom_tile(color = "black", lwd = 0.75, linetype = 1) +
    geom_text(aes(label = significance), color = "black", size = 6) +
    theme(
      text = element_text(family = "Times New Roman"),      # Global text settings
      axis.title.x = element_blank(),  # Axis titles
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),   # Axis text
      legend.title = element_text(face = "bold", size = 14),  # Legend title  # Legend text
      axis.text.x = element_text(angle = -45, hjust = 1),
      axis.text.y = element_text(face = "bold.italic"))
  ggsave(filename = "~/Documents/CHUM_git/figures/Samuel_final_il22/correlation/species_hmap.png", plot = p, bg = "white", height = 7, width = 7, dpi = 300)
}

sample_data(ps_samuel)$treatment <- factor(sample_data(ps_samuel)$treatment, levels = c("Vehicle","Putrescine"))
# Stackbar extended 
il22_exp_family <- plot_microbiota(
  ps_object = ps_samuel,
  exp_group = 'treatment',
  sample_name = 'sample_id',
  hues = c("Purples", "Blues", "Oranges", "Greens"),
  differential_analysis = T,
  sig_lab = T,
  fdr_threshold = 0.05,
  main_level = "Phylum",
  sub_level = "Family",
  n_phy = 4, # number of taxa to show
)
plot <- il22_exp_family$plot

# Replace the sample names by some number as an id
facet_levels <- unique(plot$data$treatment)
num_samples <- 13
facet_scales <- lapply(seq_along(facet_levels), function(i) {
  start_label <- 1+(i-1)*num_samples
  end_label <- num_samples*i
  scale_x_discrete(labels = as.character(start_label:end_label))
})


p <- plot + 
  facetted_pos_scales(x = facet_scales) +
  theme(
    text = element_text(family = "Times New Roman"),      # Global text settings
    strip.text = element_text(size = 14, face = "bold"),  # Facet titles
    plot.title = element_text(size = 20, face = "bold"),  # Main title
    axis.title = element_text(size = 15, face = "bold"),  # Axis titles
    axis.text = element_text(size = 12, face = "bold"),   # Axis text
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0),
    legend.title = element_text(face = "bold", size = 14),  # Legend title  # Legend text
    axis.ticks.x = element_blank()
  ) +
  labs(x = "Sample ID")
p
}

                                           