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
  library(ggpicrust2)
  library(readxl)
  
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
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/chronobiome.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/picrust2_graphs.R")


# Set working directory
setwd("~/Documents/CHUM_git/Microbiota_18_19_merged/")
asv_table <- as.data.frame(fread("from_server3/merged_asv_table.csv", sep = ";"))
rownames(asv_table) <- asv_table[,1]  # Use the first column as row names
asv_table <- asv_table[,-1]  # Drop the first column
rownames(asv_table) <- gsub("_16S", "", rownames(asv_table))

# Load taxonomical assignments
taxa <- as.matrix(fread("from_server3/merged_taxa_annotation.csv", sep = ";"))
rownames(taxa) <- taxa[,1]  # Use the first column as row names
taxa <- taxa[,-1]  # Drop the first column

# Metadata M18
{
  #loading metadata of interest
  metadata <- read.csv("../Microbiota_18/metadata/metadata.csv", sep = ";")
  
  
  # Remove the non-metadata stuff (liver measures and stuff)
  metadata <- metadata[,-c(5:8)]
  
  # Remove the letter at the end of id
  metadata$id <- substring(metadata$id, 1, 5)
  
  #adding id col as rownames too
  rownames(metadata) <- metadata$id
  
  # Remove dead mouse
  metadata <- metadata[-46,]
  
  # Extract 16S reads sample ids
  samples <- read.xlsx("../Microbiota_18/metadata/Microbiota_18_samples_2025-01-13.xlsx")
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
  metadata 
  
  # Put full_id as rownames
  rownames(metadata) <- metadata$sample_id
  
  metadata$run <- 1
  
  metadata_M18 <- metadata
}

# Metadata M19
{
  #loading metadata of interest
  metadata <- read.xlsx("../Microbiota_19/metadata/dissection.xlsx")
  metadata <- metadata[,-c(5:8)] # Remove the non-metadata stuff (liver measures and stuff)
  metadata$ID <- substring(metadata$ID, 1, 5) # Remove the letter at the end of id
  rownames(metadata) <- metadata$id #adding id col as rownames too
  metadata <- metadata[-17,] # Remove dead mouse
  colnames(metadata)[2] <- "id"
  
  # Extract 16S reads sample ids
  samples <- read.xlsx("../Microbiota_19/metadata/Microbiota_19_NextSeqReadSet_2025-06-22.xlsx")
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
  
  metadata$run <- 2
  
  metadata_M19 <- metadata
}

# Merge metadata
colnames(metadata_M19)[4] <- "cage"
metadata <- rbind(metadata_M18, metadata_M19)

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

sum(taxa_sums(ps)) # total number of reads
length(taxa_sums(ps)) # total number of ASVs

# Put as factors variables that are going to be used
sample_data(ps)$timepoint <- factor(sample_data(ps)$timepoint, levels = c("0","35","49")) # Put timepoint as factor
sample_data(ps)$week <- factor(sample_data(ps)$week, levels = c("3","8","10")) # Put week as factor
sample_data(ps)$diet <- factor(sample_data(ps)$diet, levels = c("50","500"), labels = c("50 ppm","500 ppm")) # Put diet as factor

# Trying to limit batch effects
{
  library(ConQuR)
  library(doParallel)
  
  taxa <- as.data.frame(otu_table(ps))
  taxa <- merge(taxa, metadata, by = "row.names")
  row.names(taxa) <- taxa[,1]
  taxa <- taxa[,-1]
  taxa$run <- factor(taxa$run, levels = c("1","2"))
  taxa$diet <- factor(taxa$diet, levels = c("50","500"))
  taxa$week <- factor(taxa$week, levels = c("3","8","10"))
  batchid <- taxa[,"run"]
  covar <- taxa[,c("diet","week")]
  taxa_corrected <- ConQuR(tax_tab = as.data.frame(otu_table(ps)), batchid = batchid, covariates = covar, batch_ref = 1)
  
  par(mfrow=c(1, 2))
  Plot_PCoA(TAX = otu_table(ps), factor = batchid, main="Before Correction, Bray-Curtis")
  Plot_PCoA(TAX = taxa_corrected, factor = batchid, main="After Correction, Bray-Curtis")
  
  otu_table(ps) <- otu_table(taxa_corrected, taxa_are_rows = FALSE) 
}

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
  write.tree(tree, file = "~/Documents/CHUM_git/Microbiota_18_19_merged/taxonomy/phylogenetic_tree.newick")
  
  #refinement with maximum likelihood
  {
    # Convert the alignment to a phangorn-compatible format
    alignment_matrix <- as.matrix(alignment)
    
    # Estimate the substitution model (e.g., GTR)
    fitGTR <- pml(tree, data = alignment_matrix)
    
    # Optimize the ML tree
    ml_tree <- optim.pml(fitGTR, model = "GTR", rearrangement = "stochastic")
  }
}

# Add phylogenetic tree to phyloseq object
tree <- read.tree("~/Documents/CHUM_git/Microbiota_18_19_merged/taxonomy/phylogenetic_tree.newick")
ps <- merge_phyloseq(ps, tree)
phy_tree(ps) <- midpoint(tree) # Root the tree

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

# Create phyloseq obejcts that we need for the analysis
ps_diet <- merge_phyloseq(ps_t0, ps_t35, ps_t49)
ps_flt_diet <- merge_phyloseq(ps_t0_flt, ps_t35_flt, ps_t49_flt)

#### Alpha diversity ####
existingDirCheck("../figures/Thibault_iron/alpha_diversity")

# Week as numeric
sample_data(ps_diet)$week  <- as.numeric(as.character(sample_data(ps_diet)$week))

#Estinate richness measures for dataset
richness_data <- estimate_richness(ps_diet, measures = c("Chao1", "Shannon", "InvSimpson"))
alpha_d <- cbind(as.data.frame(sample_data(ps_diet)), richness_data)
custom_colors <- c("blue","red")

graphs <- alphaDiversityTimeline(ps_diet, time = "week", group = "diet", custom_colors)

# Chao1
graphs[[1]]+
  scale_x_continuous(n.breaks = 8)+
  labs(y = "Chao1 Index", x = "Time (weeks)", color = "Diet", fill = "")+
  guides(fill = "none")+
  ylim(0,NA)+
  geom_signif( # For first timepoint
    # comparisons = list(c(groups[1],groups[4])),
    xmin = c(3),           # left box in each timepoint
    xmax = c(3),
    annotations = "n.s.",
    y_position = c(300), #
    tip_length = 0,
    color = "black",
    size = 0,
    textsize = 4,
    margin_top = 0.1, # Moves the top according to this value
    vjust = 0,
  )+geom_signif( # For first timepoint
    # comparisons = list(c(groups[1],groups[4])),
    xmin = c(8),           # left box in each timepoint
    xmax = c(8),
    annotations = "n.s.",
    y_position = c(320), #
    tip_length = 0,
    color = "black",
    size = 0,
    textsize = 4,
    margin_top = 0.1, # Moves the top according to this value
    vjust = 0,
    # fontface = "bold"
  )+geom_signif( # For first timepoint
    # comparisons = list(c(groups[1],groups[4])),
    xmin = c(10),           # left box in each timepoint
    xmax = c(10),
    annotations = "n.s.",
    y_position = c(300), #
    tip_length = 0,
    color = "black",
    size = 0,
    textsize = 4,
    margin_top = 0.1, # Moves the top according to this value
    vjust = 0,
  )

ggsave("../figures/Thibault_iron/alpha_diversity/chao1.png",
       bg = "white",height = 4, width =5, dpi = 300)

# Stats
# Within timepoints
verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "0",], group = "diet", measure = "Chao1")
t.test(Chao1 ~ diet, data = alpha_d[alpha_d$timepoint == "0",], var.equal = TRUE)

verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "35",], group = "diet", measure = "Chao1")
t.test(Chao1 ~ diet, data = alpha_d[alpha_d$timepoint == "35",], var.equal = TRUE)

verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "diet", measure = "Chao1")
t.test(Chao1 ~ diet, data = alpha_d[alpha_d$timepoint == "49",], var.equal = TRUE)

# Between t35 and t49
verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint %in% c("35","49") & alpha_d$diet == "50 ppm",], group = "timepoint", measure = "Chao1")
t.test(Chao1 ~ timepoint, data = alpha_d[alpha_d$timepoint %in% c("35","49") & alpha_d$diet == "50 ppm",], var.equal = TRUE)

verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint %in% c("35","49") & alpha_d$diet == "500 ppm",], group = "timepoint", measure = "Chao1")
t.test(Chao1 ~ timepoint, data = alpha_d[alpha_d$timepoint %in% c("35","49") & alpha_d$diet == "500 ppm",], var.equal = TRUE)

# Shannon
graphs[[2]]+
  scale_x_continuous(n.breaks = 8)+
  labs(y = "Shannon Index", x = "Time (weeks)", color = "Diet", fill = "")+
  guides(fill = "none")+
  ylim(0,NA)+
  geom_signif( # For first timepoint
    # comparisons = list(c(groups[1],groups[4])),
    xmin = c(3),           # left box in each timepoint
    xmax = c(3),
    annotations = "n.s.",
    y_position = c(3.7), #
    tip_length = 0,
    color = "black",
    size = 0,
    textsize = 4,
    margin_top = 0.1, # Moves the top according to this value
    vjust = 0,
  )+geom_signif( # For first timepoint
    # comparisons = list(c(groups[1],groups[4])),
    xmin = c(8),           # left box in each timepoint
    xmax = c(8),
    annotations = "n.s.",
    y_position = c(3), #
    tip_length = 0,
    color = "black",
    size = 0,
    textsize = 4,
    margin_top = 0.1, # Moves the top according to this value
    vjust = 0,
  )+geom_signif( # For first timepoint
    # comparisons = list(c(groups[1],groups[4])),
    xmin = c(10),           # left box in each timepoint
    xmax = c(10),
    annotations = "n.s.",
    y_position = c(2.9), #
    tip_length = 0,
    color = "black",
    size = 0,
    textsize = 4,
    margin_top = 0.1, # Moves the top according to this value
    vjust = 0,
  )

ggsave("../figures/Thibault_iron/alpha_diversity/Shannon.png",
       bg = "white",height = 4, width =5, dpi = 300)

# Stats
# Within timepoints
verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "0",], group = "diet", measure = "Shannon")
wilcox.test(Shannon ~ diet, data = alpha_d[alpha_d$timepoint == "0",])

verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "35",], group = "diet", measure = "Shannon")
t.test(Shannon ~ diet, data = alpha_d[alpha_d$timepoint == "35",], var.equal = TRUE)

verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "diet", measure = "Shannon")
t.test(Shannon ~ diet, data = alpha_d[alpha_d$timepoint == "49",], var.equal = TRUE)

# Between t35 and t49
verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint %in% c("35","49") & alpha_d$diet == "50 ppm",], group = "timepoint", measure = "Shannon")
t.test(Shannon ~ timepoint, data = alpha_d[alpha_d$timepoint %in% c("35","49") & alpha_d$diet == "50 ppm",], var.equal = TRUE)

verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint %in% c("35","49") & alpha_d$diet == "500 ppm",], group = "timepoint", measure = "Shannon")
wilcox.test(Shannon ~ timepoint, data = alpha_d[alpha_d$timepoint %in% c("35","49") & alpha_d$diet == "500 ppm",])

# Inverse Simpson
graphs[[3]]+
  scale_x_continuous(n.breaks = 8)+
  labs(y = "Inverse Simpson Index", x = "Time (weeks)", color = "Diet", fill = "")+
  guides(fill = "none")+
  ylim(0,NA)+
  geom_signif( # For first timepoint
    # comparisons = list(c(groups[1],groups[4])),
    xmin = c(3),           # left box in each timepoint
    xmax = c(3),
    annotations = "n.s.",
    y_position = c(15), #
    tip_length = 0,
    color = "black",
    size = 0,
    textsize = 4,
    margin_top = 0.1, # Moves the top according to this value
    vjust = 0,
  )+geom_signif( # For first timepoint
    # comparisons = list(c(groups[1],groups[4])),
    xmin = c(8),           # left box in each timepoint
    xmax = c(8),
    annotations = "n.s.",
    y_position = c(6), #
    tip_length = 0,
    color = "black",
    size = 0,
    textsize = 4,
    margin_top = 0.1, # Moves the top according to this value
    vjust = 0,
  )+geom_signif( # For first timepoint
    # comparisons = list(c(groups[1],groups[4])),
    xmin = c(10),           # left box in each timepoint
    xmax = c(10),
    annotations = "n.s.",
    y_position = c(5), #
    tip_length = 0,
    color = "black",
    size = 0,
    textsize = 4,
    margin_top = 0.1, # Moves the top according to this value
    vjust = 0,
  )

ggsave("../figures/Thibault_iron/alpha_diversity/InvSimpson.png",
       bg = "white",height = 4, width =5, dpi = 300)

# Stats
# Within timepoints
verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "0",], group = "diet", measure = "InvSimpson")
wilcox.test(InvSimpson ~ diet, data = alpha_d[alpha_d$timepoint == "0",])

verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "35",], group = "diet", measure = "InvSimpson")
wilcox.test(InvSimpson ~ diet, data = alpha_d[alpha_d$timepoint == "35",])

verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "diet", measure = "InvSimpson")
wilcox.test(InvSimpson ~ diet, data = alpha_d[alpha_d$timepoint == "49",])

# Between t35 and t49
verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint %in% c("35","49") & alpha_d$diet == "50 ppm",], group = "timepoint", measure = "InvSimpson")
wilcox.test(InvSimpson ~ timepoint, data = alpha_d[alpha_d$timepoint %in% c("35","49") & alpha_d$diet == "50 ppm",])

verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint %in% c("35","49") & alpha_d$diet == "500 ppm",], group = "timepoint", measure = "InvSimpson")
wilcox.test(InvSimpson ~ timepoint, data = alpha_d[alpha_d$timepoint %in% c("35","49") & alpha_d$diet == "500 ppm",])



#### Beta diversity ####
# Bray curtis filtered
betaDiversityTimepoint2Factors(ps_flt_diet, sample_id = "sample_id", timeVariable = "week",
                               varToCompare =  "diet", distMethod ="bray",
                               transform = "rel_ab", customColors = c("blue","red"),
                               font = "Arial", path = "../figures/Thibault_iron/beta_diversity/filtered/", 
                               additionnalAes = my_theme()+theme(plot.title = element_text(size = 12),
                                                                 axis.line = element_line(color = "black", size = 0.6)),
                               dim = c(3,4), displayPValue = TRUE)

# Weighted unifrac unfiltered
betaDiversityTimepoint2Factors(ps_diet, sample_id = "sample_id", timeVariable = "week",
                               varToCompare =  "diet", distMethod ="wunifrac",
                               transform = "rel_ab", customColors = c("blue","red"),
                               font = "Arial", path = "../figures/Thibault_iron/beta_diversity/filtered/", 
                               additionnalAes = my_theme()+theme(plot.title = element_text(size = 12),
                                                                 axis.line = element_line(color = "black", size = 0.6)),
                               dim = c(3,4), displayPValue = TRUE)

sample_data(ps_flt_diet)$run <- factor(sample_data(ps_flt_diet)$run, labels = c("1","2"))
sample_data(ps_diet)$run <- factor(sample_data(ps_diet)$run, labels = c("1","2"))
sample_data(ps_diet)$week <- factor(sample_data(ps_diet)$week, levels = c("3","8","10"))
betaDiversityTimepoint2Factors(ps_flt_diet, sample_id = "sample_id", timeVariable = "week",
                               varToCompare =  "run", distMethod ="bray",
                               transform = "rel_ab", customColors = c("green","yellow"),
                               font = "Arial", path = "../figures/Thibault_iron/beta_diversity/runs/", 
                               additionnalAes = my_theme()+theme(plot.title = element_text(size = 12)), dim = c(3,4), displayPValue = TRUE)

#### Relative abundance ####
# Path where to save graphs
pathToSave <- "~/Documents/CHUM_git/figures/Thibault_iron/relative_abundance/"
existingDirCheck(pathToSave)

#customColors for graph display
customColors = c("blue","red")
customPhylaColors = c("#583093","#31a354", "#3165a3")

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
  relabSingleTimepoint(ps_subset, deseq_subset, measure = "log2fold", "diet", timePoint = timePoint, taxa = "Species", threshold = 0.05, LDA = FALSE, FDR = TRUE, customColors = customColors, path = newPath, dim = c(5,5))
  
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
    relabSingleTimepoint(ps_taxa, deseq_subset, measure = "log2fold", "diet", timePoint = timePoint, taxa = txnLevel, threshold = 0.05, LDA = FALSE, FDR = TRUE, customColors = customColors, path = newPath, dim = c(5,5)) 
  }
}

#### Relative abundance timeline ####
# Path where to save graphs
pathToSave <- "~/Documents/CHUM_git/figures/Thibault_iron/relative_abundance_timeline/"
existingDirCheck(pathToSave)

#customColors for graph display
customColors = c("blue","red")
  
relabTimelineRevised(ps_flt_diet, timeVariable = "week", varToCompare = "diet", taxa = "Species", threshold = 0.05, customColors = customColors, path = pathToSave)

#At other taxonomic levels
taxonomicLevels <- c("Genus","Family","Order","Class","Phylum")

for(txnLevel in taxonomicLevels){
  
  #Creates ps subset for taxonomical level of interest
  ps_taxa <- tax_glom(ps_flt_diet, taxrank = txnLevel)
  
  #For a given taxononical levels, creates graph for each timepoint, displaying which species were found to be differentially abundant
  relabTimelineRevised(ps_taxa, timeVariable = "week", varToCompare = "diet", taxa = txnLevel, threshold = 0.05, customColors = customColors, path = pathToSave)
}

#### F/B ratio ####
ps_phyla <- tax_glom(ps_diet, taxrank = "Phylum") # Agglom ps at phylum level
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
fb_ratio_df$week <- factor(fb_ratio_df$week, levels = c("3","8","10"))
fb_ratio_df$diet <- factor(fb_ratio_df$diet, levels = c("50","500"), labels = c("50 ppm","500 ppm"))
fbRatioGraphTimeSeries(fb_ratio_df, group = "diet",time = "week", measure = "fb_ratio", custom_colors = c("blue","red"), custom_theme = my_theme())+
  ylim(0,NA)+
  labs(x = "Time (weeks)", fill = "", color = "Diet")+
  guides(fill = "none")
  
fb_ratio_df$week <- as.numeric(as.character(fb_ratio_df$week))
fbRatioGraphTimeSeries2(fb_ratio_df, group = "diet",time = "week", measure = "fb_ratio", custom_colors = c("blue","red"), custom_theme = my_theme())+
  ylim(0,NA)+
  labs(x = "Time (weeks)", fill = "", color = "Diet")+
  scale_x_continuous(n.breaks = 8)+
  guides(fill = "none")+
  geom_signif( # For first timepoint
    # comparisons = list(c(groups[1],groups[4])),
    xmin = c(3),           # left box in each timepoint
    xmax = c(3),
    annotations = "n.s.",
    y_position = c(1.3), #
    tip_length = 0,
    color = "black",
    size = 0,
    textsize = 4,
    margin_top = 0.1, # Moves the top according to this value
    vjust = -0.1,
  )+geom_signif( # For second timepoint
    # comparisons = list(c(groups[1],groups[4])),
    xmin = c(8),           # left box in each timepoint
    xmax = c(8),
    annotations = "*",
    y_position = c(3.8), #
    tip_length = 0,
    color = "black",
    size = 0,
    textsize = 5,
    margin_top = 0.1, # Moves the top according to this value
    vjust = 0,
    fontface = "bold"
  )+geom_signif( # For third timepoint
    # comparisons = list(c(groups[1],groups[4])),
    xmin = c(10),           # left box in each timepoint
    xmax = c(10),
    annotations = "n.s.",
    y_position = c(4.8), #
    tip_length = 0,
    color = "black",
    size = 0,
    textsize = 4,
    margin_top = 0.1, # Moves the top according to this value
    vjust = -0.1,
  )+geom_signif( # For third timepoint
    # comparisons = list(c(groups[1],groups[4])),
    xmin = c(8),           # left box in each timepoint
    xmax = c(10),
    annotations = "*",
    y_position = c(5.2), #
    tip_length = 0.02,
    color = "black",
    size = 0.5,
    textsize = 5,
    margin_top = 0.5, # Moves the top according to this value
    vjust = -0.1,
    fontface = "bold"
  )
ggsave("../figures/Thibault_iron/fb_ratio/fb_ratio_diet.png",
       bg = "white",height = 5, width =5, dpi = 300)
# Stats
verifyStatsAssumptions(fb_ratio_df[fb_ratio_df$timepoint == "0",], group = "diet",measure =  "fb_ratio")
wilcox.test(fb_ratio ~ diet, data = fb_ratio_df[fb_ratio_df$timepoint == "0",])

verifyStatsAssumptions(fb_ratio_df[fb_ratio_df$timepoint == "35",], group = "diet",measure =  "fb_ratio")
wilcox.test(fb_ratio ~ diet, data = fb_ratio_df[fb_ratio_df$timepoint == "35",])

verifyStatsAssumptions(fb_ratio_df[fb_ratio_df$timepoint == "49",], group = "diet",measure =  "fb_ratio")
wilcox.test(fb_ratio ~ diet, data = fb_ratio_df[fb_ratio_df$timepoint == "49",])

# Between t35 and t49
verifyStatsAssumptions(df = fb_ratio_df[fb_ratio_df$timepoint %in% c("35","49") & fb_ratio_df$diet == "50 ppm",], group = "timepoint", measure = "fb_ratio")
wilcox.test(fb_ratio ~ timepoint, data = fb_ratio_df[fb_ratio_df$timepoint %in% c("35","49") & fb_ratio_df$diet == "50 ppm",])

verifyStatsAssumptions(df = fb_ratio_df[fb_ratio_df$timepoint %in% c("35","49") & fb_ratio_df$diet == "500 ppm",], group = "timepoint", measure = "fb_ratio")
wilcox.test(fb_ratio ~ timepoint, data = fb_ratio_df[fb_ratio_df$timepoint %in% c("35","49") & fb_ratio_df$diet == "500 ppm",])



#### Stackbar ####
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
  hues = c("Blues", "Greens", "Purples", "Oranges"), # c("Purples", "Blues", "Reds", "Greens", "Oranges", "Greys", "BuPu")
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
  scale_x_discrete(labels = as.character(1:48)),
  scale_x_discrete(labels = as.character(1:48)),
  scale_x_discrete(labels = as.character(list2)),
  scale_x_discrete(labels = as.character(49:96)),
  scale_x_discrete(labels = as.character(49:96))
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
        axis.text.x = element_text(angle = 45, hjust = 0.5, margin = margin(t = -5), size = 5),
        axis.title.y = element_text(margin = margin(r = -15, unit = "pt")),
        legend.title = element_text(face = "bold", size = 14),  # Legend title  # Legend text
        axis.ticks.x = element_blank()) +
  # scale_x_discrete(labels = function(x) substr(x, 1, 5))+
  facetted_pos_scales(x = facet_scales)+
  labs(x = "Mouse number")
p

# Saving the plot and the associated stats
existingDirCheck("../figures/Thibault_iron/stackbar")
ggsave(plot = p, filename = "../figures/Thibault_iron/stackbar/diet_stackbar.png", width = 12, height = 7, dpi = 300)
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

#### Correlations ####
#Preparing dataframe for correlation
# Load all data
df <- as.data.frame(read_xlsx("~/Documents/CHUM_git/gut-microbiota-iron/experiments/data_both_exp.xlsx"))
df <- df[!is.na(df$id),] # Remove dead mice
df$id[30] <- "10958" # Correct name
df$id <- paste0(substring(df$id, 1, 5),"_T35")
colnames(df)[1] <- "sample_id"
rownames(df) <- df$sample_id
df <- df[,-c(1:3)]
df$nrm_spleen_weight <- df$spleen_weight/df$body_weight # Normalize spleen and liver weights
df$nrm_liver_weight <- df$liver_weight/df$body_weight 

ps_taxa <- tax_glom(ps_t35_flt, taxrank = "Phylum")
deseq_subset <- phyloseq_to_deseq2(ps_taxa, ~ diet)
deseq_subset <- phyloseq_to_deseq2(ps_taxa, ~ diet)
deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
# res <- results(deseq_subset, name = resultsNames(deseq_subset)[2]) #50 vs 500
# print(resultsNames(deseq_subset))

#For species level
#One heatmap for all groups
p <- correlation2Var(ps_taxa, deseq_subset, measure = "log2fold", "diet", taxa = "Phylum", displayPvalue = FALSE, threshold = 0.05, FDR = TRUE, "~/Documents/CHUM_git/figures/Thibault_dss/newTaxAnnotation/correlation_heatmaps/", df = df, global = TRUE, singleVariable = FALSE, showIndivCor = FALSE, transformation = "CLR", displayOnlySig = FALSE, returnMainFig = TRUE, displaySpeciesASVNumber = TRUE)
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
ggsave(filename = "~/Documents/CHUM_git/figures/Thibault_dss/correlation_heatmaps/species_iront35_hmap.png",
       plot = p, bg = "white", height = 6, width = 6, dpi = 300)




# Creating a cauliflower phylogenetic tree
{
  library(ggtree)
  library(ape)
  library(metacoder)
  library(RColorBrewer)
  
  tax_data <- parse_phyloseq(ps)
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
  palette_vec <- colorRampPalette(brewer.pal(min(8, length(unique_phyla)), "Set3"))(length(unique_phyla))
  phylum_pal  <- setNames(palette_vec, unique_phyla)
  tax_data$data$taxa$color <- phylum_pal[prop_phylum]
  
  p_tree <- heat_tree(
    tax_data,
    node_color       = color,
    edge_color       = color,
    node_label       = "",
    node_color_range = phylum_pal,
    edge_color_range = phylum_pal,
    node_legend_title = "Phylum",
    layout           = "davidson-harel"
  )
  
  library(patchwork)  # For combining plots
  
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
  
  ggsave(filename = "~/Documents/CHUM_git/figures/Thibault_iron/cauliflower phytree/phy_tree.png", plot = final_plot, bg = "white", height = 7, width = 11, dpi = 300)
}