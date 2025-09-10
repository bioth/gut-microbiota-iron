# Required packages
{
  require(tidyverse)
  require(ggplot2)
  require(dplyr)
  require(lme4)
  require(car)
  require(ggsignif)
  require(readxl)
  require(openxlsx)
  require(permuco)
  require(patchwork)
  require(legendry)
  require(ggbreak)
  require(lmerTest)
  require(emmeans)
  require(phyloseq)
  require(DESeq2)
  require(Biostrings)
  require(data.table)
  require(ggpubr)
  require(cowplot)
  require(vegan)
  require(DECIPHER)
  require(phangorn)
  require(ape)
  require(writexl)
  require(openxlsx)
  require(reshape2)
  require(Hmisc)
  require(plotly)
  require(StackbarExtended)
  require(ggpattern)
  require(pairwiseAdonis)
  require(caret)
  require(ggh4x)
  require(EnhancedVolcano)
  require(forcats)
}

# Load custom functions for microbiota analysis
source("path/to/repo/utils/utilities.R")
source("path/to/repo/utils/alpha_diversity_graphs_and_stats.R")
source("path/to/repo/utils/beta_diversity_graphs_and_stats.R")
source("path/to/repo/utils/plot_microbiota_extension.R")
source("path/to/repo/utils/chronobiome.R")
source("path/to/repo/utils/relab_analysis_graphs_and_stats.R")

# Load 16S data and build phyloseq objects
{
# Set working directory
asv_table <- as.data.frame(fread("path/to/asv/table", sep = ";"))
rownames(asv_table) <- asv_table[,1]  # Use the first column as row names
asv_table <- asv_table[,-1]  # Drop the first column

# Metadata handling
{
  #loading metadata of interest
  metadata <- read.csv("path/to/metadata", sep = ";")
  
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
  
  # Replace final by t112
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
taxa <- as.matrix(fread("path/to/tax/table", sep = ";"))
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

sum(taxa_sums(ps)) # total number of reads
length(taxa_sums(ps)) # total number of ASVs

# Put as factors variables that are going to be used
sample_data(ps)$gg_group2 <- factor(sample_data(ps)$gg_group2, levels = c("50:water", "500:water", "50:dss", "500:dss")) # Put gg_group2 as factor
sample_data(ps)$timepoint <- factor(sample_data(ps)$timepoint, levels = c("0","35","49","54","112")) # Put timepoint as factor
sample_data(ps)$week <- factor(sample_data(ps)$week, levels = c("0","5","7","7.7","16")) # Put week as factor
sample_data(ps)$treatment <- factor(sample_data(ps)$treatment, levels = c("water","dss")) # Put treatment as factor
sample_data(ps)$diet <- factor(sample_data(ps)$diet, levels = c("50","500")) # Put diet as factor

# Load phylogenetic tree if possible
tree <- read.tree("path/to/phylo/tree")

# Add tree to phyloseq object
ps_tree <- merge_phyloseq(ps, phy_tree(tree))
phy_tree(ps_tree) <- midpoint(tree) # Root the tree

# Create single timepoint phyloseq objects and apply filter
ps_t0 <- prune_samples(sample_data(ps)$timepoint %in% c("0"), ps)
ps_t0_flt <- prune_taxa(taxa_sums(ps_t0) > 10, ps_t0)
ps_t0_flt <- prune_taxa(colSums(otu_table(ps_t0_flt) > 0) >= (0.4 * nsamples(ps_t0_flt)), ps_t0_flt)
length(taxa_sums(ps_t0_flt))

ps_t35 <- prune_samples(sample_data(ps)$timepoint %in% c("35"), ps)
ps_t35_flt <- prune_taxa(taxa_sums(ps_t35) > 10, ps_t35)
ps_t35_flt <- prune_taxa(colSums(otu_table(ps_t35_flt) > 0) >= (0.4 * nsamples(ps_t35_flt)), ps_t35_flt)
length(taxa_sums(ps_t35_flt))

ps_t49 <- prune_samples(sample_data(ps)$timepoint %in% c("49"), ps)
ps_t49_flt <- prune_taxa(taxa_sums(ps_t49) > 10, ps_t49)
ps_t49_flt <- prune_taxa(colSums(otu_table(ps_t49_flt) > 0) >= (0.4 * nsamples(ps_t49_flt)), ps_t49_flt)
length(taxa_sums(ps_t49_flt))

ps_t54 <- prune_samples(sample_data(ps)$timepoint %in% c("54"), ps)
ps_t54_flt <- prune_taxa(taxa_sums(ps_t54) > 10, ps_t54)
ps_t54_flt <- prune_taxa(colSums(otu_table(ps_t54_flt) > 0) >= (0.4 * nsamples(ps_t54_flt)), ps_t54_flt)
length(taxa_sums(ps_t54_flt))

ps_tfinal <- prune_samples(sample_data(ps)$timepoint %in% c("112"), ps)
ps_tfinal_flt <- prune_taxa(taxa_sums(ps_tfinal) > 10, ps_tfinal)
ps_tfinal_flt <- prune_taxa(colSums(otu_table(ps_tfinal_flt) > 0) >= (0.4 * nsamples(ps_tfinal_flt)), ps_tfinal_flt)
length(taxa_sums(ps_tfinal_flt))

# Create phyloseq obejcts that we need for the analysis
ps_diet <- merge_phyloseq(ps_t0, ps_t35, ps_t49)
ps_dss_alpha <- merge_phyloseq(ps_t49, ps_t54, ps_tfinal)
ps_dss_relab_flt <- merge_phyloseq(ps_t49_flt, ps_t54_flt, ps_tfinal_flt)
ps_dss <- merge_phyloseq(ps_t54, ps_tfinal)
ps_flt_diet <- merge_phyloseq(ps_t0_flt, ps_t35_flt, ps_t49_flt)
ps_flt_dss <- merge_phyloseq(ps_t54_flt, ps_tfinal_flt)
ps_flt_all <- merge_phyloseq(ps_t0_flt, ps_t35_flt, ps_t49_flt, ps_t54_flt, ps_tfinal_flt)
}

# Alpha diversity
{
  
  # Alpha diveristy timeline (ctrl vs dss for each diet) + boxplots at last timepoint
  {
    # Week as numeric
    sample_data(ps)$timepoint  <- as.numeric(as.character(sample_data(ps)$timepoint))
    sample_data(ps)$gg_group2 
    sample_data(ps)$gg_group2 <- factor(sample_data(ps)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
    
    ps50 <- prune_samples(sample_data(ps)$diet == "50", ps)
    custom_colors <- c("#95BECF","#325BAD")
    
    # For 50 ppm
    {
      #Estinate richness measures for dataset
      richness_data <- estimate_richness(ps50, measures = c("Chao1", "Shannon", "InvSimpson"))
      alpha_d <- cbind(as.data.frame(sample_data(ps50)), richness_data)
      
      graphs <- alphaDiversityTimeline(ps50, time = "timepoint", group = "gg_group2", shape = "treatment", custom_colors, semRibbons = TRUE)
      
      # Chao1
      p1 <- graphs[[1]]+
        scale_x_continuous(
          breaks = seq(min(sample_data(ps)$timepoint), max(sample_data(ps)$timepoint), by = 7) # breaks every 7 days                    # "T" before each label
        ) +
        labs(y = "Index", x = "Days", color = NULL, fill = "", title = 'Chao1')+
        guides(shape = "none", fill = "none")+
        theme(axis.ticks.x = element_line(),
              legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5),
              axis.text.x = element_text(size = 6)
        )+
        coord_cartesian(ylim = c(0, 1300))+
        geom_text(
          data = data.frame(
            x = c(54, 112),
            y = c(900, 1050),
            label = c("***", "*")
          ),
          aes(x = x, y = y, label = label),
          size = 4,
          inherit.aes = FALSE)
      
      p1 <- p1+
        annotate("rect", xmin = 49, xmax = 54,
                 fill = "gray56", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)+
        annotate("rect", xmin = -Inf, xmax = 49,
                 fill = "#95BECF", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)+
        annotate("rect", xmin = 49, xmax = Inf,
                 fill = "#95BECF", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)
      
      p1 <- putGgLastLayerBack(p1, nLayers = 3)
      
      p1 <- p1+
        scale_y_continuous(labels = scales::label_comma(big.mark = ","))+
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
      p1
      
      # Stats
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "0",], group = "gg_group2", measure = "Chao1")
      wilcox.test(Chao1 ~ treatment, data = alpha_d[alpha_d$timepoint == "0",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "35",], group = "gg_group2", measure = "Chao1")
      wilcox.test(Chao1 ~ treatment, data = alpha_d[alpha_d$timepoint == "35",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "gg_group2", measure = "Chao1")
      wilcox.test(Chao1 ~ treatment, data = alpha_d[alpha_d$timepoint == "49",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "54",], group = "gg_group2", measure = "Chao1")
      wilcox.test(Chao1 ~ treatment, data = alpha_d[alpha_d$timepoint == "54",]) # ***
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "112",], group = "gg_group2", measure = "Chao1")
      wilcox.test(Chao1 ~ treatment, data = alpha_d[alpha_d$timepoint == "112",]) # *
      
      # Shannon
      p2 <- graphs[[2]]+
        scale_x_continuous(
          breaks = seq(min(sample_data(ps)$timepoint), max(sample_data(ps)$timepoint), by = 7) # breaks every 7 days                    # "T" before each label
        ) +
        labs(y = "Index", x = "Days", color = NULL, fill = "", title = 'Shannon')+
        guides(shape = "none", fill = "none")+
        theme(axis.ticks.x = element_line(),
              legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5),
              axis.text.x = element_text(size = 6)
        )+
        coord_cartesian(ylim = c(0, 5))+
        geom_text(
          data = data.frame(
            x = c(54),
            y = c(3.5),
            label = c("**")
          ),
          aes(x = x, y = y, label = label),
          size = 4,
          inherit.aes = FALSE)
      
      p2 <- p2+
        annotate("rect", xmin = 49, xmax = 54,
                 fill = "gray56", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)+
        annotate("rect", xmin = -Inf, xmax = 49,
                 fill = "#95BECF", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)+
        annotate("rect", xmin = 49, xmax = Inf,
                 fill = "#95BECF", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)
      
      p2 <- putGgLastLayerBack(p2, nLayers = 3)
      
      p2 <- p2+
        scale_y_continuous(labels = scales::label_comma(big.mark = ","))+
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
      p2
      
      # Stats
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "0",], group = "gg_group2", measure = "Shannon")
      wilcox.test(Shannon ~ treatment, data = alpha_d[alpha_d$timepoint == "0",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "35",], group = "gg_group2", measure = "Shannon")
      wilcox.test(Shannon ~ treatment, data = alpha_d[alpha_d$timepoint == "35",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "gg_group2", measure = "Shannon")
      wilcox.test(Shannon ~ treatment, data = alpha_d[alpha_d$timepoint == "49",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "54",], group = "gg_group2", measure = "Shannon")
      wilcox.test(Shannon ~ treatment, data = alpha_d[alpha_d$timepoint == "54",]) # **
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "112",], group = "gg_group2", measure = "Shannon")
      wilcox.test(Shannon ~ treatment, data = alpha_d[alpha_d$timepoint == "112",]) # n.s.
      
      # Inverse Simpson
      p3 <- graphs[[3]]+
        scale_x_continuous(
          breaks = seq(min(sample_data(ps)$timepoint), max(sample_data(ps)$timepoint), by = 7) # breaks every 7 days                  # "T" before each label
        ) +
        labs(y = "Index", x = "Days", color = NULL, fill = "", title = 'Inverse Simpson')+
        guides(shape = "none", fill = "none")+
        theme(axis.ticks.x = element_line(),
              legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5),
              axis.text.x = element_text(size = 6))+
        coord_cartesian(ylim = c(0, 28))+
        geom_text(
          data = data.frame(
            x = c(54),
            y = c(11),
            label = c("***")
          ),
          aes(x = x, y = y, label = label),
          size = 4,
          inherit.aes = FALSE)
      
      p3 <- p3+
        annotate("rect", xmin = 49, xmax = 54,
                 fill = "gray56", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)+
        annotate("rect", xmin = -Inf, xmax = 49,
                 fill = "#95BECF", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)+
        annotate("rect", xmin = 49, xmax = Inf,
                 fill = "#95BECF", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)
      
      p3 <- putGgLastLayerBack(p3, nLayers = 3)
      
      p3 <- p3+
        scale_y_continuous(labels = scales::label_comma(big.mark = ","))+
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
      p3
      
      # Stats
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "0",], group = "gg_group2", measure = "InvSimpson")
      wilcox.test(InvSimpson ~ treatment, data = alpha_d[alpha_d$timepoint == "0",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "35",], group = "gg_group2", measure = "InvSimpson")
      wilcox.test(InvSimpson ~ treatment, data = alpha_d[alpha_d$timepoint == "35",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "gg_group2", measure = "InvSimpson")
      wilcox.test(InvSimpson ~ treatment, data = alpha_d[alpha_d$timepoint == "49",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "54",], group = "gg_group2", measure = "InvSimpson")
      wilcox.test(InvSimpson ~ treatment, data = alpha_d[alpha_d$timepoint == "54",]) # ***
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "112",], group = "gg_group2", measure = "InvSimpson")
      wilcox.test(InvSimpson ~ treatment, data = alpha_d[alpha_d$timepoint == "112",]) # n.s.
      
      alpha_d50 <- p1 / p2 /p3 +plot_layout(guides = "collect") 
    }
    
    ps500 <- prune_samples(sample_data(ps)$diet == "500", ps)
    custom_colors <- c("#F2AA84","#B22222")
    
    # For 500 ppm
    {
      #Estinate richness measures for dataset
      richness_data <- estimate_richness(ps500, measures = c("Chao1", "Shannon", "InvSimpson"))
      alpha_d <- cbind(as.data.frame(sample_data(ps500)), richness_data)
      
      graphs <- alphaDiversityTimeline(ps500, time = "timepoint", group = "gg_group2", shape = "treatment", custom_colors, semRibbons = TRUE)
      
      # Chao1
      p1 <- graphs[[1]]+
        scale_x_continuous(
          breaks = seq(min(sample_data(ps)$timepoint), max(sample_data(ps)$timepoint), by = 7) # breaks every 7 days                    # "T" before each label
        ) +
        labs(y = "Index", x = "Days", color = NULL, fill = "", title = 'Chao1')+
        guides(shape = "none", fill = "none")+
        theme(axis.ticks.x = element_line(),
              legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5),
              axis.text.x = element_text(size = 6)
        )+
        coord_cartesian(ylim = c(0, 1300))+
        geom_text(
          data = data.frame(
            x = c(54, 112),
            y = c(1020, 1050),
            label = c("***", "*")
          ),
          aes(x = x, y = y, label = label),
          size = 4,
          inherit.aes = FALSE)
      
      p1 <- p1+
        annotate("rect", xmin = 35, xmax = 49, # 50 ppm
                 fill = "#95BECF", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)+
        annotate("rect", xmin = 54, xmax = Inf, # 50 ppm
                 fill = "#95BECF", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)+
        annotate("rect", xmin = -Inf, xmax = 35, # Exposure
                 fill = "#F2AA84", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)+
        annotate("rect", xmin = 49, xmax = 54, # DSS
                 fill = "gray56", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)
      
      p1 <- putGgLastLayerBack(p1, nLayers = 4)
      
      p1 <- p1+
        scale_y_continuous(labels = scales::label_comma(big.mark = ","))+
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
      p1
      
      # Stats
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "0",], group = "gg_group2", measure = "Chao1")
      wilcox.test(Chao1 ~ treatment, data = alpha_d[alpha_d$timepoint == "0",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "35",], group = "gg_group2", measure = "Chao1")
      wilcox.test(Chao1 ~ treatment, data = alpha_d[alpha_d$timepoint == "35",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "gg_group2", measure = "Chao1")
      wilcox.test(Chao1 ~ treatment, data = alpha_d[alpha_d$timepoint == "49",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "54",], group = "gg_group2", measure = "Chao1")
      wilcox.test(Chao1 ~ treatment, data = alpha_d[alpha_d$timepoint == "54",]) # ***
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "112",], group = "gg_group2", measure = "Chao1")
      wilcox.test(Chao1 ~ treatment, data = alpha_d[alpha_d$timepoint == "112",]) # *
      
      # Shannon
      p2 <- graphs[[2]]+
        scale_x_continuous(
          breaks = seq(min(sample_data(ps)$timepoint), max(sample_data(ps)$timepoint), by = 7) # breaks every 7 days                    # "T" before each label
        ) +
        labs(y = "Index", x = "Days", color = NULL, fill = "", title = 'Shannon')+
        guides(shape = "none", fill = "none")+
        theme(axis.ticks.x = element_line(),
              legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5),
              axis.text.x = element_text(size = 6)
        )+
        coord_cartesian(ylim = c(0, 5))+
        geom_text(
          data = data.frame(
            x = c(54, 112),
            y = c(3.75, 3.5),
            label = c("*", "*")
          ),
          aes(x = x, y = y, label = label),
          size = 4,
          inherit.aes = FALSE)
      
      p2 <- p2+
        annotate("rect", xmin = 35, xmax = 49, # 50 ppm
                 fill = "#95BECF", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)+
        annotate("rect", xmin = 54, xmax = Inf, # 50 ppm
                 fill = "#95BECF", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)+
        annotate("rect", xmin = -Inf, xmax = 35, # Exposure
                 fill = "#F2AA84", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)+
        annotate("rect", xmin = 49, xmax = 54, # DSS
                 fill = "gray56", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)
      
      p2 <- putGgLastLayerBack(p2, nLayers = 4)
      
      p2 <- p2+
        scale_y_continuous(labels = scales::label_comma(big.mark = ","))+
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
      p2
      
      # Stats
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "0",], group = "gg_group2", measure = "Shannon")
      wilcox.test(Shannon ~ treatment, data = alpha_d[alpha_d$timepoint == "0",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "35",], group = "gg_group2", measure = "Shannon")
      wilcox.test(Shannon ~ treatment, data = alpha_d[alpha_d$timepoint == "35",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "gg_group2", measure = "Shannon")
      wilcox.test(Shannon ~ treatment, data = alpha_d[alpha_d$timepoint == "49",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "54",], group = "gg_group2", measure = "Shannon")
      wilcox.test(Shannon ~ treatment, data = alpha_d[alpha_d$timepoint == "54",]) # *
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "112",], group = "gg_group2", measure = "Shannon")
      wilcox.test(Shannon ~ treatment, data = alpha_d[alpha_d$timepoint == "112",]) # *
      
      # Inverse Simpson
      p3 <- graphs[[3]]+
        scale_x_continuous(
          breaks = seq(min(sample_data(ps)$timepoint), max(sample_data(ps)$timepoint), by = 7) # breaks every 7 days                  # "T" before each label
        ) +
        labs(y = "Index", x = "Days", color = NULL, fill = "", title = 'Inverse Simpson')+
        guides(shape = "none", fill = "none")+
        theme(axis.ticks.x = element_line(),
              legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5),
              axis.text.x = element_text(size = 6))+
        coord_cartesian(ylim = c(0, 28))+
        geom_text(
          data = data.frame(
            x = c(54, 112),
            y = c(15, 7),
            label = c("**", "*")
          ),
          aes(x = x, y = y, label = label),
          size = 4,
          inherit.aes = FALSE)
      
      p3 <- p3+
        annotate("rect", xmin = 35, xmax = 49, # 50 ppm
                 fill = "#95BECF", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)+
        annotate("rect", xmin = 54, xmax = Inf, # 50 ppm
                 fill = "#95BECF", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)+
        annotate("rect", xmin = -Inf, xmax = 35, # Exposure
                 fill = "#F2AA84", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)+
        annotate("rect", xmin = 49, xmax = 54, # DSS
                 fill = "gray56", alpha = 0.2,
                 ymin = -Inf, ymax = Inf)
      
      p3 <- putGgLastLayerBack(p3, nLayers = 4)
      
      p3 <- p3+
        scale_y_continuous(labels = scales::label_comma(big.mark = ","))+
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
      p3
      
      # Stats
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "0",], group = "gg_group2", measure = "InvSimpson")
      wilcox.test(InvSimpson ~ treatment, data = alpha_d[alpha_d$timepoint == "0",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "35",], group = "gg_group2", measure = "InvSimpson")
      wilcox.test(InvSimpson ~ treatment, data = alpha_d[alpha_d$timepoint == "35",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "gg_group2", measure = "InvSimpson")
      wilcox.test(InvSimpson ~ treatment, data = alpha_d[alpha_d$timepoint == "49",]) # n.s.
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "54",], group = "gg_group2", measure = "InvSimpson")
      wilcox.test(InvSimpson ~ treatment, data = alpha_d[alpha_d$timepoint == "54",]) # **
      
      verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "112",], group = "gg_group2", measure = "InvSimpson")
      wilcox.test(InvSimpson ~ treatment, data = alpha_d[alpha_d$timepoint == "112",]) # *
      
      alpha_d500 <- p1 / p2 /p3 +plot_layout(guides = "collect") 
    }
    
    # Boxplots at last timepoint
    {
      sample_data(ps_tfinal)$gg_group2 <- factor(sample_data(ps_tfinal)$gg_group2, levels = c(c("50:water","50:dss","500:water","500:dss")),
                                                 labels = c("Ctrl 50","DSS 50","Ctrl 500","DSS 500"))
      
      # Estinate richness measures for dataset
      richness_data <- estimate_richness(ps_tfinal, measures = c("Chao1", "Shannon", "InvSimpson"))
      alpha_d <- cbind(as.data.frame(sample_data(ps_tfinal)), richness_data)
      
      graphs <- alphaDiversityBoxplot(ps = ps_tfinal, group = "gg_group2", shape = "treatment", custom_colors = c("#95BECF","#325BAD","#F2AA84","#B22222"))
      
      # Chao1 
      p1 <- graphs[[1]]+
        geom_signif(
          comparisons = list(c("Ctrl 50","DSS 50"),c("Ctrl 500","DSS 500")),
          annotations = c("*"),
          y_position = c(1300, 1300),
          tip_length = 0.03,
          color = "black",
          size = 0.4,
          textsize = 5,
          margin_top = 0.1)+
        guides(x = legendry::guide_axis_nested())+
        geom_vline(xintercept = 2.5,
                   linetype    = "dashed",
                   colour      = "black",
                   size        = 0.5)+
        scale_y_continuous(labels = scales::label_comma(big.mark = ","), limits = c(0,1550))
        
      # Stats
      verifyStatsAssumptions(df = alpha_d, group = "gg_group2", measure = "Chao1")
      TukeyHSD(aov(Chao1 ~ gg_group2 , data = alpha_d))
      
      # Shannon
      p2 <- graphs[[2]]+
        ylim(0,5)+
        geom_signif(
          comparisons = list(c("Ctrl 50","DSS 50")),
          annotations = c("n.s."),
          y_position = c(4),
          tip_length = 0.03,
          color = "black",
          size = 0.4,
          textsize = 3,
          margin_top = 0.1)+
        geom_signif(
          comparisons = list(c("Ctrl 500","DSS 500")),
          annotations = c("*"),
          y_position = c(4),
          tip_length = 0.03,
          color = "black",
          size = 0.4,
          textsize = 5,
          margin_top = 0.1)+
        guides(x = legendry::guide_axis_nested())+
        geom_vline(xintercept = 2.5,
                   linetype    = "dashed",
                   colour      = "black",
                   size        = 0.5)
      
      # Stats
      verifyStatsAssumptions(df = alpha_d, group = "gg_group2", measure = "Chao1")
      TukeyHSD(aov(Chao1 ~ gg_group2 , data = alpha_d))
      
      # Inverse Simpson
      p3 <- graphs[[3]]+
        ylim(0,18)+
        geom_signif(
          comparisons = list(c("Ctrl 50","DSS 50")),
          annotations = c("n.s."),
          y_position = c(15),
          tip_length = 0.03,
          color = "black",
          size = 0.4,
          textsize = 3,
          margin_top = 0.1)+
        geom_signif(
          comparisons = list(c("Ctrl 500","DSS 500")),
          annotations = c("*"),
          y_position = c(10),
          tip_length = 0.03,
          color = "black",
          size = 0.4,
          textsize = 5,
          margin_top = 0.1)+
        guides(x = legendry::guide_axis_nested())+
        geom_vline(xintercept = 2.5,
                   linetype    = "dashed",
                   colour      = "black",
                   size        = 0.5)
      
      tfinal_stack <- p1/p2/p3
      
    }
    
    # Combining the two
    alpha_d <- (alpha_d50 | alpha_d500 | tfinal_stack)+
      plot_layout(widths = c(2, 2, 1))
    alpha_d
  
  }
  
}

# Beta diversity
{
  # Generate graphs
  {
    set.seed(100)
    
    # For diet at t0 t35 and t49
    phy_tree(ps_flt_diet) <- midpoint(tree) # Add rooted tree to phyloseq object
    
    # For diet timepoints
    sample_data(ps_flt_diet)$diet <- factor(sample_data(ps_flt_diet)$diet, labels = c("50 ppm","500 ppm"))
    
    # Weighted unifrac
    diet_beta_d_plot <- betaDiversityTimepoint2Factors(ps_flt_diet, sample_id = "sample_id", timeVariable = "timepoint",
                                                       varToCompare =  "diet", shape = "treatment", distMethod ="wunifrac",
                                                       transform = "rel_ab", customColors = c("#95BECF","#F2AA84"),
                                                       font = "Arial", path = "path/to/save/graph", 
                                                       additionnalAes = my_theme()+theme(plot.title = element_text(size = 10),
                                                                                         legend.position = "none"), dim = c(4,12), displayPValue = FALSE,
                                                       combineGraphs = TRUE, returnFig = TRUE, customTitles = 
                                                         c("Weaning","End of diets","2 weeks of\niron sufficent diet"),
                                                       positionPvalue = "right")
    
    # For diet + treatment at t54 and tfinal
    phy_tree(ps_flt_dss) <- midpoint(tree) # Add rooted tree to phyloseq object
    sample_data(ps_flt_dss)$gg_group2 <- factor(sample_data(ps_flt_dss)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl", "50 ppm DSS", "500 ppm DSS"))
    # Weighted unifrac => not very informative
    diet_dss_beta_d_plot <- betaDiversityTimepoint2Factors(ps_flt_dss, sample_id = "sample_id", timeVariable = "timepoint",
                                                           varToCompare =  "gg_group2", distMethod ="wunifrac", shape = "treatment",
                                                           customColors = c("#95BECF","#F2AA84","#325BAD","#B22222"),
                                                           font = "Arial", path = "path/to/save/graph",
                                                           additionnalAes = my_theme()+theme(plot.title = element_text(size = 10)), dim = c(4,8), combineGraphs = TRUE, returnFig = TRUE,
                                                           customTitles = c("End of DSS treatment","End of experiment"), pairwiseAdonis = FALSE, displayPValue = FALSE, hideLegend = TRUE)
    
    # Compute pvalues of pairwise comparisons with FDR correction
    {
      # For diet + treatment at t54 and tfinal
      phy_tree(ps_flt_dss) <- midpoint(tree) # Add rooted tree to phyloseq object
      sample_data(ps_flt_dss)$gg_group2 <- factor(sample_data(ps_flt_dss)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl", "50 ppm DSS", "500 ppm DSS"))
      # Weighted unifrac => not very informative
      betaDiversityTimepoint2Factors(ps_flt_dss, sample_id = "sample_id", timeVariable = "timepoint",
                                     varToCompare =  "gg_group2", distMethod ="wunifrac", shape = "treatment",
                                     customColors = c("#95BECF","#F2AA84","#325BAD","#B22222"),
                                     font = "Arial", path = "path/to/save/graph",
                                     additionnalAes = my_theme()+theme(plot.title = element_text(size = 10)), dim = c(4,12), combineGraphs = TRUE, returnFig = FALSE,
                                     customTitles = c("End of DSS treatment","End of experiment"), pairwiseAdonis = TRUE, displayPValue = FALSE, hideLegend = TRUE,printPairwiseSigtable = TRUE)
      
    }
    
    # Weighted unifrac for 50 ctrl vs 500 ctrl - last timepoint
    ps_sub <- prune_samples(sample_data(ps_tfinal_flt)$treatment == "water", ps_tfinal_flt)
    sample_data(ps_sub)$diet <- factor(sample_data(ps_sub)$diet, labels = c("50 ppm Ctrl","500 ppm Ctrl"))
    phy_tree(ps_sub) <- midpoint(tree) # Add rooted tree to phyloseq object
    Ctrl50vs500_beta_d_plot <- betaDiversityTimepoint2Factors(ps_sub , sample_id = "sample_id", timeVariable = "timepoint",
                                                              varToCompare =  "diet", distMethod ="wunifrac", shape = "treatment",
                                                              customColors = c("#95BECF","#F2AA84"),
                                                              font = "Arial", path = "path/to/save/graph",
                                                              additionnalAes = theme_beta_d(), dim = c(2.5,4), returnFig = TRUE, displayPValue = FALSE,
                                                              customTitles = c("50 ppm Ctrl\nvs 500 ppm Ctrl"), combineGraphs = FALSE, hideLegend = TRUE, pvalueToDisplay = 0.268, subGraph = TRUE,
                                                              size = 1, stroke = 0.15)
    
    # Weighted unifrac for dss groups - last timepoint
    ps_sub <- prune_samples(sample_data(ps_tfinal_flt)$treatment == "dss", ps_tfinal_flt)
    sample_data(ps_sub)$diet <- factor(sample_data(ps_sub)$diet, labels = c("50 ppm DSS","500 ppm DSS"))
    phy_tree(ps_sub) <- midpoint(tree) # Add rooted tree to phyloseq object
    DSS_50v500_beta_d_plot <- betaDiversityTimepoint2Factors(ps_sub , sample_id = "sample_id", timeVariable = "timepoint",
                                                             varToCompare =  "diet", distMethod ="wunifrac", shape = "treatment",
                                                             customColors = c("#325BAD","#B22222"),
                                                             font = "Arial", path = "path/to/save/graph",
                                                             additionnalAes = theme_beta_d(), dim = c(2.5,4), returnFig = TRUE, displayPValue = FALSE,
                                                             customTitles = c("50 ppm DSS\nvs 500 ppm DSS"), combineGraphs = FALSE, hideLegend = TRUE, pvalueToDisplay = 0.015, selectShape = c(22), subGraph = TRUE,
                                                             size = 1, stroke = 0.15)
    
    # Weighted unifrac for 50 ctrl vs 50 dss - last timepoint
    ps_sub <- prune_samples(sample_data(ps_tfinal_flt)$diet == "50", ps_tfinal_flt)
    sample_data(ps_sub)$treatment <- factor(sample_data(ps_sub)$treatment, labels = c("50 ppm Ctrl","50 ppm DSS"))
    phy_tree(ps_sub) <- midpoint(tree) # Add rooted tree to phyloseq object
    Ctrl50vDSS50_beta_d_plot <- betaDiversityTimepoint2Factors(ps_sub , sample_id = "sample_id", timeVariable = "timepoint",
                                                               varToCompare =  "treatment", distMethod ="wunifrac", shape = "treatment",
                                                               customColors = c("#95BECF","#325BAD"),
                                                               font = "Arial", path = "path/to/save/graph",
                                                               additionnalAes = theme_beta_d(), dim = c(2.5,4), returnFig = TRUE, displayPValue = FALSE,
                                                               customTitles = c("50 ppm Ctrl\nvs 50 ppm DSS"), combineGraphs = FALSE, hideLegend = TRUE, pvalueToDisplay = 0.003, subGraph = TRUE, positionPvalue = "right",
                                                               size = 1, stroke = 0.15)
    
    # Weighted unifrac for 500 ctrl vs 500 dss - last timepoint
    ps_sub <- prune_samples(sample_data(ps_tfinal_flt)$diet == "500", ps_tfinal_flt)
    sample_data(ps_sub)$treatment <- factor(sample_data(ps_sub)$treatment, labels = c("500 ppm Ctrl","500 ppm DSS"))
    phy_tree(ps_sub) <- midpoint(tree) # Add rooted tree to phyloseq object
    Ctrl500vDSS500_beta_d_plot <- betaDiversityTimepoint2Factors(ps_sub , sample_id = "sample_id", timeVariable = "timepoint",
                                                                 varToCompare =  "treatment", distMethod ="wunifrac", shape = "treatment",
                                                                 customColors = c("#F2AA84","#B22222"),
                                                                 font = "Arial", path = "path/to/save/graph",
                                                                 additionnalAes =theme_beta_d(), dim = c(2.5,4), returnFig = TRUE, displayPValue = FALSE,
                                                                 customTitles = c("500 ppm Ctrl\nvs 500 ppm DSS"), combineGraphs = FALSE, hideLegend = TRUE, pvalueToDisplay = 0.048, subGraph = TRUE, positionPvalue = "right",
                                                                 size = 1, stroke = 0.15)
    
  }
  
  # Combine all graphs
  stack <- Ctrl50vs500_beta_d_plot+DSS_50v500_beta_d_plot + Ctrl50vDSS50_beta_d_plot+Ctrl500vDSS500_beta_d_plot &
    theme(plot.margin = margin(3, 8, 3, 8))
  main <- (diet_beta_d_plot | diet_dss_beta_d_plot)+
    plot_layout(widths = c(1.5,1))
  

}

# Stackbar figure at last day of DSS
{
  # Timepoints and groups must be ordered properly and as factors
  sample_data(ps_t54_flt)$diet 
  sample_data(ps_t54_flt)$treatment 
  sample_data(ps_t54_flt)$gg_group2
  sample_data(ps_t54_flt)$timepoint
  
  diet_dss_phyla_fam <- plot_microbiota_2Fac(
    ps_object = ps_t54_flt,
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
    legend_size = 8,
    fdr_threshold = 0.05,
    main_level = "Phylum",
    sub_level = "Family",
    n_phy = 4, # number of taxa to show 
    mult_comp = F, # pairwise comparaisons for diff ab analysis
    selected_comparisons = list(c("50:water", "500:water"),
                                c("50:dss", "500:dss"),
                                c("50:water", "50:dss"),
                                c("500:water", "500:dss")
                                ),
    showOnlySubLegend = FALSE
  )
  
  print(diet_dss_phyla_fam$plot)
  print(diet_dss_phyla_fam$significant_table_main)
  print(diet_dss_phyla_fam$significant_table_sub)
  
  # Prepare Mouse IDs to display for each facet
  facet_scales <- list(
    scale_x_discrete(labels = as.character(1:12)),
    scale_x_discrete(labels = as.character(13:24)),
    scale_x_discrete(labels = as.character(25:35)),
    scale_x_discrete(labels = as.character(36:45))
  )
  
  # Customize plot
  stackbar_t54 <- diet_dss_phyla_fam$plot + 
    facet_wrap2(~ gg_group2, 
                scales  = "free_x", nrow = 2, ncol = 2,
                strip = strip_themed(background_x = elem_list_rect(fill = c("#95BECF","#F2AA84","#325BAD","#B22222"))), 
                labeller = as_labeller(c("50:water" = "50 ppm Ctrl",
                                         "500:water" = "500 ppm Ctrl",
                                         "50:dss" = "50 ppm DSS",
                                         "500:dss" = "500 ppm DSS")))+
    theme(text = element_text(family = "Arial"),      # Global text settings
          strip.text = element_text(size = 14, face = "bold", color = "white"),  # Facet titles
          plot.title = element_text(size = 20, face = "bold"),  # Main title
          axis.title = element_text(size = 15, face = "bold"),  # Axis titles
          axis.text = element_text(size = 12, face = "bold"),   # Axis text
          axis.title.y = element_text(margin = margin(r = -15)),
          legend.title = element_text(face = "bold", size = 14)  # Legend title  # Legend text
    ) +
    facetted_pos_scales(x = facet_scales)+
    labs(x = "Mouse ID")
  stackbar_t54
  
  # Saving the stats data
  existingDirCheck("path/to/save/stats/table")
  writeStackbarExtendedSigTable(main_table = diet_dss_phyla_fam$significant_table_main, includeSubTable = TRUE, sub_table = diet_dss_phyla_fam$significant_table_sub, filepath = "path/to/save/stats/table")
  
  # pvalues heatmap for the main lvl stats
  pvalHmapPhyla <- pvaluesHmap(stats = as.data.frame(readxl::read_excel("path/to/save/stats/table")),
                               selected_comparisons = c("50:water_vs_500:water","50:dss_vs_500:dss","50:water_vs_50:dss", "500:water_vs_500:dss"), displayChangeArrows = TRUE, displayPValues = FALSE,
                               txn_lvl="Phylum", lvl = "main", taxons = diet_dss_phyla_fam$main_names, group = "gg_group2", path, verticalTilesSpacing = 0.9, lineWidth = 0.4)
  pvalHmapPhyla <- pvalHmapPhyla+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          plot.margin = margin(0, 0, 0, 0))+
    guides(fill = "none")
  pvalHmapPhyla
  
  # pvalues heatmap for the sub lvl stats
  pvalHmapFamily <- pvaluesHmap(stats = as.data.frame(readxl::read_excel("path/to/save/stats/table")),
                                selected_comparisons = c("50:water_vs_500:water","50:dss_vs_500:dss","50:water_vs_50:dss", "500:water_vs_500:dss"),
                                txn_lvl="Family", lvl = "sub", taxons =  diet_dss_phyla_fam$sub_names, group = "gg_group2", displayPValues = FALSE, displayChangeArrows = TRUE, path, verticalTilesSpacing = 0.9, lineWidth = 0.4) # You can add [!grepl("Others", x = iron_exp_family$sub_names)] to remove "others"
  pvalHmapFamily <- pvalHmapFamily+scale_x_discrete(labels = c("50 VS 500 / CTRL", "50 VS 500 / DSS", "CTRL VS DSS / 50", "CTRL VS DSS / 500"))+
    guides(fill = guide_legend(ncol = 2))+
    theme(text = element_text(family = "Arial"),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          plot.margin = margin(0, 0, 0, 0),
          legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5))
  pvalHmapFamily
  
  # Combine stats hmaps
  statHmap <- pvalHmapPhyla / pvalHmapFamily
}

# Relative abundance by phyla and family over the time course of the experiment
{
  # Generate graph
  {
    # Custom theme
    theme_chronobiome <- function() {
      theme_bw(base_size = 18) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black", face = "bold"),
          axis.title = element_text(size = 22, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA),
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(size = 10),
          panel.spacing = unit(0, "lines"),
          legend.title = element_text(face = "bold", size = 18),
          strip.text = element_text(face = "bold", color = "white", size = 22)
        )
    }
    
    sample_data(ps_flt_all)$gg_group2 <- factor(sample_data(ps_flt_all)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
    
    p <- plot_timeline_2_groups(
      ps_object = ps_flt_all,
      exp_group =  "gg_group2", # must be as factor
      time_group = "timepoint", # must be as factor
      sample_name = "sample_id",
      main_level = 'Phylum',
      sub_level = 'Family',
      threshold = 0.8,
      average_relab_per_group = TRUE,
      smoothing = FALSE,
      n_phy = 4,
      hues = c("Blues", "Greens", "Purples", "Oranges"),
      color_bias = 2,
      custom_theme = theme_chronobiome(),
    )
    p$main_fig
    
    chronobiome <- p$main_fig+
      facet_wrap2(~ gg_group2, 
                  scales  = "free_x", nrow = 2, ncol = 2,
                  strip = strip_themed(background_x = elem_list_rect(fill = c("#95BECF","#F2AA84","#325BAD","#B22222"))))+
      scale_x_continuous(breaks = seq(min(as.numeric(levels(sample_data(ps_flt_all)$timepoint))), max(as.numeric(levels(sample_data(ps_flt_all)$timepoint))), by = 7))+
      labs(x = "Days")+
      coord_cartesian(xlim = c(0, 112), ylim = c(0, 100), expand = FALSE)+
      theme(panel.spacing = unit(0.75, "cm"))
  }
  
  # Compute stats at last timepoint (with Deseq2) - at the Family level
  {
  ps_taxa <- tax_glom(ps_tfinal_flt, taxrank = "Family", NArm = FALSE)
  deseq <- phyloseq_to_deseq2(ps_taxa, ~ treatment+diet+diet:treatment)
  deseq <- DESeq(deseq, test="Wald", fitType = "parametric")
  print(resultsNames(deseq))
  
  #Partition results for specific pairwise comparaisons
  res_subset1 <- results(deseq, contrast = list(resultsNames(deseq)[3])) # 50vs500 ctrl
  sigtab_1 <- cbind(as(res_subset1, "data.frame"), as(tax_table(ps)[rownames(res_subset1), ], "matrix"))
  sigtab_1$comparaison <- 1
  
  res_subset2 <- results(deseq, contrast = list(c(resultsNames(deseq)[3], resultsNames(deseq)[4]))) # 50vs500 dss
  sigtab_2 <- cbind(as(res_subset2, "data.frame"), as(tax_table(ps)[rownames(res_subset2), ], "matrix"))
  sigtab_2$comparaison <- 2
  
  res_subset3 <- results(deseq, contrast = list(resultsNames(deseq)[2])) # 50 ctrl vs 50 dss
  sigtab_3 <- cbind(as(res_subset3, "data.frame"), as(tax_table(ps)[rownames(res_subset3), ], "matrix"))
  sigtab_3$comparaison <- 3
  
  res_subset4 <- results(deseq, contrast = list(c(resultsNames(deseq)[2], resultsNames(deseq)[4]))) # 500 ctrl vs 500 dss
  sigtab_4 <- cbind(as(res_subset4, "data.frame"), as(tax_table(ps)[rownames(res_subset4), ], "matrix"))
  sigtab_4$comparaison <- 4
  
  #Append the sigtabs together
  sigtab <- bind_rows(sigtab_1, sigtab_2, sigtab_3, sigtab_4)
  
  #Replacing NA padj by 1 (they correspond to this anyways)
  sigtab$padj[is.na(sigtab$padj)] <- 1
  
  #Add column that adds symbols for the significance 
  # Define significance levels
  sigtab$significance <- as.character(cut(sigtab[,"padj"],
                                          breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                          labels = c("***", "**", "*", "NS")))
  
  comparisonsNames <- c("50:water_vs_500:water","50:dss_vs_500:dss","50:water_vs_50:dss", "500:water_vs_500:dss")
  sig_df_list <- list() 
  for(vs in unique(sigtab$comparaison)){
    
    df <- sigtab[sigtab$comparaison == vs,]
    rownames(df) <- gsub("\\..*", "", rownames(df))
    sig_df_list[[vs]] <- df
  }
  names(sig_df_list) <- comparisonsNames
  sig_df_list_family <- sig_df_list
  }
  
  # Compute stats at last timepoint (with Deseq2) - at the Phylum level
  {
  ps_taxa <- tax_glom(ps_tfinal_flt, taxrank = "Phylum", NArm = FALSE)
  deseq <- phyloseq_to_deseq2(ps_taxa, ~ treatment+diet+diet:treatment)
  deseq <- DESeq(deseq, test="Wald", fitType = "parametric")
  print(resultsNames(deseq))
  
  #Partition results for specific pairwise comparaisons
  res_subset1 <- results(deseq, contrast = list(resultsNames(deseq)[3])) #wt putrescine vs vehicle // 50vs500 ctrl
  sigtab_1 <- cbind(as(res_subset1, "data.frame"), as(tax_table(ps)[rownames(res_subset1), ], "matrix"))
  sigtab_1$comparaison <- 1
  
  res_subset2 <- results(deseq, contrast = list(c(resultsNames(deseq)[3], resultsNames(deseq)[4]))) #il22 ko putrescine vs vehicle // 50vs500 dss
  sigtab_2 <- cbind(as(res_subset2, "data.frame"), as(tax_table(ps)[rownames(res_subset2), ], "matrix"))
  sigtab_2$comparaison <- 2
  
  res_subset3 <- results(deseq, contrast = list(resultsNames(deseq)[2])) #vehicle wt vs il22 ko // 50 ctrl vs 50 dss
  sigtab_3 <- cbind(as(res_subset3, "data.frame"), as(tax_table(ps)[rownames(res_subset3), ], "matrix"))
  sigtab_3$comparaison <- 3
  
  res_subset4 <- results(deseq, contrast = list(c(resultsNames(deseq)[2], resultsNames(deseq)[4]))) #putrescine wt vs il22 ko // 500 ctrl vs 500 dss
  sigtab_4 <- cbind(as(res_subset4, "data.frame"), as(tax_table(ps)[rownames(res_subset4), ], "matrix"))
  sigtab_4$comparaison <- 4
  
  #Append the sigtabs together
  sigtab <- bind_rows(sigtab_1, sigtab_2, sigtab_3, sigtab_4) #sigtab_interaction
  
  #Replacing NA padj by 1 (they correspond to this anyways)
  sigtab$padj[is.na(sigtab$padj)] <- 1
  
  #Add column that adds symbols for the significance 
  # Define significance levels
  sigtab$significance <- as.character(cut(sigtab[,"padj"],
                                          breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                          labels = c("***", "**", "*", "NS")))
  
  comparisonsNames <- c("50:water_vs_500:water","50:dss_vs_500:dss","50:water_vs_50:dss", "500:water_vs_500:dss")
  sig_df_list <- list() 
  for(vs in unique(sigtab$comparaison)){
    
    df <- sigtab[sigtab$comparaison == vs,]
    rownames(df) <- gsub("\\..*", "", rownames(df))
    sig_df_list[[vs]] <- df
  }
  names(sig_df_list) <- comparisonsNames
  sig_df_list_phylum <- sig_df_list
  }
  
  # Saving the plot and the associated stats
  existingDirCheck("path/to/save/stats/table")
  writeStackbarExtendedSigTable(main_table =  sig_df_list_phylum, includeSubTable = TRUE, sub_table =  sig_df_list_family, filepath = "path/to/save/stats/table")
  
  # pvalues heatmap for the main lvl stats
  pvalHmapPhyla <- pvaluesHmap(stats = as.data.frame(readxl::read_excel("path/to/save/stats/table")),
                               selected_comparisons = c("50:water_vs_500:water","50:dss_vs_500:dss","50:water_vs_50:dss", "500:water_vs_500:dss"), displayChangeArrows = TRUE, displayPValues = FALSE,
                               txn_lvl="Phylum", lvl = "main", taxons = p$main_names, group = "gg_group2", path, verticalTilesSpacing = 0.9, lineWidth = 0.4)
  pvalHmapPhyla <- pvalHmapPhyla+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          plot.margin = margin(0, 0, 0, 0))+
    guides(fill = "none")
  pvalHmapPhyla
  
  # pvalues heatmap for the sub lvl stats
  pvalHmapFamily <- pvaluesHmap(stats = as.data.frame(readxl::read_excel("path/to/save/stats/table")),
                                selected_comparisons = c("50:water_vs_500:water","50:dss_vs_500:dss","50:water_vs_50:dss", "500:water_vs_500:dss"),
                                txn_lvl="Family", lvl = "sub", taxons = p$sub_names, group = "gg_group2", displayPValues = FALSE, displayChangeArrows = TRUE, path, verticalTilesSpacing = 0.9, lineWidth = 0.4) # You can add [!grepl("Others", x = iron_exp_family$sub_names)] to remove "others"
  pvalHmapFamily <- pvalHmapFamily+scale_x_discrete(labels = c("50 VS 500 / CTRL", "50 VS 500 / DSS", "CTRL VS DSS / 50", "CTRL VS DSS / 500"))+
    guides(fill = guide_legend(ncol = 2))+
    theme(text = element_text(family = "Arial"),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          plot.margin = margin(0, 0, 0, 0),
          legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5))
  pvalHmapFamily
  
  # Combine stats hmaps
  statHmap <- pvalHmapPhyla / pvalHmapFamily
  
  # Compute stats timepoint t0 t35 and t49 between 50 ppm and 500 ppm diets (DEseq2) - Family level 
  {
  sigtab_all_tp <- list()
  for(tp in levels(sample_data(ps_flt_diet)$timepoint)){
    ps_sub <-  prune_samples(sample_data(ps_flt_diet)$timepoint ==  tp, ps_flt_diet)
    ps_taxa <- tax_glom(ps_sub, taxrank = "Family", NArm = FALSE) 
    deseq <- phyloseq_to_deseq2(ps_taxa, ~ diet)
    deseq <- DESeq(deseq, test="Wald", fitType = "parametric")
    print(resultsNames(deseq))
    res <- results(deseq, contrast = list(resultsNames(deseq)[2])) #wt putrescine vs vehicle // 50vs500 ctrl
    sigtab <- cbind(as(res, "data.frame"), as(tax_table(ps_taxa)[rownames(res), ], "matrix"))
    sigtab <- sigtab[sigtab$Family %in% p$sub_names,]
    sigtab <- sigtab[c("log2FoldChange","padj","Family")]
    sigtab$timepoint <- tp
    sigtab_all_tp[[as.character(tp)]] <- sigtab
  }
    
    sigtab_all_tp <- dplyr::bind_rows(sigtab_all_tp)
    
    # Build arrow + p-value label 
    sigtab_all_tp <- sigtab_all_tp %>%
      mutate(
        arrow = if_else(log2FoldChange < 0, "↓",
                        if_else(log2FoldChange > 0, "↑", "↔")),
        p_lab = case_when(
          is.na(padj)      ~ "n.s.",
          padj < 0.001     ~ "P<0.001",
          padj < 0.01      ~ "P<0.01",
          padj < 0.05      ~ "P<0.05",
          TRUE             ~ "n.s."
        )
      )
    
    sigtab_all_tp$dir_p = ifelse(sigtab_all_tp$p_lab == "n.s.", sigtab_all_tp$p_lab, paste0(sigtab_all_tp$arrow, ", ", sigtab_all_tp$p_lab))
    
    #  Order by levels(p$sub_names) 
    levs <- if (is.factor(p$sub_names)) levels(p$sub_names) else as.character(p$sub_names)
    
    sigtab_all_tp <- sigtab_all_tp %>%
      mutate(Family = factor(Family, levels = levs)) %>%
      arrange(Family, timepoint)
    
    # Respect your custom Family order and the timepoint order from the ps object
    fam_levels <- if (is.factor(p$sub_names)) levels(p$sub_names) else as.character(p$sub_names)
    tp_levels  <- levels(sample_data(ps_flt_diet)$timepoint)
    
    out_wide <- sigtab_all_tp %>%
      mutate(
        Family    = factor(Family, levels = fam_levels),
        timepoint = factor(timepoint, levels = tp_levels)
      ) %>%
      select(Family, timepoint, dir_p) %>%
      distinct(Family, timepoint, .keep_all = TRUE) %>%        # safety in case of duplicates
      pivot_wider(names_from = timepoint, values_from = dir_p) %>%
      arrange(Family)
    
    write.xlsx(out_wide, "~/CHUM_git/figures/memoire/dss/Annotated/supplementary/family_sig_diet_timepoints.xlsx")
    

  }
  
}

# Differentially abundant bacterial species at last timepoint
{
  # Volcano plot - 50 DSS VS 500 DSS last timepoint
  {
      deseq <- phyloseq_to_deseq2(ps_tfinal_flt, ~ treatment+diet+diet:treatment)
      deseq <- DESeq(deseq, test="Wald", fitType = "parametric")
      print(resultsNames(deseq))
      
      # Manually add species names that we identified with BLASTn (ASV 8,38,19,6 had no species assignments) 
      tax_table(ps_tfinal_flt)["ASV8","Species"] <- "murinus (ASV8)"
      tax_table(ps_tfinal_flt)["ASV38","Species"] <- "murinus (ASV38)"
      tax_table(ps_tfinal_flt)["ASV19","Species"] <- "intestinale"
      tax_table(ps_tfinal_flt)["ASV6","Species"] <- "rodentium (ASV6)"
      tax_table(ps_tfinal_flt)["ASV1","Species"] <- "rodentium (ASV1)"

      # Volcano plot
      volcanoPlot50vs500dss <- volcanoPlot2GroupsMultifactorDesign(ps = ps_tfinal_flt, deseq = deseq, varToCompare = "diet",
                                          taxa = "Species", threshold = 0.05, FCcutoff = 0.49, customColors = NULL,
                                          FDR = TRUE, includeUnknownSpecies = TRUE, selectedComparison = 2,
                                          title = "50 ppm DSS VS 500 ppm DSS\nat end of experiment")+
        theme(
            plot.title = element_text(size = 16, face = "bold"),
            axis.text.x = element_text(color = "black", face = "bold", size = 12),
            axis.text.y = element_text(color = "black", face = "bold", size = 12),
            axis.title.x = element_text(size = 13, face = "bold"),
            axis.title.y = element_text(size = 13, face = "bold"),
            axis.line = element_line(color = "black", size = 1),
            legend.text = element_text(face = "bold", size = 8, margin = margin(l = 6, unit = "pt")),
            legend.title = element_text(face = "bold", size = 10),
            legend.spacing.y = unit(0.1, 'cm'),
            legend.spacing.x = unit(0.1, 'cm'),
            legend.position = "right",
            legend.margin     = margin(3, 3, 3, 3),
            legend.box.spacing = unit(0.3, "cm"),
            legend.key.width = unit(0.05, "lines"),
            legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x = element_line()
          )+ 
        guides(
            color = guide_legend(ncol = 1, byrow = TRUE),
            fill  = guide_legend(ncol = 1, byrow = TRUE)  # if applicable
          )+
        labs(color = "Significance")+
        scale_color_manual(
          breaks = c("n.s.", "Unknown species", "Up", "Down"),
          values = c('grey','black','#B22222','#325BAD'))
      volcanoPlot50vs500dss
    
  }
  
  # Individual graphs for each species displayed as full timeline - for each diet group
  {
    sample_data(ps_flt_all)$gg_group2 <- factor(sample_data(ps_flt_all)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
    sample_data(ps_flt_all)$timepoint <- as.numeric(as.character(sample_data(ps_flt_all)$timepoint))
    
    # Manually add species names that we identified with BLASTn
    tax_table(ps_flt_all)["ASV8","Species"] <- "murinus (ASV8)"
    tax_table(ps_flt_all)["ASV38","Species"] <- "murinus (ASV38)"
    tax_table(ps_flt_all)["ASV19","Species"] <- "intestinale"
    tax_table(ps_flt_all)["ASV6","Species"] <- "rodentium (ASV6)"
    tax_table(ps_flt_all)["ASV1","Species"] <- "rodentium (ASV1)"
    
    # Custom theme
    timeline_theme <- theme(
      axis.text.x    = element_text(size = 7),
      axis.title.x   = element_text(size = 10, margin = margin(t = 10)),
      axis.title.y   = element_text(size = 8, margin = margin(r = 10)),
      plot.title     = element_text(size = 10, margin = margin(b = 10)),
      plot.margin    = margin(5, 5, 5, 5),
      legend.position = "none"
    )
    
    ps_sub50 <- prune_samples(sample_data(ps_flt_all)$diet == "50", ps_flt_all)
    ps_sub500 <- prune_samples(sample_data(ps_flt_all)$diet == "500", ps_flt_all)
    
    # Function that produces graph for full timeline asv, adapts background depending on 50 or 500
    asvRelAbDistributionTimelineExtension <- function(ps, asv, taxon, group, time, custom_colors, diet){
      
      p <- asvRelAbDistributionTimeline(ps = ps, asv = asv, taxon = taxon,
                                    group = group, time = time,custom_colors = custom_colors ,
                                   displayASVNumber = FALSE)+
        labs(color = NULL, x = "Days")+
        scale_x_continuous(
          breaks = seq(min(sample_data(ps)[[time]]), max(sample_data(ps)[[time]]), by = 7)                   # "T" before each label
        )+
        timeline_theme
      
      if(diet == "50"){
        p <- p+
          annotate("rect", xmin = -Inf, xmax = 49, # 50 ppm
                   fill = "#95BECF", alpha = 0.2,
                   ymin = -Inf, ymax = Inf)+
          annotate("rect", xmin = 54, xmax = Inf, # 50 ppm
                   fill = "#95BECF", alpha = 0.2,
                   ymin = -Inf, ymax = Inf)+
          annotate("rect", xmin = 49, xmax = 54, # DSS
                   fill = "gray56", alpha = 0.2,
                   ymin = -Inf, ymax = Inf)
        p <- putGgLastLayerBack(p, nLayers = 3)
      }else if(diet == "500"){
        p <- p+
          annotate("rect", xmin = 35, xmax = 49, # 50 ppm
                   fill = "#95BECF", alpha = 0.2,
                   ymin = -Inf, ymax = Inf)+
          annotate("rect", xmin = 54, xmax = Inf, # 50 ppm
                   fill = "#95BECF", alpha = 0.2,
                   ymin = -Inf, ymax = Inf)+
          annotate("rect", xmin = -Inf, xmax = 35, # 500 ppm
                   fill = "#F2AA84", alpha = 0.2,
                   ymin = -Inf, ymax = Inf)+
          annotate("rect", xmin = 49, xmax = 54, # DSS
                   fill = "gray56", alpha = 0.2,
                   ymin = -Inf, ymax = Inf)
        p <- putGgLastLayerBack(p, nLayers = 4)
      }

      return(p)
      
    }
    
    # For 50
    # Species M. intestinale ASV19
    asv19 <- asvRelAbDistributionTimelineExtension(ps_sub50, "ASV19", taxon = "Species",
                                          "gg_group2", "timepoint", c("#95BECF","#325BAD"),
                                          diet = "50")
    asv19
    
    # Species ligilactobacillus murinus ASV 8
    asv8 <- asvRelAbDistributionTimelineExtension(ps_sub50, "ASV8", taxon = "Species",
                                                  "gg_group2", "timepoint", c("#95BECF","#325BAD"),
                                                  diet = "50")
    
    # Species ligilactobacillus murinus ASV 38
    asv38 <- asvRelAbDistributionTimelineExtension(ps_sub50, "ASV38", taxon = "Species",
                                                   "gg_group2", "timepoint", c("#95BECF","#325BAD"),
                                                   diet = "50")
    
    # F rodentium - ASV1
    asv1 <- asvRelAbDistributionTimelineExtension(ps_sub50, "ASV1", taxon = "Species",
                                                  "gg_group2", "timepoint", c("#95BECF","#325BAD"),
                                                  diet = "50")
    
    # F rodentium - ASV6
    asv6 <- asvRelAbDistributionTimelineExtension(ps_sub50, "ASV6", taxon = "Species",
                                                  "gg_group2", "timepoint", c("#95BECF","#325BAD"),
                                                  diet = "50")
    asv6
    
    
    # Combine individual graphs into one with merged legend
    species_relab_timeline50 <- (asv1 / asv6 / asv19 / asv8 / asv38)+
      plot_layout(guides = "collect")
    species_relab_timeline50
    
    
    # For 500
    # Species M. intestinale ASV19
    asv19 <- asvRelAbDistributionTimelineExtension(ps_sub500, "ASV19", taxon = "Species",
                                                   "gg_group2", "timepoint", c("#F2AA84","#B22222"),
                                                   diet = "500")
    asv19
    
    # Species ligilactobacillus murinus ASV 8
    asv8 <- asvRelAbDistributionTimelineExtension(ps_sub500, "ASV8", taxon = "Species",
                                                  "gg_group2", "timepoint", c("#F2AA84","#B22222"),
                                                  diet = "500")
    
    # Species ligilactobacillus murinus ASV 38
    asv38 <- asvRelAbDistributionTimelineExtension(ps_sub500, "ASV38", taxon = "Species",
                                                   "gg_group2", "timepoint", c("#F2AA84","#B22222"),
                                                   diet = "500")
    
    # F rodentium - ASV1
    asv1 <- asvRelAbDistributionTimelineExtension(ps_sub500, "ASV1", taxon = "Species",
                                                  "gg_group2", "timepoint", c("#F2AA84","#B22222"),
                                                  diet = "500")
    
    # F rodentium - ASV6
    asv6 <- asvRelAbDistributionTimelineExtension(ps_sub500, "ASV6", taxon = "Species",
                                                  "gg_group2", "timepoint", c("#F2AA84","#B22222"),
                                                  diet = "500")
    asv6
    
    
    # Combine individual graphs into one with merged legend
    species_relab_timeline500 <- (asv1 / asv6 / asv19 / asv8 / asv38)+
      plot_layout(guides = "collect")
    species_relab_timeline500
    
  }
  
  # Individual graphs for each species displayed for last timepoint
  {
    sample_data(ps_tfinal_flt)$gg_group2 <- factor(sample_data(ps_tfinal_flt)$gg_group2, labels = c("Ctrl 50","DSS 50","Ctrl 500","DSS 500"),
                                                   levels = c("50:water","50:dss","500:water","500:dss"))
    sample_data(ps_tfinal_flt)$treatment
    
    # Manually add species names that we identified with BLASTn
    tax_table(ps_tfinal_flt)["ASV8","Species"] <- "murinus (ASV8)"
    tax_table(ps_tfinal_flt)["ASV38","Species"] <- "murinus (ASV38)"
    tax_table(ps_tfinal_flt)["ASV19","Species"] <- "intestinale"
    tax_table(ps_tfinal_flt)["ASV6","Species"] <- "rodentium (ASV6)"
    tax_table(ps_tfinal_flt)["ASV1","Species"] <- "rodentium (ASV1)"
    
    # Species M. intestinale ASV19
    asv19 <- asvRelAbDistribution(ps_tfinal_flt, "ASV19",
                                  "gg_group2","treatment",
                                c("#95BECF","#325BAD","#F2AA84","#B22222"),
                                test_results = c("**","n.s.","n.s.","*"),
                                relativeAbundance = FALSE,
                                text_sizes = c(4,2,2,4), stats = TRUE, vjustList = c(0.5,0.05,0.05,0.5))+
      guides(x = legendry::guide_axis_nested())+
      ylim(0,NA)+
      theme(axis.title.y = element_text(size = 8))
    asv19

    
    # Species ligilactobacillus murinus ASV 8
    asv8 <- asvRelAbDistribution(ps_tfinal_flt, "ASV8",
                                 "gg_group2","treatment",
                                 c("#95BECF","#325BAD","#F2AA84","#B22222"),
                                 test_results = c("n.s.","n.s.","n.s.","*"),
                                 relativeAbundance = FALSE,
                                 text_sizes = c(2,2,2,4), stats =  TRUE, vjustList = c(0.05,0.05,0.05,0.5))+
      guides(x = legendry::guide_axis_nested())+
      ylim(0,NA)+
      theme(axis.title.y = element_text(size = 8))
    asv8
    
    # Species ligilactobacillus murinus ASV 38
    asv38 <- asvRelAbDistribution(ps_tfinal_flt, "ASV38",
                                  "gg_group2","treatment",
                                  c("#95BECF","#325BAD","#F2AA84","#B22222"),
                                  test_results = c("n.s.","n.s.","n.s.","*"),
                                  text_sizes = c(2,2,2,4), stats = TRUE, vjustList = c(0.05,0.05,0.05,0.5))+
      guides(x = legendry::guide_axis_nested())+
      ylim(0,NA)+
      theme(axis.title.y = element_text(size = 8))
    asv38
    
    # F rodentium - ASV1
    asv1 <- asvRelAbDistribution(ps_tfinal_flt, "ASV1",
                                 "gg_group2","treatment",
                                 c("#95BECF","#325BAD","#F2AA84","#B22222"),
                                 test_results = c("n.s.","n.s.","n.s.","*"),
                                 text_sizes = c(2,2,2,4), stats = TRUE, vjustList = c(0.05,0.05,0.05,0.5))+
      guides(x = legendry::guide_axis_nested())+
      ylim(0,NA)+
      theme(axis.title.y = element_text(size = 8))
    asv1
    
    # F rodentium - ASV6
    asv6 <- asvRelAbDistribution(ps_tfinal_flt, "ASV6",
                                 "gg_group2","diet",
                                 c("#95BECF","#325BAD","#F2AA84","#B22222"),
                                 test_results = c("n.s.","n.s.","n.s.","*"),
                                 text_sizes = c(2,2,2,4), stats = TRUE, vjustList = c(0.05,0.05,0.05,0.5))+
      guides(x = legendry::guide_axis_nested())+
      ylim(0,NA)+
      theme(axis.title.y = element_text(size = 8))
    asv6
    
    # Combine individual graphs into one with merged legend
    species_relab_time_final <- (asv1 / asv6 / asv19 / asv8 / asv38)
    species_relab_time_final
    
  }
  
  # Combine graphs and make final figure

  species_rel_ab_timeline <- (wrap_elements(full = species_relab_timeline50)+ wrap_elements(full = species_relab_timeline500) + wrap_elements(full = species_relab_time_final))+
    plot_layout(widths = c(2,2,1))
  
  
}

# Correlations heatmap between qPCRs results and variables
{
  # Load species qPCR measured abundance 
  df_pcr <- as.data.frame(read_excel("path/to/qPCR/data"))
  colnames(df_pcr) <- df_pcr[3,]
  df_pcr <- df_pcr[!is.na(df_pcr$`Muribaculum intestinale`),]
  df_pcr <- df_pcr[-1,]
  df_pcr$lm <- as.numeric(df_pcr$`Lm/16S`)
  df_pcr$mi <- as.numeric(df_pcr$`Mi/16S`)
  df_pcr$fr <- as.numeric(df_pcr$`Fr/16S`)
  df_pcr <- df_pcr[,-c(5,10,15)]
  row.names(df_pcr) <- df_pcr$`Mouse ID`
  df_pcr <- df_pcr[,c(17:19)]
  
  # Load variables 
  df_var <- as.data.frame(read_excel("path/to/variables/data"))
  df_var$id <- substring(df_var$id, first = 1, last = 5)
  df_var <- df_var[df_var$exp == "dss",]
  df_var <- df_var[-46,]
  df_var <- df_var[df_var$treatment == "dss",]
  df_var$diet <- factor(df_var$diet, levels = c(50,500))
  rownames(df_var) <- df_var$id
  df_var <- df_var[,-c(1:7,9,12:14,16)]
  colnames(df_var) <- c("Stool iron at\nend of exposure","Spleen iron","Liver iron","DAI")
  
  # Extract relative abundance of bacterial families of interest
  {
    ps_counts = transformCounts(tax_glom(ps_t35_flt, taxrank = "Family", NArm = FALSE), transformation = "rel_ab")
    counts <- t(otu_table(ps_counts))[,c("ASV2","ASV3","ASV9")]
    row.names(counts) <- str_replace(row.names(counts), pattern = "_T35", replacement = "")
    df_var<- merge(df_var, counts, by = "row.names")
    rownames(df_var) <- df_var$Row.names
    df_var <- df_var[,-1]
  }
  
  # Merge dataframes and perform spearman correlations
  df_cor <- merge(df_pcr, df_var, by = "row.names")
  row.names(df_cor) <- df_cor[,1] 
  df_cor <- df_cor[,-1]
  cor_results <- rcorr(as.matrix(df_cor), type = "spearman")
  cor_matrix <- cor_results$r  # Extract correlation coefficients
  cor_matrix <- cor_matrix[4:10,1:3]
  p_values <- cor_results$P    # Extract p-values
  p_values <- p_values[4:10,1:3]
  
  # Add LCN2 data
  df_var <- as.data.frame(read_excel("path/to/LCN2/data"))
  df_var$id <- substring(df_var$id, first = 1, last = 5)
  df_var <- df_var[!is.na(df_var$lcn2),]
  row.names(df_var) <- df_var[,1]
  df_var <- df_var[,-1]
  df_var <- df_var["lcn2"]
  colnames(df_var) <- "LCN2"

  df_cor <- merge(df_pcr, df_var, by = "row.names")
  row.names(df_cor) <- df_cor[,1] 
  df_cor <- df_cor[,-1]
  
  cor_results_lcn2 <- rcorr(as.matrix(df_cor), type = "spearman")
  cor_matrix_lcn2 <- cor_results_lcn2$r  # Extract correlation coefficients
  cor_matrix_lcn2 <- cor_matrix_lcn2[4,1:3]
  p_values_lcn2 <- cor_results_lcn2$P    # Extract p-values
  p_values_lcn2 <- p_values_lcn2[4,1:3]
  
  cor_matrix <- rbind(cor_matrix, cor_matrix_lcn2)
  p_values <- rbind(p_values, p_values_lcn2)
  row.names(cor_matrix)[8] <- "LCN2"
  row.names(p_values)[8] <- "LCN2"
  
  cor_matrix <- cor_matrix[c(1,2,3,4,8,5,6,7),]

  p_adjusted <- matrix(
    p.adjust(as.vector(p_values), method = "fdr"),
    nrow = nrow(p_values),
    ncol = ncol(p_values),
    dimnames = dimnames(p_values))
  
  p_adjusted <-  p_adjusted[c(1,2,3,4,8,5,6,7),]
  
  cor_melt <- melt(cor_matrix, varnames = c("Variables", "Species"))
  p_melt <- melt(p_adjusted, varnames = c("Variables","Species"))
  cor_melt$p_value <- p_melt$value # Combine melted data
  cor_melt <- cor_melt[!is.na(cor_melt$value), ]  # Remove NA rows
  cor_melt$significance <- ifelse(cor_melt$p_value < 0.001, "***", 
                                  ifelse(cor_melt$p_value < 0.01, "**", 
                                         ifelse(cor_melt$p_value < 0.05, "*", "")))
  
  correlation_hmap <- ggplot(cor_melt, aes(x = Variables, y = Species, fill = value)) + 
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#d1762d", high = "#639381", mid = "white",
                         midpoint = 0, limit = c(-1, 1), space = "Lab",
                         name = "Spearman's\nCorrelation") +
    geom_text(aes(label = significance), color = "black") +
    scale_y_discrete(labels  = c("L. murinus\nabundance\nat T112","M. intestinale\nabundance\nat T112","F. rodentium\nabundance\nat T112"))+
    scale_x_discrete(labels  = c("Stool iron at\nend of exposure","Spleen\niron","Liver\niron",
                                 "DAI","Fecal\nLCN2","Lactobacillaceae\nabundance\nat T35",
                                 "Muribaculaceae\nabundance\nat T35","Bacteroidaceae\nabundance\nat T35"))+
    my_theme()+
    theme(axis.line = element_blank(),
          legend.background = element_blank())
  
  correlation_hmap
  
}

# Supplementary figures
{
  # B. acidifaciens abundance at t49 and t54 (ASV70)
  {
  sample_data(ps_t49_flt)$diet <- factor(sample_data(ps_t49_flt)$diet, labels = c("50 ppm", "500 ppm"))
  sample_data(ps_t49_flt)$treatment
  
  # t54
  asv70_t54 <- asvRelAbDistribution(ps_t54_flt, "ASV70",displayTitle = TRUE,displayAsvNumber = TRUE,
                                    "diet","treatment",
                                    c("#95BECF","#325BAD","#F2AA84","#B22222"),
                                    relativeAbundance = TRUE, stats = TRUE,
                                    test_results = c("***","n.s.","n.s.","P = 0.14"),
                                    text_sizes = c(4,2,2,3), vjustList = c(0.5,0.05,0.05,0.05))+
    ylim(0,NA)+
    guides(x = legendry::guide_axis_nested())+
    theme(axis.title.y = element_text(size = 8))
  asv70_t54
  
  asv53_t54 <- asvRelAbDistribution(ps_t54_flt, "ASV53",displayTitle = TRUE,displayAsvNumber = TRUE,
                                    "diet","treatment",
                                    c("#95BECF","#325BAD","#F2AA84","#B22222"),
                                    relativeAbundance = TRUE, stats = TRUE,
                                    test_results = c("***","n.s.","n.s.","P = 0.18"),
                                    text_sizes = c(4,2,2,3), vjustList = c(0.5,0.05,0.05,0.05))+
    ylim(0,NA)+
    guides(x = legendry::guide_axis_nested())+
    theme(axis.title.y = element_text(size = 8))
  asv53_t54
  
  baT54 <- (asv70_t54 | asv53_t54)
  }
  
  # Relative abundance of species genera over time course of the experiment (For genus Bacteroides and family Enterobacteriaceae)
  {
    
    theme_chronobiome <- function() {
      theme_bw(base_size = 18) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black", face = "bold"),
          axis.title = element_text(size = 22, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA),
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(size = 10),
          panel.spacing = unit(0, "lines"),
          legend.title = element_text(face = "bold", size = 18),
          strip.text = element_text(face = "bold", color = "white", size = 22)
        )}
    
    sample_data(ps_flt_all)$gg_group2 <- factor(sample_data(ps_flt_all)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
    View(tax_table(tax_glom(ps_flt_all,taxrank = "Genus",NArm = FALSE)))
    ps_chrb <- ps_flt_diet
    sample_data(ps_chrb)$diet <- factor(sample_data(ps_chrb)$diet, labels = c("50 ppm","500 ppm"))
    tax_table(ps_chrb)[,"Species"] <- ifelse(is.na(tax_table(ps_chrb)[,"Species"]),
                                             paste0(tax_table(ps_chrb)[,"Genus"]," unkwown", " (",row.names(tax_table(ps_chrb)),")"),
                                             paste0(tax_table(ps_chrb)[,"Genus"], " ",tax_table(ps_chrb)[,"Species"], " (", row.names(tax_table(ps_chrb)),")"))
    # For b acidifaciens 
    p <- plot_timeline_selected_taxa(
      ps_object = ps_chrb,
      exp_group =  "diet", # must be as factor
      time_group = "timepoint", # must be as factor
      sample_name = "sample_id",
      main_level = 'Family',
      sub_level = 'Species',
      selected_taxa = "Enterobacteriaceae",
      threshold = 0.5,
      average_relab_per_group = TRUE,
      smoothing = FALSE,
      n_phy = 4,
      hues = c("Greens"),
      color_bias = 1,
      custom_theme = theme_chronobiome()
    )
    p$main_fig
    
    b <- p$main_fig+
      facet_wrap2(~ diet, 
                  scales  = "free_x", nrow = 2, ncol = 1,
                  strip = strip_themed(background_x = elem_list_rect(fill = c("#95BECF","#F2AA84"))))+
      scale_x_continuous(breaks = seq(min(as.numeric(levels(sample_data(ps_chrb)$timepoint))), max(as.numeric(levels(sample_data(ps_chrb)$timepoint))), by = 7))+
      labs(x = "Days")
    
    
    b <- b+
      coord_cartesian(xlim = c(0, 49), ylim = c(0, 3.5), expand = FALSE)+
      theme(panel.spacing = unit(0.75, "cm"))
    b 
    
    
    # For genus Eschiera Shigella 
    sample_data(ps_flt_all)$gg_group2 <- factor(sample_data(ps_flt_all)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
    ps_chrb <- ps_flt_all
    tax_table(ps_chrb)[,"Species"] <- ifelse(is.na(tax_table(ps_chrb)[,"Species"]),
                                             paste0(tax_table(ps_chrb)[,"Genus"]," unkwown", " (",row.names(tax_table(ps_chrb)),")"),
                                             paste0(tax_table(ps_chrb)[,"Genus"], " ",tax_table(ps_chrb)[,"Species"], " (", row.names(tax_table(ps_chrb)),")"))
    
    p <- plot_timeline_selected_taxa(
      ps_object = ps_chrb,
      exp_group =  "gg_group2", # must be as factor
      time_group = "timepoint", # must be as factor
      sample_name = "sample_id",
      main_level = 'Family',
      sub_level = 'Genus',
      selected_taxa = "Enterobacteriaceae",
      threshold = 0.01,
      average_relab_per_group = TRUE,
      smoothing = FALSE,
      n_phy = 4,
      hues = c("Oranges"),
      color_bias = 1,
      custom_theme = theme_chronobiome()
    )
    p$main_fig
    
    e <- p$main_fig+
      facet_wrap2(~ gg_group2, 
                  scales  = "free_x", nrow = 2, ncol = 2,
                  strip = strip_themed(background_x = elem_list_rect(fill = c("#95BECF","#F2AA84","#325BAD","#B22222"))))+
      scale_x_continuous(breaks = seq(min(as.numeric(levels(sample_data(ps_chrb)$timepoint))), max(as.numeric(levels(sample_data(ps_chrb)$timepoint))), by = 7))+
      labs(x = "Days")
    
    
    e <- e+
      coord_cartesian(xlim = c(0, 112), ylim = c(0, 6), expand = FALSE)+
      theme(panel.spacing = unit(0.75, "cm"))
    e 
    
    
    
  }
}
