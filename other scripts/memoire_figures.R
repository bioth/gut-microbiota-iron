#### Biological data related figures ####
# loading libraries
{
  library("tidyverse") # loading bunch of packages
  library("ggplot2") # come on, everyone knows what it is used for
  library("dplyr") # arranging and manipulating data easily
  library("lme4") # library for loading ANOVA
  library("car") # for anova too
  library("ggsignif") # adding significance bars to ggplots
  library("readxl") # Read and write excel files
  library(permuco)
  library(patchwork)
}

# Setting working directory
setwd("~/CHUM_git/gut-microbiota-iron/")

# Loading functions for data manipulation
source("other scripts/dataManipFunctions.R")
source("pipeline/microbiota_analysis/utilities.R")

# DSS project
# Figure 1 - iron measurements in young mice 
{
  # Set working directory
  setwd("experiments/finished exp/young-DSS-exp3")
  
  #T35 ferrozine assay for stools
  df <- read_xlsx("young48_dss_ferrozine_t35.xlsx")
  df <- df[,c(2,14:16)]
  colnames(df) <- df[2,]
  colnames(df)[1:2] <- c("id","iron_concentration")
  df <- df[-c(1,2,27),]
  df$gg_group <- paste(df$treatment, "+", df$diet, sep = "")
  df$gg_group <- factor(df$gg_group, levels = c("water+50","dss+50","water+500","dss+500"))
  df$diet <- factor(df$diet, levels = c("50","500"))
  df$gg_group <- df$diet
  df$iron_concentration <- as.numeric(df$iron_concentration)
  
  p = ironBoxplot(df, "iron_concentration", group = "diet", title = "Iron concentration in stools at day 35", y_axis_title = "yg of iron per g of stools", custom_colors = c("blue","red"),
                  stats = TRUE, test_results = c("***"), text_sizes = c(5))
  p1 <- p+
    scale_x_discrete(labels = c("50 ppm","500 ppm"))+
    labs(y = "µg Fe/g of stool", title = "Feces iron concentration\nat end of diet exposure (T35)")+
    guides(color = "none")
  p1
  
  # Stats 
  verifyStatsAssumptions(df, "diet" , "iron_concentration")
  wilcox.test(iron_concentration ~ diet, data = df)
  
  # Tf ferrozine assay for stools at last day
  sheets <- read_excel_allsheets("young48_dss_ferrozine_tF.xlsx")
  df <- as.data.frame(sheets["Ferrozine Stool Tfinal"])
  df <- df[4:53,c(3,14:16)]
  colnames(df) <- df[1,]
  colnames(df)[1:2] <- c("id","iron_concentration")
  df <- df[-c(1,26,48),]
  df$gg_group <- paste0(df$diet, ":", df$treatment)
  df$gg_group <- factor(df$gg_group, levels = c("50:water","500:water","50:dss","500:dss"))
  df$iron_concentration <- as.numeric(df$iron_concentration)
  
  p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Feces iron concentration\nat end timepoint (T112)", y_axis_title = "yg of iron/g of stools", custom_colors = c("blue","red","darkblue", "darkred"),
                  stats = TRUE, all.ns = TRUE)
  p2 <- p+
    scale_x_discrete(labels = c("50 ppm\nCtrl","500 ppm\nCtrl","50 ppm\nDSS","500 ppm\nDSS"))+
    labs(y = "µg Fe/g of stool")+
    guides(color = "none")
  p2
  
  # Stats 
  verifyStatsAssumptions(df, "gg_group" , "iron_concentration")
  pairwise.wilcox.test(df$iron_concentration, df$gg_group, p.adjust.method = "fdr")
  
  
  df$diet <- factor(df$diet, levels = c("50","500"))
  df$treatment <- factor(df$treatment, levels = c("water","dss"))
  pairwise_permuco(df, "iron_concentration", "diet", "treatment")
  
  
  # Tf ferrozine assay for liver
  df <- as.data.frame(sheets["Ferrozine Liver"])
  df <- df[1:51,c(3,6,8,14:18)]
  colnames(df) <- df[2,]
  colnames(df)[1:4] <- c("id","wet_weight","dry_weight","iron_concentration")
  df <- df[-c(1,2,27),]
  df$gg_group <- paste0(df$diet, ":", df$treatment)
  df$gg_group <- factor(df$gg_group, levels = c("50:water","500:water","50:dss","500:dss"))
  df$iron_concentration <- as.numeric(df$iron_concentration)
  df$liver_weight <- as.numeric(df$liver_weight)
  df$wet_weight <- as.numeric(df$wet_weight)
  df$dry_weight <- as.numeric(df$dry_weight)
  df$dry_to_wet_ratio <- df$dry_weight/df$wet_weight
  df <- df[-46,]
  
  
  # To calculate total iron in organ = iron concentration per g of dry weight*total organ weight
  # need to take into account the wet to dry ratio! (not so sure about that, need to check again)
  df$total_iron <- df$iron_concentration*df$liver_weight*df$dry_to_wet_ratio
  
  p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Liver iron concentration\nat end timepoint (T112)", y_axis_title = "yg of iron", custom_colors = c("blue","red","darkblue", "darkred"),
                  stats = TRUE, test_results = c("**","**","n.s.","n.s."), text_sizes = c(5,5,3,3), upper_margin = 30)
  p3 <- p+
    scale_x_discrete(labels = c("50 ppm\nCtrl","500 ppm\nCtrl","50 ppm\nDSS","500 ppm\nDSS"))+
    labs(y = "µg Fe/g dry weight")+
    guides(color = "none")
  p3
  
  # Stats 
  verifyStatsAssumptions(df, "gg_group" , "iron_concentration")
  pairwise.wilcox.test(df$iron_concentration, df$gg_group, p.adjust.method = "BH")
  
  df$diet <- factor(df$diet, levels = c("50","500"))
  df$treatment <- factor(df$treatment, levels = c("water","dss"))
  pairwise_permuco(df, "iron_concentration", "diet", "treatment")
  
  #Tf ferrozine assay for spleen
  df <- as.data.frame(sheets["Ferrozine Spleen"])
  df <- df[1:51,c(3,6,8,14:18)]
  colnames(df) <- df[2,]
  colnames(df)[1:4] <- c("id","wet_weight","dry_weight","iron_concentration")
  df <- df[-c(1,2,27),]
  df$gg_group <- paste0(df$diet, ":", df$treatment)
  df$gg_group <- factor(df$gg_group, levels = c("50:water","500:water","50:dss","500:dss"))
  df$iron_concentration <- as.numeric(df$iron_concentration)
  df$spleen_weight <- as.numeric(df$spleen_weight)
  df$wet_weight <- as.numeric(df$wet_weight)
  df$dry_weight <- as.numeric(df$dry_weight)
  df$dry_to_wet_ratio <- df$wet_weight/df$dry_weight
  df <- df[-46,]
  
  #To calculate total iron in organ = iron concentration per g of dry weight*total organ weight
  #need to take into account the wet to dry ratio!
  df$total_iron <- df$iron_concentration*df$spleen_weight*df$dry_to_wet_ratio
  
  p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Spleen iron concentration\nat end timepoint (T112)", y_axis_title = "yg of iron", custom_colors = c("blue","red","darkblue", "darkred"),
                  stats = TRUE, test_results = c("*","n.s.","*","*"), text_sizes = c(5,3,5,5), upper_margin = 400)
  p4 <- p+
    scale_x_discrete(labels = c("50 ppm\nCtrl","500 ppm\nCtrl","50 ppm\nDSS","500 ppm\nDSS"))+
    labs(y = "µg Fe/g dry weight")+
    guides(color = "none")
  p4

  # Stats 
  verifyStatsAssumptions(df, "gg_group" , "iron_concentration")
  anova <- aov(iron_concentration ~ diet * treatment, data = df) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  
  kruskal.test(iron_concentration ~ gg_group, data = df)
  dunnTest(iron_concentration ~ gg_group, data = df, method = "bonferroni")
  
  df$diet <- factor(df$diet, levels = c("50","500"))
  df$treatment <- factor(df$treatment, levels = c("water","dss"))
  pairwise_permuco(df, "iron_concentration", "diet", "treatment")
  
  # Combine plots
  fig1 <- p1 + p2 +p3 +p4 +
    plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 25))
    
  fig1
  
  existingDirCheck("~/CHUM_git/figures/memoire/dss/")
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig1.png", plot = fig1, width = 8, height = 9, dpi = 500)
}

# Figure 2 - dss induced colitis model
{
  # Set working directory
  setwd("experiments/finished exp/young-DSS-exp3")

  young_weight <- read.csv("young48_weight_cageChanges.csv", header = TRUE, sep = ";")
  young_weight <- young_weight[,-c(6:10)] # Remove weight measures before start of DSS
  young_weight <- weightDataManipulation(young_weight, groupInfoCols = 4, fromDay0 = TRUE)
  
  # Scatter plot with the four different treatments (diet combined with dss or control) - all timepoints
  young_weight_plot <- weightPlot(young_weight, percentage = TRUE, diet_only = FALSE, title = "Body weight")
  young_weight_plot <- young_weight_plot+
    scale_x_continuous(
      breaks = seq(min(young_weight$time_numeric), max(young_weight$time_numeric), by = 7), # breaks every 7 days
      labels = function(x) paste0("T", x)                    # "T" before each label
    )+
    labs(x = "Timepoint", color = "Group")+
    theme(axis.text.x = element_text(size = 7))
  
  # Mice DAI evolution 
  young_dss_followup <- read.csv("young48_dss_followup.csv", header = TRUE, sep=";")
  young_dss_followup <- dssFollowupManipulation(df = young_dss_followup,groupInfoCols = 4,dateStart = "2024-05-29",nbrDays = 5, negativeOnly = TRUE)
  young_dssflwup_plot <- dssDiseaseIndexPlot(young_dss_followup, statBarLast = TRUE, signifAnnotation = "*")
  young_dssflwup_plot <- young_dssflwup_plot+
    theme(axis.ticks.x = element_line(),
          legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5))
  
  # Statistics for last day of DAI
  verifyStatsAssumptions(df = young_dss_followup[young_dss_followup$time_numeric == "5",], group = "gg_group", measure = "index")
  wilcox.test(index ~ gg_group, data = young_dss_followup[young_dss_followup$time_numeric == "5",])
  
  # Mice dissection data
  dissec_young <- read.csv("young48_dss_dissection.csv", sep = ";", header = TRUE)
  dissec_young <- dissectionDataManipulation(dissec_young, groupInfoCols = 4, numerical = FALSE)
  dissec_young$gg_group <- factor(dissec_young$gg_group, levels = c("50 water","500 water","50 dss","500 dss"))
  dissec_young$colon_length_nrm <- dissec_young$colon_length/dissec_young$body_weight
  
  #boxplot for colon length (non std)
  young_colon <- dissecBoxplot(dissec_young,"colon", stats = TRUE, all.ns = TRUE)
  young_colon <- young_colon+
    labs(title = "Colon length at\nend timepoint (T112)")
  
  # Statistics
  verifyStatsAssumptions(df = dissec_young, group = "gg_group", measure = "colon_length")
  anova <- aov(colon_length ~ diet * treatment, data = dissec_young) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  
  # Combine plots
  fig2 <- young_weight_plot + young_dssflwup_plot + young_colon +
    plot_layout(widths = c(2, 1,1))+
    plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 25))
  
  existingDirCheck("~/CHUM_git/figures/memoire/dss/")
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig2.png", plot = fig2, width = 16, height = 4, dpi = 500)
}

# Additionnal figure - dss induced colitis model but in the adults
{
  setwd("experiments/finished exp/adults-all-exp/")
  
  adult_weight <- as.data.frame(read_xlsx("adult-all/body_weight2.xlsx"))
  # adult_weight <- adult_weight[!(adult_weight$batch == "C" & adult_weight$treatment == "abx"), ] # Remove failed abx mice from batch C
  adult_weight <- adult_weight[!adult_weight$batch == "D" & !adult_weight$treatment == "abx",] # Keep only DSS data and batch A B C
  adult_weight <- adult_weight[-12,] # Remove dead mouse/mice
  adult_weight <- weightDataManipulation(adult_weight, groupInfoCols = 5, fromDay0 = TRUE, arrangeDateColsFormat = FALSE)
  
  # Scatter plot with the four different treatments (diet combined with dss or control) - all timepoints
  adult_weight_plot <- weightPlot(adult_weight, percentage = TRUE, diet_only = FALSE, title = "Body weight")
  adult_weight_plot <- adult_weight_plot+
    scale_x_continuous(
      breaks = seq(min(adult_weight$time_numeric), max(adult_weight$time_numeric), by = 7), # breaks every 7 days
      labels = function(x) paste0("T", x)                    # "T" before each label
    )+
    labs(x = "Timepoint", color = "Group")+
    theme(axis.text.x = element_text(size = 7))
  adult_weight_plot
  
  # Mice DAI evolution
  adult_dss_followup <- read.csv("combined_adult_dss_followup.csv", header = TRUE, sep = ";")
  adult_dss_followup <- dssFollowupManipulation(df = adult_dss_followup,groupInfoCols = 4,dateStart = "2024-05-29",nbrDays = 5, negativeOnly = TRUE)
  adult_dssflwup_plot <- dssDiseaseIndexPlot(adult_dss_followup, statBarLast = TRUE, signifPositionShift = 0.5, signifSize = 5)
  adult_dssflwup_plot
  
  # Stats
  verifyStatsAssumptions(df = adult_dss_followup[adult_dss_followup$time_numeric == "5",], group = "gg_group", measure = "index")
  wilcox.test(index ~ gg_group, data = adult_dss_followup[adult_dss_followup$time_numeric == "5",])
  
  # Dissection data
  dissec_adult <- read_excel("dss-groups-dissection.xlsx")
  dissec_adult <- dissec_adult[dissec_adult$treatment != "abx",] # Remove abx mice
  dissec_adult <- dissec_adult[-12,] # Remove dead mouse
  dissec_adult[dissec_adult$ID == "36697B", 'liver weight'] <- 1 # For 36697B we are missing liver and spleen measures, replace by 1 and then exlcude if from analysis for spleen and liver measures
  dissec_adult[dissec_adult$ID == "36697B", 'spleen weight'] <- 1 
  colnames(dissec_adult)[6:9] <- str_replace(colnames(dissec_adult)[6:9], pattern = " ", replacement = "_") # Replace spaces by underscore for variables colnames
  dissec_adult <- dissectionDataManipulation(dissec_adult, groupInfoCols = 4)
  dissec_adult$gg_group <- factor(dissec_adult$gg_group, levels = c("50 water","500 water","50 dss","500 dss"))
  dissec_adult$colon_length_nrm <- dissec_adult$colon_length/dissec_adult$body_weight
  
  # Boxplot for colon length (non std)
  adult_dissec_cln <- dissecBoxplot(dissec_adult,"colon", all.ns = TRUE) 
  adult_dissec_cln <- adult_dissec_cln+
    labs(title = "Colon length at\nend timepoint (T112)")
  
  # Stats
  verifyStatsAssumptions(df = dissec_adult, group = "gg_group", measure = "colon_length")
  pairwise.wilcox.test(dissec_adult$colon_length, dissec_adult$gg_group, p.adjust.method = "fdr")
  
  # Combine plots
  fig2_adults <- adult_weight_plot + adult_dssflwup_plot + adult_dissec_cln +
    plot_layout(widths = c(2, 1,1))+
    plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 25))
  
  existingDirCheck("~/CHUM_git/figures/memoire/dss/")
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig2_adults.png", plot = fig2_adults, width = 16, height = 4, dpi = 500)
}

#### Microbiota related figures ####
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
  library(EnhancedVolcano)
  
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

sum(taxa_sums(ps)) # total number of reads
length(taxa_sums(ps)) # total number of ASVs

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

# Figure 3 - alpha diversity and beta diversity
{
  
  # Alpha diveristy timeline
  {
    existingDirCheck("../figures/Thibault_dss/alpha_diversity_timeline/")
    
    # Week as numeric
    sample_data(ps)$timepoint  <- as.numeric(as.character(sample_data(ps)$timepoint))
    sample_data(ps)$gg_group2 
    sample_data(ps)$gg_group2 <- factor(sample_data(ps)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
    
    #Estinate richness measures for dataset
    richness_data <- estimate_richness(ps, measures = c("Chao1", "Shannon", "InvSimpson", "Observed"))
    alpha_d <- cbind(as.data.frame(sample_data(ps)), richness_data)
    
    custom_colors <- c("blue","red","darkblue", "darkred")
    graphs <- alphaDiversityTimeline(ps, time = "timepoint", group = "gg_group2", custom_colors, semRibbons = TRUE)
    
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
    
    verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "49",], group = "gg_group2", measure = "observed")
    TukeyHSD(aov(observed ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "49",]))
    
    verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "54",], group = "gg_group2", measure = "observed")
    TukeyHSD(aov(observed ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "54",]))
    
    verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "final",], group = "gg_group2", measure = "observed")
    TukeyHSD(aov(observed ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "final",]))
    
    verifyStatsAssumptions(df = alpha_d[alpha_d$timepoint == "54",], group = "gg_group2", measure = "Observed")
    TukeyHSD(aov(Observed ~ gg_group2 , data = alpha_d[alpha_d$timepoint == "54",]))
    
    df <- data.frame(Comparison=c("50 vs 500 / Ctrl", "50 vs 500 / DSS","50 Ctrl vs 50 DSS","500 Ctrl vs 500 DSS"),
                     T0='',T35='',T49='',T54='',T112='')
    
    
    df <- df %>%
      mutate(across(2:4, ~ as.character(cut(.,
                                            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                            labels = c("***", "**", "*", "n.s.")))))
    
    ggtexttable(df, rows =NULL, theme = ttheme("light"))
    
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
    
    alpha_d_plot <- p1 / p2 /p3 +plot_layout(guides = "collect")
    alpha_d_plot
  }
  
  # Beta diversity
  {
    set.seed(100)
    
    phy_tree(ps_flt_diet) <- midpoint(tree) # Add rooted tree to phyloseq object
    
    # For diet timepoints
    sample_data(ps_flt_diet)$diet <- factor(sample_data(ps_flt_diet)$diet, labels = c("50 ppm","500 ppm"))
    
    # Weighted unifrac
    diet_beta_d_plot <- betaDiversityTimepoint2Factors(ps_flt_diet, sample_id = "sample_id", timeVariable = "timepoint",
                                   varToCompare =  "diet", distMethod ="wunifrac",
                                   transform = "rel_ab", customColors = c("blue","red"),
                                   font = "Arial", path = "../figures/Thibault_dss/beta_diversity/diet/", 
                                   additionnalAes = my_theme()+theme(plot.title = element_text(size = 16)), dim = c(4,12), displayPValue = TRUE,
                                   combineGraphs = TRUE, returnFig = TRUE, customTitles = 
                                     c("Weaning\n(T0)","End of diet\n(T35)","After 2 weeks of\niron sufficent diet (T49)"),
                                   positionPvalue = "right")
    
    # For diet + treatment at t54 and tfinal
    phy_tree(ps_flt_dss) <- midpoint(tree) # Add rooted tree to phyloseq object
    sample_data(ps_flt_dss)$gg_group2 <- factor(sample_data(ps_flt_dss)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl", "50 ppm DSS", "500 ppm DSS"))
    # Weighted unifrac => not very informative
    diet_dss_beta_d_plot <- betaDiversityTimepoint2Factors(ps_flt_dss, sample_id = "sample_id", timeVariable = "timepoint",
                                   varToCompare =  "gg_group2", distMethod ="wunifrac",
                                   customColors = c("blue","red","darkblue","darkred"),
                                   font = "Arial", path = "../figures/Thibault_dss/beta_diversity/diet_dss/",
                                   additionnalAes = my_theme()+theme(plot.title = element_text(size = 16)), dim = c(4,12), combineGraphs = TRUE, returnFig = TRUE,
                                   customTitles = c("End of DSS treatment\n(T54)","End timepoint\n(T112)"), pairwiseAdonis = FALSE, displayPValue = TRUE, hideLegend = TRUE)
    
    diet_dss_beta_d_plot 
    
    # Build table to display stats
    {
      df_1 <- read.xlsx("~/CHUM_git/figures/Thibault_dss/beta_diversity/diet_dss/week_49/wunifrac_week_49.xlsx")
      colnames(df_1)[2] <- "pvalue_t49"
      df_2 <- read.xlsx("~/CHUM_git/figures/Thibault_dss/beta_diversity/diet_dss/week_54/wunifrac_week_54.xlsx")
      colnames(df_2)[2] <- "pvalue_t54"
      df_3 <- read.xlsx("~/CHUM_git/figures/Thibault_dss/beta_diversity/diet_dss/week_112/wunifrac_week_112.xlsx")
      colnames(df_3)[2] <- "pvalue_t112"
      df <- cbind(df_1, df_2[2], df_3[2])
      
      df <- df[c(1,4,2,3),]
      colnames(df) <- c("Comparison","T49","T54","T112")
      df$Comparison <- c("50 vs 500 / Ctrl", "50 vs 500 / DSS","50 Ctrl vs 50 DSS","500 Ctrl vs 500 DSS")
      df <- df %>%
        mutate(across(2:4, ~ as.character(cut(.,
                                              breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                              labels = c("***", "**", "*", "n.s.")))))
      
      ggtexttable(df, rows =NULL, theme = ttheme("light"))
    }
    
    # Weighted unifrac for dss groups - last timepoint
    ps_sub <- prune_samples(sample_data(ps_tfinal_flt)$treatment == "dss", ps_tfinal_flt)
    sample_data(ps_sub)$diet <- factor(sample_data(ps_sub)$diet, labels = c("50 ppm DSS","500 ppm DSS"))
    phy_tree(ps_sub) <- midpoint(tree) # Add rooted tree to phyloseq object
    DSS_50v500_beta_d_plot <- betaDiversityTimepoint2Factors(ps_sub , sample_id = "sample_id", timeVariable = "timepoint",
                                   varToCompare =  "diet", distMethod ="wunifrac",
                                   customColors = c("darkblue","darkred"),
                                   font = "Arial", path = "../figures/Thibault_dss/beta_diversity/dss_only/",
                                   additionnalAes = my_theme()+theme(plot.title = element_text(face = "italic", size = 12)), dim = c(2.5,4), returnFig = TRUE, displayPValue = TRUE,
                                   customTitles = c("50 ppm DSS vs 500 ppm DSS"), combineGraphs = FALSE, hideLegend = TRUE)
    
    # Weighted unifrac for 50 ctrl vs 50 dss - last timepoint
    ps_sub <- prune_samples(sample_data(ps_tfinal_flt)$diet == "50", ps_tfinal_flt)
    sample_data(ps_sub)$treatment <- factor(sample_data(ps_sub)$treatment, labels = c("50 ppm Ctrl","50 ppm DSS"))
    phy_tree(ps_sub) <- midpoint(tree) # Add rooted tree to phyloseq object
    Ctrl50vDSS50_beta_d_plot <- betaDiversityTimepoint2Factors(ps_sub , sample_id = "sample_id", timeVariable = "timepoint",
                                                             varToCompare =  "treatment", distMethod ="wunifrac",
                                                             customColors = c("blue","darkblue"),
                                                             font = "Arial", path = "../figures/Thibault_dss/beta_diversity/dss_only/",
                                                             additionnalAes = my_theme()+theme(plot.title = element_text(face = "italic", size = 12)), dim = c(2.5,4), returnFig = TRUE, displayPValue = TRUE,
                                                             customTitles = c("50 ppm Ctrl vs 50 ppm DSS"), combineGraphs = FALSE, hideLegend = TRUE)
    
    
  }
  
  # # Combine all graphs
  # alpha_d_plot <- p1 / p2 /p3 +plot_layout(guides = "collect")
  # 
  # alpha_d_plot
  # diet_beta_d_plot
  # diet_dss_beta_d_plot
  # DSS_50v500_beta_d_plot
  # Ctrl50vDSS50_beta_d_plot
  
  dss_stack <- DSS_50v500_beta_d_plot / Ctrl50vDSS50_beta_d_plot
  
  # Bottom row: diet_dss_plot next to the stack
  bottom_row <- (diet_dss_beta_d_plot | dss_stack) +
    plot_layout(widths = c(2.4, 1))
  
  # Left column: top = diet_beta_d_plot, bottom = bottom_row, plus padding between them
  left_block <- diet_beta_d_plot /
    bottom_row +
    plot_layout(heights = c(1, 1))  # second row gets 1.2× the height of the first
  
  final_plot <- left_block | alpha_d_plot
  
  # Shrink alpha diversity column if desired
  fig3 <- final_plot +
    plot_layout(widths = c(3, 1))
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig3.png", plot = fig3, width = 18, height = 9, dpi = 500)
  
  
}
