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
  library(legendry)
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
  df$treatment <- factor(df$treatment, levels = c("water","dss"), labels = c("Water","DSS"))
  df$gg_group <- df$diet
  df$iron_concentration <- as.numeric(df$iron_concentration)
  
  p = ironBoxplot(df, "iron_concentration", group = "diet", shape = "treatment", title = "Iron concentration in stools at day 35", y_axis_title = "yg of iron per g of stools", custom_colors = c("#95BECF","#F2AA84"),
                  stats = TRUE, test_results = c("***"), text_sizes = c(5))
  p1 <- p+
    scale_x_discrete(labels = c("50 ppm","500 ppm"))+
    labs(y = "µg Fe/g dry weight", title = "Stool iron\nat end of exposure (T35)")+
    guides(color = "none")+
    theme(plot.title = element_text(size = 10),
          axis.title.y = element_text(size = 10))
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
  df$gg_group <- interaction(df$diet, df$treatment, sep=":")
  df$gg_group <- factor(df$gg_group, levels = c("50:water","500:water","50:dss","500:dss"), labels = c(c("50:Water","500:Water","50:DSS","500:DSS")))
  df$iron_concentration <- as.numeric(df$iron_concentration)
  df$treatment <- factor(df$treatment, levels = c("water","dss"), labels = c("Water","DSS"))
  
  p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Stool iron concentration\nat end of experiment", y_axis_title = "yg of iron/g dry weight", custom_colors = c("#95BECF","#F2AA84","#325B6C","#B22222"),
                  stats = TRUE, test_results = c("n.s.","n.s.","n.s.","n.s."), upper_margin = 200, text_sizes = c(3,3,3,3), shape = "treatment")
  p2 <- p +
    guides(x = "axis_nested")+
    labs(y = "µg Fe/g dry weight")+
    theme(plot.title = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  p2
  
  # Stats 
  verifyStatsAssumptions(df, "gg_group" , "iron_concentration")
  pairwise.wilcox.test(df$iron_concentration, df$gg_group, p.adjust.method = "fdr")
  
  
  df$diet <- factor(df$diet, levels = c("50","500"))
  df$treatment <- factor(df$treatment, levels = c("water","dss"))
  pairwise_permuco(df, "iron_concentration", "diet", "treatment")
  
  # Combined stool iron graph
  {
    #T35 ferrozine assay for stools
    df <- read_xlsx("young48_dss_ferrozine_t35.xlsx")
    df <- df[,c(2,14:16)]
    colnames(df) <- df[2,]
    colnames(df)[1:2] <- c("id","iron_concentration")
    df <- df[-c(1,2,27),]
    df$gg_group <- paste(df$treatment, "+", df$diet, sep = "")
    df$gg_group <- factor(df$gg_group, levels = c("water+50","dss+50","water+500","dss+500"))
    df$diet <- factor(df$diet, levels = c("50","500"), labels = c("50 ppm","500 ppm"))
    df$treatment <- factor(df$treatment, levels = c("water","dss"), labels = c("Water","DSS"))
    df$gg_group <- df$diet
    df$iron_concentration <- as.numeric(df$iron_concentration)
    dfT35 <- df
    
    # Tf ferrozine assay for stools at last day
    sheets <- read_excel_allsheets("young48_dss_ferrozine_tF.xlsx")
    df <- as.data.frame(sheets["Ferrozine Stool Tfinal"])
    df <- df[4:53,c(3,14:16)]
    colnames(df) <- df[1,]
    colnames(df)[1:2] <- c("id","iron_concentration")
    df <- df[-c(1,26,48),]
    df$gg_group <- interaction(df$diet, df$treatment, sep=":")
    df$gg_group <- factor(df$gg_group, levels = c("50:water","500:water","50:dss","500:dss"), labels = c(c("50:Water","500:Water","50:DSS","500:DSS")))
    df$iron_concentration <- as.numeric(df$iron_concentration)
    df$treatment <- factor(df$treatment, levels = c("water","dss"), labels = c("Water","DSS"))
    dfTfinal <- df
    
    df_merged <- rbind(dfT35, dfTfinal)
    df_merged$gg_group 
    
    p = ironBoxplot(df_merged, "iron_concentration", group = "gg_group", title = "Stool iron", custom_colors = c("#95BECF","#F2AA84","#95BECF","#F2AA84","#325B6C","#B22222"),
                    stats = FALSE, shape = "treatment")
    p_combined <- p +
      guides(x = "axis_nested")+
      labs(y = "µg Fe/g dry weight")+
      theme(plot.title = element_text(size = 10),
            axis.title.y = element_text(size = 10))+
      geom_vline(xintercept = 2.5,
                 color      = "black",
                 linetype   = "dashed",
                 size       = 0.4)
    p_combined
    
    # Get factor levels in order
    groups <- levels(df_merged$gg_group)
    
    comparisons_list <- list(
      c(groups[1], groups[2]),  # "50 ppm" vs "500 ppm"
      c(groups[3], groups[4]),  # "50 water  vs "500 water"
      c(groups[5], groups[6]),  # "50 dss" vs "500 dss"
      c(groups[3], groups[5]), # "50 water" vs "50 dss"
      c(groups[4], groups[6])) # "500 water" vs "500 dss" 
  
  y_positions <- c(
    max(df_merged$iron_concentration),
    min(df_merged$iron_concentration)+ 3*min(df_merged$iron_concentration),
    min(df_merged$iron_concentration)+ 6*min(df_merged$iron_concentration),
    min(df_merged$iron_concentration)+ 9*min(df_merged$iron_concentration),
    min(df_merged$iron_concentration)+ 12*min(df_merged$iron_concentration))
  
  test_results <- c("***","n.s","n.s","n.s","n.s")
  text_sizes <- c(5,3,3,3,3)
  vjustList <- c(0,0,0,0,0)
  
  # Apply each geom_signif layer individually
  for (i in seq_along(comparisons_list)) {
    p_combined <- p_combined +
      geom_signif(
        comparisons = list(comparisons_list[[i]]),
        annotations = test_results[i],
        y_position = y_positions[i],
        tip_length = 0.03,
        color = "black",
        size = 0.4,
        textsize = text_sizes[i],
        margin_top = 0.1,
        vjust = vjustList[i])}
    
  }
  p_combined
  
  # Tf ferrozine assay for liver
  df <- as.data.frame(sheets["Ferrozine Liver"])
  df <- df[1:51,c(3,6,8,14:18)]
  colnames(df) <- df[2,]
  colnames(df)[1:4] <- c("id","wet_weight","dry_weight","iron_concentration")
  df <- df[-c(1,2,27),]
  df$gg_group <- interaction(df$diet, df$treatment, sep=":")
  df$gg_group <- factor(df$gg_group, levels = c("50:water","500:water","50:dss","500:dss"), labels = c(c("50:Water","500:Water","50:DSS","500:DSS")))
  df$iron_concentration <- as.numeric(df$iron_concentration)
  df$liver_weight <- as.numeric(df$liver_weight)
  df$wet_weight <- as.numeric(df$wet_weight)
  df$dry_weight <- as.numeric(df$dry_weight)
  df$dry_to_wet_ratio <- df$dry_weight/df$wet_weight
  df <- df[-46,]
  df$treatment <- factor(df$treatment, levels = c("water","dss"), labels = c("Water","DSS"))
  
  # To calculate total iron in organ = iron concentration per g of dry weight*total organ weight
  # need to take into account the wet to dry ratio! (not so sure about that, need to check again)
  df$total_iron <- df$iron_concentration*df$liver_weight*df$dry_to_wet_ratio
  
  p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Liver iron", y_axis_title = "yg of iron", custom_colors = c("#95BECF","#F2AA84","#325B6C","#B22222"),
                  stats = TRUE, test_results = c("**","**","n.s.","n.s."), text_sizes = c(5,5,3,3), upper_margin = 30, vjustList = c(0.5,0.5,0,0), shape = "treatment")
  p3 <- p+
    guides(x = "axis_nested")+
    labs(y = "µg Fe/g dry weight")+
    theme(plot.title = element_text(size = 10),
          axis.title.y = element_text(size = 10))
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
  df$gg_group <- interaction(df$diet, df$treatment, sep=":")
  df$gg_group <- factor(df$gg_group, levels = c("50:water","500:water","50:dss","500:dss"), labels = c(c("50:Water","500:Water","50:DSS","500:DSS")))
  df$iron_concentration <- as.numeric(df$iron_concentration)
  df$spleen_weight <- as.numeric(df$spleen_weight)
  df$wet_weight <- as.numeric(df$wet_weight)
  df$dry_weight <- as.numeric(df$dry_weight)
  df$dry_to_wet_ratio <- df$wet_weight/df$dry_weight
  df <- df[-46,]
  df$treatment <- factor(df$treatment, levels = c("water","dss"), labels = c("Water","DSS"))
  
  
  #To calculate total iron in organ = iron concentration per g of dry weight*total organ weight
  #need to take into account the wet to dry ratio!
  df$total_iron <- df$iron_concentration*df$spleen_weight*df$dry_to_wet_ratio
  
  p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Spleen iron", y_axis_title = "yg of iron", custom_colors = c("#95BECF","#F2AA84","#325B6C","#B22222"),
                  stats = TRUE, test_results = c("*","n.s.","*","*"), text_sizes = c(5,3,5,5), upper_margin = 400, vjustList = c(0.5,0,0.5,0.5), shape = "treatment")
  p4 <- p+
    guides(x = "axis_nested")+
    labs(y = "µg Fe/g dry weight")+
    theme(plot.title = element_text(size = 10),
          axis.title.y = element_text(size = 10))
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
  fig1 <- wrap_elements(p1 + p2 +p3 +p4) +
    plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 20))
    
  fig1
  
  existingDirCheck("~/CHUM_git/figures/memoire/dss/")
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig1.png", plot = fig1, width = 5, height = 6, dpi = 500)
  
  # Another version with grouped stool iron
  fig1_v2 <- p_combined/(p3+p4)+
    plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 20))
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig1v2.png", plot = fig1_v2, width = 4, height = 5, dpi = 500)
}

# Figure 2 - dss induced colitis model
{
  # Set working directory
  setwd("experiments/finished exp/young-DSS-exp3")

  young_weight <- read.csv("young48_weight_cageChanges.csv", header = TRUE, sep = ";")
  young_weight <- young_weight[,-c(6:10)] # Remove weight measures before start of DSS
  young_weight <- weightDataManipulation(young_weight, groupInfoCols = 4, fromDay0 = TRUE)
  
  # Scatter plot with the four different treatments (diet combined with dss or control) - all timepoints
  young_weight_plot <- weightPlot(young_weight, percentage = FALSE, diet_only = FALSE, title = "Body weight")
  young_weight_plot <- young_weight_plot+
    scale_x_continuous(
      breaks = seq(min(young_weight$time_numeric), max(young_weight$time_numeric), by = 7) # breaks every 7 days
    )+
    labs(x = "Days", color = NULL)+
    theme(axis.text.x = element_text(size = 7))
  young_weight_plot
  
  # Mice DAI evolution 
  young_dss_followup <- read.csv("young48_dss_followup.csv", header = TRUE, sep=";")
  young_dss_followup <- dssFollowupManipulation(df = young_dss_followup,groupInfoCols = 4,dateStart = "2024-05-29",nbrDays = 5, negativeOnly = TRUE)
  young_dssflwup_plot <- dssDiseaseIndexPlot(young_dss_followup, statBarLast = TRUE, signifAnnotation = "*", signifPositionShift = 0.2)
  young_dssflwup_plot <- young_dssflwup_plot+
    theme(axis.ticks.x = element_line(),
          legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5))
  young_dssflwup_plot
  
  # Statistics for last day of DAI
  verifyStatsAssumptions(df = young_dss_followup[young_dss_followup$time_numeric == "5",], group = "gg_group", measure = "index")
  wilcox.test(index ~ gg_group, data = young_dss_followup[young_dss_followup$time_numeric == "5",])
  
  # Mice dissection data
  dissec_young <- read.csv("young48_dss_dissection.csv", sep = ";", header = TRUE)
  dissec_young <- dissectionDataManipulation(dissec_young, groupInfoCols = 4, numerical = FALSE)
  dissec_young$gg_group <- interaction(dissec_young$diet, dissec_young$treatment, sep = ":")
  dissec_young$gg_group <- factor(dissec_young$gg_group, levels = c("50:water","500:water","50:dss","500:dss"), labels = c("50:Water","500:Water","50:DSS","500:DSS"))
  dissec_young$colon_length_nrm <- dissec_young$colon_length/dissec_young$body_weight
  dissec_young$treatment <- factor(dissec_young$treatment, levels = c("water","dss"))
  
  #boxplot for colon length (non std)
  young_colon <- dissecBoxplot(dissec_young,"colon", stats = TRUE, all.ns = FALSE, test_results = c("n.s.","n.s.","n.s.","n.s."), text_sizes = c(3,3,3,3), upper_margin = 0.5)
  young_colon <- young_colon+
    labs(title = "Colon length")+
    guides(x = "axis_nested")
  young_colon
  
  # Statistics
  verifyStatsAssumptions(df = dissec_young, group = "gg_group", measure = "colon_length")
  anova <- aov(colon_length ~ diet * treatment, data = dissec_young) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  
  # Combine plots
  fig2 <- young_weight_plot + young_dssflwup_plot + young_colon +
    plot_layout(widths = c(1.75, 1,1))+
    plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 20))
  
  existingDirCheck("~/CHUM_git/figures/memoire/dss/")
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig2.png", plot = fig2, width =12, height = 3, dpi = 500)
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
  adult_weight_plot <- weightPlot(adult_weight, percentage = FALSE, diet_only = FALSE, title = "Body weight")
  adult_weight_plot <- adult_weight_plot+
    scale_x_continuous(
      breaks = seq(min(adult_weight$time_numeric), max(adult_weight$time_numeric), by = 7) # breaks every 7 days                   # "T" before each label
    )+
    labs(x = "Days", color = "Group")+
    theme(axis.text.x = element_text(size = 7))+
    ylim(10,NA)
  adult_weight_plot
  
  # Mice DAI evolution
  adult_dss_followup <- read.csv("combined_adult_dss_followup.csv", header = TRUE, sep = ";")
  adult_dss_followup <- dssFollowupManipulation(df = adult_dss_followup,groupInfoCols = 4,dateStart = "2024-05-29",nbrDays = 5, negativeOnly = TRUE)
  adult_dssflwup_plot <- dssDiseaseIndexPlot(adult_dss_followup, statBarLast = TRUE, signifPositionShift = 0.6, signifSize = 4)
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
  dissec_adult$gg_group <- interaction(dissec_adult$diet, dissec_adult$treatment, sep = ":")
  dissec_adult$gg_group <- factor(dissec_adult$gg_group, levels = c("50:water","500:water","50:dss","500:dss"), labels = c("50:Water","500:Water","50:DSS","500:DSS"))
  dissec_adult$colon_length_nrm <- dissec_adult$colon_length/dissec_adult$body_weight
  
  # Boxplot for colon length (non std)
  adult_dissec_cln <- dissecBoxplot(dissec_adult,"colon", all.ns = FALSE, stats = TRUE, test_results = c("n.s.","n.s.","n.s.","n.s."), text_sizes = c(3,3,3,3), upper_margin = 0.5) 
  adult_dissec_cln <- adult_dissec_cln+
    labs(title = "Colon length")+
    guides(x = "axis_nested")
  adult_dissec_cln
  
  # Stats
  verifyStatsAssumptions(df = dissec_adult, group = "gg_group", measure = "colon_length")
  pairwise.wilcox.test(dissec_adult$colon_length, dissec_adult$gg_group, p.adjust.method = "fdr")
  
  # Combine plots
  fig2_adults <- adult_weight_plot + adult_dssflwup_plot + adult_dissec_cln +
    plot_layout(widths = c(1.75, 1,1))+
    plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 20))
  
  existingDirCheck("~/CHUM_git/figures/memoire/dss/")
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig2_adults.png", plot = fig2_adults, width = 12, height = 3, dpi = 500)
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
    
    custom_colors <- c("#95BECF","#F2AA84","#325B6C","#B22222")
    graphs <- alphaDiversityTimeline(ps, time = "timepoint", group = "gg_group2", custom_colors, semRibbons = TRUE)
    
    # Observed
    p1 <- graphs[[4]]+
      scale_x_continuous(
        breaks = seq(min(sample_data(ps)$timepoint), max(sample_data(ps)$timepoint), by = 7) # breaks every 7 days                    # "T" before each label
      ) +
      labs(y = "Number of species", x = "Days", color = "Group", fill = "", title = 'Species observed')+
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
        breaks = seq(min(sample_data(ps)$timepoint), max(sample_data(ps)$timepoint), by = 7) # breaks every 7 days                   # "T" before each label
      ) +
      labs(y = "Shannon index", x = "Days", color = "Group", fill = "", title = 'Shannon')+
      guides(fill = "none")+
      theme(axis.ticks.x = element_line(),
            legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5),
            axis.text.x = element_text(size = 6))+
      ylim(0,NA)
    
    # Inverse Simpson
    p3 <- graphs[[3]]+
      scale_x_continuous(
        breaks = seq(min(sample_data(ps)$timepoint), max(sample_data(ps)$timepoint), by = 7) # breaks every 7 days                  # "T" before each label
      ) +
      labs(y = "Inverse Simpson index", x = "Days", color = "Group", fill = "", title = 'Inverse Simpson')+
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
                                   transform = "rel_ab", customColors = c("#95BECF","#F2AA84"),
                                   font = "Arial", path = "../figures/Thibault_dss/beta_diversity/diet/", 
                                   additionnalAes = my_theme()+theme(plot.title = element_text(size = 10)), dim = c(4,12), displayPValue = TRUE,
                                   combineGraphs = TRUE, returnFig = TRUE, customTitles = 
                                     c("Weaning","End of diets","2 weeks of\niron sufficent diet"),
                                   positionPvalue = "right")
    
    # For diet + treatment at t54 and tfinal
    phy_tree(ps_flt_dss) <- midpoint(tree) # Add rooted tree to phyloseq object
    sample_data(ps_flt_dss)$gg_group2 <- factor(sample_data(ps_flt_dss)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl", "50 ppm DSS", "500 ppm DSS"))
    # Weighted unifrac => not very informative
    diet_dss_beta_d_plot <- betaDiversityTimepoint2Factors(ps_flt_dss, sample_id = "sample_id", timeVariable = "timepoint",
                                   varToCompare =  "gg_group2", distMethod ="wunifrac",
                                   customColors = c("#95BECF","#F2AA84","#325B6C","#B22222"),
                                   font = "Arial", path = "../figures/Thibault_dss/beta_diversity/diet_dss/",
                                   additionnalAes = my_theme()+theme(plot.title = element_text(size = 10)), dim = c(4,12), combineGraphs = TRUE, returnFig = TRUE,
                                   customTitles = c("End of DSS treatment","End of experiment"), pairwiseAdonis = FALSE, displayPValue = TRUE, hideLegend = TRUE)
    
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
                                   customColors = c("#325B6C","#B22222"),
                                   font = "Arial", path = "../figures/Thibault_dss/beta_diversity/dss_only/",
                                   additionnalAes = my_theme()+theme(plot.title = element_text(face = "bold", size = 8)), dim = c(2.5,4), returnFig = TRUE, displayPValue = TRUE,
                                   customTitles = c("50 ppm DSS vs 500 ppm DSS"), combineGraphs = FALSE, hideLegend = TRUE)
    
    # Weighted unifrac for 50 ctrl vs 50 dss - last timepoint
    ps_sub <- prune_samples(sample_data(ps_tfinal_flt)$diet == "50", ps_tfinal_flt)
    sample_data(ps_sub)$treatment <- factor(sample_data(ps_sub)$treatment, labels = c("50 ppm Ctrl","50 ppm DSS"))
    phy_tree(ps_sub) <- midpoint(tree) # Add rooted tree to phyloseq object
    Ctrl50vDSS50_beta_d_plot <- betaDiversityTimepoint2Factors(ps_sub , sample_id = "sample_id", timeVariable = "timepoint",
                                                             varToCompare =  "treatment", distMethod ="wunifrac",
                                                             customColors = c("#95BECF","#325B6C"),
                                                             font = "Arial", path = "../figures/Thibault_dss/beta_diversity/dss_only/",
                                                             additionnalAes = my_theme()+theme(plot.title = element_text(face = "bold", size = 8)), dim = c(2.5,4), returnFig = TRUE, displayPValue = TRUE,
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
  
  final_plot <- left_block
  
  # Shrink alpha diversity column if desired
  fig3beta <- final_plot +
    plot_layout(widths = c(3, 1))
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig3beta.png", plot = fig3beta, width = 9, height = 6, dpi = 500)
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig3alpha.png", plot = alpha_d_plot, width = 5, height = 6, dpi = 500)
  
}

# Figure 4 - bacteria relative abundance and overall composition of the gut microbiota community
{
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
    palette_vec <- c("#483C46","#70AE6E","#FFC300","#900C3F","#C70039","#3C6E71","#BEEE62","#FF5733","black")
    palette_vec <- c("#c7522a","#e5c185","#f0daa5","#fbf2c4","#b8cdab","#74a892","#008585","#004343")
    palette_vec <- c("#003f5c","#58508d","#8a508f","#bc5090","#de5a79","#ff6361","#ff8531","#ffa600")
    palette_vec <- c()
    palette_vec <- c("#715660","#846470","#566071","#727f96","#c39f72","#d5bf95","#667762","#879e82")
    
    phylum_pal  <- setNames(palette_vec, unique_phyla)
    # Shuffle (randomize) the values, keeping the same names
    phylum_pal[] <- sample(phylum_pal)
    tax_data$data$taxa$color <- phylum_pal[prop_phylum]
    
    p_tree <- heat_tree(
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
    
    # Data frame for legend
    legend_df <- data.frame(
      Phylum = names(phylum_pal),
      Color  = phylum_pal
    )
    
    library(forcats)  # For fct_rev
    
    legend_df$Phylum <- factor(legend_df$Phylum, levels = rev(legend_df$Phylum)) # For top-down order
    
    p_legend <- ggplot(legend_df) +
      geom_point(aes(x = 0, y = Phylum, color = Phylum), size = 3) +
      geom_text(aes(x = 0, y = Phylum, label = Phylum), hjust = -0.1, size = 3) + # Play with x/hjust for spacing
      scale_color_manual(values = phylum_pal) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)
      ) +
      labs(title = NULL)
    
    cauliflower_tree <- p_tree +
      inset_element(
        p_legend +
          theme(plot.background = element_rect(fill = NA, color = NA)),
        left = -0.4, bottom = 0.2, right = 0.2, top = 0.8, align_to = "full"
      )
    cauliflower_tree
  }
  
  # Stackbar at last day of dss
  {
    # First, timepoints and groups must be ordered properly and as factors
    sample_data(ps_t54_flt)$diet 
    sample_data(ps_t54_flt)$treatment 
    sample_data(ps_t54_flt)$gg_group2 <- factor(sample_data(ps_t54_flt)$gg_group2, levels = c("50:water", "50:dss", "500:water","500:dss"))
    sample_data(ps_t54_flt)$timepoint
    
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
      legend_size = 5,
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
    
    facet_scales <- list(
      scale_x_discrete(labels = as.character(1:12)),
      scale_x_discrete(labels = as.character(25:35)),
      scale_x_discrete(labels = as.character(13:24)),
      scale_x_discrete(labels = as.character(36:45))
    )
    
    # Custom the plot
    stackbar_t54 <- diet_dss_phyla_fam$plot + 
      facet_wrap2(~ gg_group2, 
                  scales  = "free_x", nrow = 2, ncol = 2,
                  strip = strip_themed(background_x = elem_list_rect(fill = c("#95BECF","#325B6C","#F2AA84","#B22222"))), 
                  labeller = as_labeller(c("50:water" = "50 ppm Ctrl",
                                           "50:dss" = "50 ppm DSS",
                                           "500:water" = "500 ppm Ctrl",
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
    
    # Saving the plot and the associated stats
    existingDirCheck("../figures/Thibault_dss/stackbar")
    writeStackbarExtendedSigTable(main_table = diet_dss_phyla_fam$significant_table_main, includeSubTable = TRUE, sub_table = diet_dss_phyla_fam$significant_table_sub, filepath = "../figures/Thibault_dss/stackbar/diet_dss_stackbar_stats.xlsx")
    
    # pvalues heatmap for the main lvl stats
    pvalHmapPhyla <- pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_dss/stackbar/diet_dss_stackbar_stats.xlsx")),
                selected_comparisons = c("50:water_vs_50:dss", "500:water_vs_500:dss","50:water_vs_500:water","50:dss_vs_500:dss"), displayChangeArrows = TRUE, displayPValues = FALSE,
                txn_lvl="Phylum", lvl = "main", taxons = diet_dss_phyla_fam$main_names, group = "gg_group2", path)
    pvalHmapPhyla <- pvalHmapPhyla+
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank())+
      guides(fill = "none")
    
    # pvalues heatmap for the sub lvl stats
    pvalHmapFamily <- pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_dss/stackbar/diet_dss_stackbar_stats.xlsx")),
                    selected_comparisons = c("50:water_vs_50:dss", "500:water_vs_500:dss","50:water_vs_500:water","50:dss_vs_500:dss"),
                    txn_lvl="Family", lvl = "sub", taxons =  diet_dss_phyla_fam$sub_names, group = "gg_group2", displayPValues = FALSE, displayChangeArrows = TRUE, path) # You can add [!grepl("Others", x = iron_exp_family$sub_names)] to remove "others"
    pvalHmapFamily <- pvalHmapFamily+scale_x_discrete(labels = c("50 Ctrl vs 50 DSS", "500 Ctrl vs 500 DSS", "50 Ctrl vs 500 Ctrl", "50 DSS vs 500 DSS"))+
      theme(text = element_text(family = "Arial"),
            axis.text.y = element_blank(),
            axis.text.x = element_blank())
    
    # Add stats hmap to stackbar figure
    full_stackbar <- (stackbar_t54 | wrap_elements((pvalHmapPhyla / pvalHmapFamily)))+
      plot_layout(widths = c(1, 0.5), ncol = 2)
    full_stackbar
    
    ggsave(filename = "~/CHUM_git/figures/memoire/dss/stackbar.png", width = 8, height = 6, dpi = 500)
    
    
    # Final stackbar with stats heatmap
    # Extract stackbar legend
    stackbar_lgd <- get_legend(stackbar_t54)
    
    # Main plot without legend:
    stackbar_t54 <- stackbar_t54 + theme(legend.position = "none")
    
    hmap <- (pvalHmapPhyla / pvalHmapFamily)
    full_stackbar <- ((stackbar_t54 | wrap_elements(full = hmap) | wrap_elements(full = stackbar_lgd)))+
      plot_layout(widths = c(2, 1, 1), ncol = 3)
    
    
    # Combine them side-by-side
    (stackbar_t54 | (pvalHmapPhyla / pvalHmapFamily) | wrap_elements(stackbar_lgd))+
      plot_layout(widths = c(4, 1, 1), ncol = 3)
    stackbar_t54_full <- a
    
    
    (stackbar_t54 + plot_spacer() + (pvalHmapPhyla / pvalHmapFamily))+
      plot_layout(widths = c(4, -12, 4), guides = "collect")
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
    
    sample_data(ps_flt_all)$gg_group2 <- factor(sample_data(ps_flt_all)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
    
    p <- plot_timeline_2_groups(
      ps_object = ps_flt_all,
      exp_group =  "gg_group2", # must be as factor
      time_group = "timepoint", # must be as factor
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
    
    chronobiome <- p+
      facet_wrap2(~ gg_group2, 
                  scales  = "free_x", nrow = 2, ncol = 2,
                  strip = strip_themed(background_x = elem_list_rect(fill = c("#95BECF","#F2AA84","#325B6C","#B22222"))))+
      scale_x_continuous(breaks = seq(min(as.numeric(levels(sample_data(ps_flt_all)$timepoint))), max(as.numeric(levels(sample_data(ps_flt_all)$timepoint))), by = 7))+
      theme(axis.text.x = element_text(size = 6, face = "bold"))+
      labs(x = "Days")+
        guides(fill = "none", alpha = "none")
  }
  
  # Final figure 4
  fig4 <- ((cauliflower_tree | stackbar_t54) / chronobiome)+
    plot_layout()
  
  fig4 <- ((cauliflower_tree | stackbar_t54) / chronobiome) +
    plot_layout(
      heights = c(1, 1.5),
      widths  = c(1, 1),
      guides  = "collect"
    )+
    plot_annotation(tag_levels = list(c("A", "","B","C"))) &
    theme(panel.spacing = unit(0, "pt"),
          plot.tag = element_text(face = "bold", size = 25))
  
  ((cauliflower_tree | full_stackbar))+
    plot_layout(widths = c(1,2))
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig4_top.png", width = 16, height = 9, dpi = 500)
  
  
  (wrap_elements(full = top) / chronobiome)+
    plot_layout(heights = c(2,1))
  
  
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig4.png", plot = fig4, width = 11, height = 9, dpi = 500)
  
  
}

# Figure 5 - bacterial species relative abundance
{
  # Volcano plots - 50 VS 500 DSS last timepoint
  {
      deseq <- phyloseq_to_deseq2(ps_tfinal_flt, ~ treatment+diet+diet:treatment)
      deseq <- DESeq(deseq, test="Wald", fitType = "parametric")
      print(resultsNames(deseq))
      
      # Manually add species names that we identified with BLASTn
      tax_table(ps_tfinal_flt)["ASV8","Species"] <- "murinus (ASV8)"
      tax_table(ps_tfinal_flt)["ASV38","Species"] <- "murinus (ASV38)"
      tax_table(ps_tfinal_flt)["ASV19","Species"] <- "intestinale"
      tax_table(ps_tfinal_flt)["ASV6","Species"] <- "rodentium (ASV6)"
      tax_table(ps_tfinal_flt)["ASV1","Species"] <- "rodentium (ASV1)"

      #For a given taxononical levels, creates graph for each timepoint, displaying a volcano plot for taxononical level of interest
      volcanoPlot50vs500dss <- volcanoPlot2GroupsMultifactorDesign(ps = ps_tfinal_flt, deseq = deseq, varToCompare = "diet",
                                          taxa = "Species", threshold = 0.05, FCcutoff = 0.49, customColors = NULL,
                                          FDR = TRUE, includeUnknownSpecies = TRUE, selectedComparison = 2,
                                          title = "50 ppm DSS VS 500 ppm DSS\nat end timepoint (T112)")+
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

            axis.ticks.x = element_line()
          )+ 
        guides(
            color = guide_legend(ncol = 1, byrow = TRUE),
            fill  = guide_legend(ncol = 1, byrow = TRUE)  # if applicable
          )
      
      volcanoPlot50vs500dss
    
  }
  
  # Individual graphs for each species displayed 
  {
    
    # Returns graph of ASV relative abundance, use a ps with multiple timepoints
    # with annotations at the taxon level of interest and choose asv of interest
    # Does not include statistics
    asvRelAbDistributionTimeline <- function(ps, asv, taxon, group, time, custom_colors,displayASVNumber = TRUE){
      
      # Relative abundance of otu_table
      ps <- transformCounts(ps, transformation = "rel_ab")
      asv_ab <- t(otu_table(ps)[asv])
      asv_ab <- merge(asv_ab, sample_data(ps), by = 'row.names') # Bind metadata information
      taxa <- as.data.frame(tax_table(ps))
      p <- ggplot(data = asv_ab, aes(x = .data[[time]], y = .data[[asv]], group = .data[[group]], color = .data[[group]])) +
        stat_summary(fun.data = mean_se, geom = "errorbar", width = 1, aes(color = .data[[group]]))+
        stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
        stat_summary(fun = mean, geom = "point", size = 1, color = "black") +
        scale_color_manual(values = custom_colors)+
        scale_fill_manual(values = custom_colors)+
        labs(y = "Relative abundance (%)", title = paste(ifelse(taxon == "Species",paste(taxa[asv,"Genus"],taxa[asv,"Species"]), taxa[asv,taxon]), ifelse(displayASVNumber, asv, "")))
      print(p+my_theme())
      
    }
    sample_data(ps_flt_all)$gg_group2 <- factor(sample_data(ps_flt_all)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
    sample_data(ps_flt_all)$timepoint <- as.numeric(as.character(sample_data(ps_flt_all)$timepoint))
    
    # Manually add species names that we identified with BLASTn
    tax_table(ps_flt_all)["ASV8","Species"] <- "murinus (ASV8)"
    tax_table(ps_flt_all)["ASV38","Species"] <- "murinus (ASV38)"
    tax_table(ps_flt_all)["ASV19","Species"] <- "intestinale"
    tax_table(ps_flt_all)["ASV6","Species"] <- "rodentium (ASV6)"
    tax_table(ps_flt_all)["ASV1","Species"] <- "rodentium (ASV1)"
    
    timeline_theme <- theme(
      axis.text.x    = element_text(size = 7),
      axis.title.x   = element_text(size = 10, margin = margin(t = 10)),
      axis.title.y   = element_text(size = 10, margin = margin(r = 10)),
      plot.title     = element_text(size = 10, margin = margin(b = 10)),
      plot.margin    = margin(5, 5, 5, 5)    # optional: equalize outer margins
    )

    
    # Species M. intestinal ASV19
    asv19 <- asvRelAbDistributionTimeline(ps_flt_all, "ASV19", taxon = "Species",
                                 "gg_group2", "timepoint", c("#95BECF","#F2AA84","#325B6C","#B22222"),
                                 displayASVNumber = FALSE)+
      labs(color = "Group", x = "Days")+
      scale_x_continuous(
        breaks = seq(min(sample_data(ps_flt_all)$timepoint), max(sample_data(ps_flt_all)$timepoint), by = 7)                   # "T" before each label
      )+timeline_theme

    
    # Species ligilactobacillus murinus ASV 8/ ASV 38
    asv8 <- asvRelAbDistributionTimeline(ps_flt_all, "ASV8", taxon = "Species",
                                         "gg_group2", "timepoint", c("#95BECF","#F2AA84","#325B6C","#B22222"),
                                         displayASVNumber = FALSE)+
      labs(color = "Group", x = "Days")+
      scale_x_continuous(
        breaks = seq(min(sample_data(ps_flt_all)$timepoint), max(sample_data(ps_flt_all)$timepoint), by = 7)                   # "T" before each label
      )+timeline_theme
    
    asvRelAbDistributionTimeline(ps_flt_all, "ASV38", taxon = "Species",
                                 "gg_group2", "timepoint", c("#95BECF","#F2AA84","#325B6C","#B22222"),
                                 displayASVNumber = FALSE)+
      labs(color = "Group", x = "Days")+
      scale_x_continuous(
        breaks = seq(min(sample_data(ps_flt_all)$timepoint), max(sample_data(ps_flt_all)$timepoint), by = 7)                   # "T" before each label
      )+timeline_theme
    
    # F rodentium - ASV1
    asv1 <- asvRelAbDistributionTimeline(ps_flt_all, "ASV1", taxon = "Species",
                                        "gg_group2", "timepoint", c("#95BECF","#F2AA84","#325B6C","#B22222"),
                                        displayASVNumber = FALSE)+
      labs(color = "Group", x = "Days")+
      scale_x_continuous(
        breaks = seq(min(sample_data(ps_flt_all)$timepoint), max(sample_data(ps_flt_all)$timepoint), by = 7)                   # "T" before each label
      )+timeline_theme
    
    
    # Combine individual graphs into one with merged legend
    species_relab_timeline <- (asv1 / asv19 / asv8)+
      plot_layout(guides = "collect")
    species_relab_timeline
      
  }
  
  # Create figure 5
  fig5 <- (volcanoPlot50vs500dss+
             theme(legend.position = c(0.85, 0.8))) + wrap_elements(full = species_relab_timeline) +
    plot_layout(widths = c(1, 1))
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig5.png", plot = fig5, width = 12, height = 6, dpi = 500)
  
}


"#1f78b4"
"#e31a1c"
"#084594"
"#99000d"

