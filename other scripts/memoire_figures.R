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
  
  # Mice DAI evolution
  adult_dss_followup <- read.csv("combined_adult_dss_followup.csv", header = TRUE, sep = ";")
  adult_dss_followup <- dssFollowupManipulation(df = adult_dss_followup,groupInfoCols = 4,dateStart = "2024-05-29",nbrDays = 5, negativeOnly = TRUE)
  adult_dssflwup_plot <- dssDiseaseIndexPlot(adult_dss_followup, statBarLast = TRUE, signifPositionShift = 0.5)
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
    labs(title = "Colon length at end timepoint (T112)")
  
  # Stats
  verifyStatsAssumptions(df = dissec_adult, group = "gg_group", measure = "colon_length")
  pairwise.wilcox.test(dissec_adult$colon_length, dissec_adult$gg_group, p.adjust.method = "fdr")
}
