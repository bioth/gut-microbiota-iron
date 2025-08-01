# loading libraries
{
  library("tidyverse") # loading bunch of packages
  library("ggplot2") # come on, everyone knows what it is used for
  library("dplyr") # arranging and manipulating data easily
  # library("geepack") # library for loading GEE tests
  library("lme4") # library for loading ANOVA
  library("car") # for anova too
  library("ggsignif") # adding significance bars to ggplots
  library("readxl") # Read and write excel files
  # library("PMCMRplus") # gamesHowellTest (post hoc test for welch test when variances are unequal)
  # library("FSA")
  library(permuco)
}


# Setting working directory
setwd("~/Documents/CHUM_git/gut-microbiota-iron/")
setwd("I:/Chercheurs/Santos_Manuela/Thibault M/gut-microbiota-iron/")
setwd("~/CHUM_git/gut-microbiota-iron/")

# Loading functions for data manipulation
source("other scripts/dataManipFunctions.R")
source("pipeline/microbiota_analysis/utilities.R")

# Young mice + DSS experiment
{
  setwd("experiments/finished exp/young-DSS-exp3")
  young_weight <- read.csv("young48_weight_cageChanges.csv", header = TRUE, sep = ";")
  young_weight <- young_weight[,-c(6:10,18:ncol(young_weight))] # Remove repeated weight measures at start of experiment + measures after dss
  young_weight <- weightDataManipulation(young_weight, groupInfoCols = 4, fromDay0 = TRUE)
  young_weight$time_numeric[young_weight$time_numeric == 6] <- 7
  
  # Scatter plot with weight evolution over time for diet exposure
  young_weight_plot <- weightPlot(young_weight, percentage = TRUE, diet_only = TRUE, title = "Body weight change during diet exposure")
  young_weight_plot+
    scale_x_continuous(
      breaks = seq(min(young_weight$time_numeric), max(young_weight$time_numeric), by = 7), # breaks every 7 days
      labels = function(x) paste0("T", x)                    # "T" before each label
    )+
    labs(x = "Timepoint")
  ggsave("../../../figures/youngDSS/diet_weight.png", bg = "white",height = 5, width =6, dpi = 300)
  
  young_weight <- read.csv("young48_weight_cageChanges.csv", header = TRUE, sep = ";")
  young_weight <- young_weight[,-c(5:16)] # Remove weight measures before start of DSS
  young_weight <- weightDataManipulation(young_weight, groupInfoCols = 4, fromDay0 = TRUE)
  
  # Scatter plot with the four different treatments (diet combined with dss or control)
  young_weight_plot <- weightPlot(young_weight, percentage = TRUE, diet_only = FALSE, title = "Body weight")
  young_weight_plot+
    scale_x_continuous(
    breaks = seq(min(young_weight$time_numeric), max(young_weight$time_numeric), by = 7), # breaks every 7 days
    labels = function(x) paste0("T", x)                    # "T" before each label
  )+
    labs(x = "Timepoint")+
    theme(axis.ticks.x = element_line(),
          legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5),
          axis.text.x = element_text(size = 6))
  ggsave("../../../figures/youngDSS/weight.png", bg = "white",height = 5, width =6, dpi = 300)
  
  
  young_weight <- read.csv("young48_weight_cageChanges.csv", header = TRUE, sep = ";")
  young_weight <- young_weight[,-c(6:10)] # Remove weight measures before start of DSS
  young_weight <- weightDataManipulation(young_weight, groupInfoCols = 4, fromDay0 = TRUE)
  
  # Scatter plot with the four different treatments (diet combined with dss or control) - all timepoints
  young_weight_plot <- weightPlot(young_weight, percentage = TRUE, diet_only = FALSE, title = "Body weight", nbrDaysStart = 50)
  young_weight_plot
  ggsave("../../../figures/youngDSS/dss_weight.png",
         young_weight_plot, bg = "white",height = 5, width =8, dpi = 300)
  
  # Mice DAI evolution 
  young_dss_followup <- read.csv("young48_dss_followup.csv", header = TRUE, sep=";")
  young_dss_followup <- dssFollowupManipulation(df = young_dss_followup,groupInfoCols = 4,dateStart = "2024-05-29",nbrDays = 5, negativeOnly = TRUE)
  young_dssflwup_plot <- dssDiseaseIndexPlot(young_dss_followup, statBarLast = TRUE, signifAnnotation = "*")
  young_dssflwup_plot+
    theme(axis.ticks.x = element_line(),
          legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5))
  ggsave("../../../figures/youngDSS/dai.png", bg = "white",height = 4, width =5, dpi = 300)
  
# Statistics for last day of DAI
  verifyStatsAssumptions(df = young_dss_followup[young_dss_followup$time_numeric == "5",], group = "gg_group", measure = "index")
  wilcox.test(index ~ gg_group, data = young_dss_followup[young_dss_followup$time_numeric == "5",])
  
  # Mice body weight changes during DSS (scatter plot)
  young_dss_weight <- read.csv("young48_weights_dss.csv", header = TRUE, sep = ";")
  young_dss_weight <- weightDataManipulation(young_dss_weight, groupInfoCols = 4, fromDay0 = TRUE)
  weightPlot(young_dss_weight, percentage = TRUE, diet_only = TRUE)
  
  # Mice dissection data
  dissec_young <- read.csv("young48_dss_dissection.csv", sep = ";", header = TRUE)
  dissec_young <- dissectionDataManipulation(dissec_young, groupInfoCols = 4, numerical = FALSE)
  dissec_young$gg_group <- factor(dissec_young$gg_group, levels = c("50 water","500 water","50 dss","500 dss"))
  dissec_young$colon_length_nrm <- dissec_young$colon_length/dissec_young$body_weight
  
  #boxplot for body weight
  young_dissec_bw <- dissecBoxplot(dissec_young,"body", stats = TRUE,
                                   test_results = c("*","n.s.","n.s.","p=0.059"),
                                   text_sizes = c(5,3,3,3), upper_margin = 2) 
  ggsave("../../../figures/youngDSS/bw_final.png", bg = "white",height = 4, width =4, dpi = 300)
  
  # Statistics
  verifyStatsAssumptions(df = dissec_young, group = "gg_group", measure = "body_weight")
  anova <- aov(body_weight ~ diet * treatment, data = dissec_young) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  
  #boxplot for std liver weight
  dissecBoxplot(dissec_young,"liver",stats = TRUE, all.ns = TRUE) 
  ggsave("../../../figures/youngDSS/lvr_final.png", bg = "white",height = 4, width =4, dpi = 300)
  
  # Statistics
  verifyStatsAssumptions(df = dissec_young, group = "gg_group", measure = "std_liver_weigth") 
  anova <- aov(std_liver_weigth ~ diet * treatment, data = dissec_young) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  
  # boxplot for std spleen weight
  dissecBoxplot(dissec_young,"spleen", stats = TRUE, all.ns = TRUE) 
  ggsave("../../../figures/youngDSS/spln_final.png", bg = "white",height = 4, width =4, dpi = 300)
  
  # Statistics
  verifyStatsAssumptions(df = dissec_young, group = "gg_group", measure = "std_spleen_weigth")
  pairwise.wilcox.test(dissec_young$std_spleen_weigth, dissec_young$gg_group, p.adjust.method = "BH")
  
  #boxplot for colon length (non std)
  dissecBoxplot(dissec_young,"colon", stats = TRUE, all.ns = TRUE) 
  ggsave("../../../figures/youngDSS/cln_final.png", bg = "white",height = 4, width =4, dpi = 300)
  
  # Statistics
  verifyStatsAssumptions(df = dissec_young, group = "gg_group", measure = "colon_length")
  anova <- aov(colon_length ~ diet * treatment, data = dissec_young) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  
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
  p <- p+
  scale_x_discrete(labels = c("50 ppm","500 ppm"))+
  labs(y = "µg iron/g of stool", title = "Feces iron concentration\nat end of diet exposure")+
  guides(color = "none")
  p
  ggsave("../../../figures/youngDSS/t35_iron.png",
         p, bg = "white",height = 4, width =4, dpi = 300)
  
  # Stats 
  verifyStatsAssumptions(df, "diet" , "iron_concentration")
  wilcox.test(iron_concentration ~ diet, data = df)
  
  #Tf ferrozine assay for stools
  sheets <- read_excel_allsheets("young48_dss_ferrozine_tF.xlsx")
  df <- as.data.frame(sheets["Ferrozine Stool Tfinal"])
  df <- df[4:53,c(3,14:16)]
  colnames(df) <- df[1,]
  colnames(df)[1:2] <- c("id","iron_concentration")
  df <- df[-c(1,26,48),]
  df$gg_group <- paste0(df$diet, ":", df$treatment)
  df$gg_group <- factor(df$gg_group, levels = c("50:water","500:water","50:dss","500:dss"))
  df$iron_concentration <- as.numeric(df$iron_concentration)
  
  p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Feces iron concentration\nat final day", y_axis_title = "yg of iron per g of stools", custom_colors = c("blue","red","darkblue", "darkred"),
                  stats = TRUE, all.ns = TRUE)
  p <- p+
  scale_x_discrete(labels = c("50 ppm\nCtrl","500 ppm\nCtrl","50 ppm\nDSS","500 ppm\nDSS"))+
  labs(y = "µg iron per g of liver")+
  guides(color = "none")
  p
  ggsave("../../../figures/youngDSS/tFinal_iron.png",
         p, bg = "white",height = 4, width =4, dpi = 300)
  
  # Stats 
  verifyStatsAssumptions(df, "gg_group" , "iron_concentration")
  pairwise.wilcox.test(df$iron_concentration, df$gg_group, p.adjust.method = "BH")
  
  
  df$diet <- factor(df$diet, levels = c("50","500"))
  df$treatment <- factor(df$treatment, levels = c("water","dss"))
  pairwise_permuco(df, "iron_concentration", "diet", "treatment")
  
  
  #Tf ferrozine assay for liver
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
  
  p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Liver iron concentration\nat final day", y_axis_title = "yg of iron", custom_colors = c("blue","red","darkblue", "darkred"),
                  stats = TRUE, test_results = c("**","**","n.s.","n.s."), text_sizes = c(5,5,3,3), upper_margin = 30)
  p <- p+
  scale_x_discrete(labels = c("50 ppm\nCtrl","500 ppm\nCtrl","50 ppm\nDSS","500 ppm\nDSS"))+
  labs(y = "µg Fe/g dry weight")+
  guides(color = "none")
  p
  ggsave("../../../figures/youngDSS/lvr_iron.png",
         p, bg = "white",height = 4, width =4, dpi = 300)
  
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
  
  p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Spleen iron concentration\nat final day", y_axis_title = "yg of iron", custom_colors = c("blue","red","darkblue", "darkred"),
                  stats = TRUE, test_results = c("*","n.s.","*","*"), text_sizes = c(5,3,5,5), upper_margin = 290)
  p+
  scale_x_discrete(labels = c("50 ppm\nCtrl","500 ppm\nCtrl","50 ppm\nDSS","500 ppm\nDSS"))+
  labs(y = "µg Fe/g dry weight")+
  guides(color = "none")
  ggsave("../../../figures/youngDSS/spln_iron.png", bg = "white",height = 4, width =4, dpi = 300)
  
  
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
}

# PCR validation of B. pseudolongum and F rodentium in young mice + DSS experiment
{
  # Load pcr df
  df <- read_excel("experiments/finished exp/young-DSS-exp3/pcr_validation/RT-PCR B pseudolongum Primers2.xlsx")
  colnames(df) <- df[3,]
  colnames(df)[c(1,4,6,7,8)] <- c("id","frod","bpse2","pseudolongum","bpse") # T.C. pair of primers B. psed
  df <- df[-c(1:3),]
  
  # Load metadata and bind information to dataframe
  meta <- read.csv("experiments/finished exp/young-DSS-exp3/young48_dss_dissection.csv", sep = ";")
  meta$id <- substring(meta$id, first = 1, last = 5)
  df <- merge(df, meta, by = "id")
  df$diet <- factor(df$diet, levels = c("50", "500"))
  df$treatment <- factor(df$treatment, levels = c("water", "dss"))
  df$gg_group <- factor(paste0(df$diet, ":", df$treatment), levels = c("50:water", "500:water", "50:dss", "500:dss"))
  df$frod <- as.numeric(df$frod)
  df$bpse <- as.numeric(df$bpse)
  df$bpse2 <- as.numeric(df$bpse2)
  
  p <- ironBoxplot(df, "frod", group = "gg_group", title = "F. rodentium qPCR abundance", y_axis_title = "F. rodentium/16S", custom_colors = c("blue","red","darkblue","darkred"))
  p
  
  pairwise_permuco(df, "frod", "diet", "treatment")
  pairwise_permuco(df, "bpse", "diet", "treatment")
  pairwise_permuco(df, "bpse2", "diet", "treatment")
  
  p <- ironBoxplot( df[df$treatment == "dss",], "frod", group = "diet", title = "F. rodentium qPCR abundance", y_axis_title = "F. rodentium/16S", custom_colors = c("darkblue","darkred"))
  p
  
  # Stats
  verifyStatsAssumptions(df = df, group = "gg_group", measure = "frod") 
  anova <- aov(frod ~ diet * treatment, data = df) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  pairwise.wilcox.test(df$frod, df$gg_group, p.adjust.method = "BH")
  
  # DSS groups only
  p <- ironBoxplot(df[df$treatment == "dss",], "frod", group = "diet", title = "F. rodentium qPCR abundance", y_axis_title = "F. rodentium /16S", custom_colors = c("darkblue","darkred"), stats = FALSE)
  p
  
  verifyStatsAssumptions(df = df[df$treatment == "dss",], group = "diet", measure = "frod") 
  t.test(frod ~ diet, data = df[df$treatment == "dss",], var.equal = TRUE)
  
  p <- ironBoxplot(df[df$treatment == "dss",], "bpse", group = "diet", title = "B. pseudolongum qPCR abundance", y_axis_title = "B. pseudolongum /16S", custom_colors = c("darkblue","darkred"), stats = FALSE)
  p
  
  verifyStatsAssumptions(df = df[df$treatment == "dss",], group = "diet", measure = "bpse") 
  wilcox.test(bpse ~ diet, data = df[df$treatment == "dss",])
  
  p <- ironBoxplot(df[df$treatment == "dss",], "bpse2", group = "diet", title = "B. pseudolongum qPCR abundance", y_axis_title = "B. pseudolongum /16S", custom_colors = c("darkblue","darkred"), stats = FALSE)
  p
  
  verifyStatsAssumptions(df = df[df$treatment == "dss",], group = "diet", measure = "bpse2") 
  wilcox.test(bpse2 ~ diet, data = df[df$treatment == "dss",])
  
}

# Young mice + Abx experiment
{
  setwd("experiments/finished exp/young-abx-exp6/")

  # Body weight measurements for diet exposure
  young_weight <- read.csv("young48-abx-bw.csv", sep = ";")
  young_weight <- young_weight[,c(1:12)]
  young_weight <- weightDataManipulation(young_weight, groupInfoCols = 4, fromDay0 = TRUE)
  
  # Scatter plot with weight evolution over time for diet exposure
  young_weight_plot <- weightPlot(young_weight, percentage = TRUE, diet_only = TRUE)
  young_weight_plot
  
  # Body weight measurements after abx exposure
  young_weight <- read.csv("young48-abx-bw.csv", sep = ";")
  young_weight <- young_weight[,c(1:4,12:21)]
  young_weight <- weightDataManipulation(young_weight, groupInfoCols = 4, fromDay0 = TRUE)
  
  # Scatter plot with weight evolution over time for abx exposure
  young_weight_plot <- weightPlot(young_weight, percentage = FALSE, diet_only = FALSE, abxExp = TRUE)
  young_weight_plot
  
  # Stats
  young_weight$gg_group <- factor(young_weight$gg_group, levels = c("50 water", "500 water", "50 abx", "500 abx"))
  timePoints <- unique(young_weight$time_numeric)
  timePoints <- timePoints[-1] # Remove baseline timepoint because all values are the same
  for(t in timePoints){
    verifyStatsAssumptions(df = young_weight[young_weight$time_numeric == t,], group = "gg_group", measure = "weight")
  }
  
  for(t in timePoints){
    print(pairwise.wilcox.test(young_weight[young_weight$time_numeric == t,]$weight, young_weight$gg_group, p.adjust.method = "BH"))
    
  }
  
  
  verifyStatsAssumptions(df = young_weight, group = "gg_group", measure = "body_weight")
  
  # Mice dissection data
  dissec_young <- read_excel("dissection.xlsx")
  dissec_young <- dissec_young[-17,] # Remove empty rows
  colnames(dissec_young)[5:8] <- str_replace(colnames(dissec_young)[5:8], pattern = " ", replacement = "_") # Replace spaces by underscore for variables colnames
  dissec_young <- dissectionDataManipulation(dissec_young,groupInfoCols = 3)
  dissec_young$gg_group <- factor(dissec_young$gg_group, levels = c("50 water","500 water","50 abx","500 abx"))
  dissec_young$colon_length_nrm <- dissec_young$colon_length/dissec_young$body_weight
  
  #boxplot for body weight
  young_dissec_bw <- dissecBoxplot(dissec_young,"body", abxExp = TRUE) 
  young_dissec_bw
  
  # Stats
  verifyStatsAssumptions(df = dissec_young, group = "gg_group", measure = "body_weight")
  anova <- aov(body_weight ~ diet * treatment, data = dissec_young) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  
  #boxplot for std liver weight
  young_dissec_lvr <- dissecBoxplot(dissec_young,"liver", abxExp = TRUE) 
  young_dissec_lvr
  
  # Stats
  verifyStatsAssumptions(df = dissec_young, group = "gg_group", measure = "std_liver_weigth")
  anova <- aov(std_liver_weigth ~ diet * treatment, data = dissec_young) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  
  # Boxplot for std spleen weight
  young_dissec_spln <- dissecBoxplot(dissec_young,"spleen", abxExp = TRUE) 
  young_dissec_spln
  
  # Stats
  verifyStatsAssumptions(df = dissec_young, group = "gg_group", measure = "std_spleen_weigth")
  anova <- aov(std_spleen_weigth ~ diet * treatment, data = dissec_young) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  
  # Boxplot for colon length (non std)
  young_dissec_cln <- dissecBoxplot(dissec_young,"colon", abxExp = TRUE) 
  young_dissec_cln
  
  # Stats
  verifyStatsAssumptions(df = dissec_young, group = "gg_group", measure = "colon_length")
  anova <- aov(colon_length ~ diet * treatment, data = dissec_young) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  
  # Ferrozine measurements
  # Load metadata to retrieve group information
  meta <- read_xlsx("dissection.xlsx")
  meta <- meta[,c(1:4)]
  meta$ID <- substring(meta$ID, first = 1, last = 5)
  meta <- meta[-17,] # Remove dead mouse
  
  #T35 ferrozine assay for stools
  df <- read_xlsx("t35_stool_ferrozine.xlsx")
  df <- df[,c(1,12)]
  colnames(df) <- df[1,]
  colnames(df)[1:2] <- c("ID","iron_concentration")
  df <- df[-c(1),] 
  df <- df[!is.na(df$iron_concentration),]  # Remove dead mice
  df$iron_concentration <- as.numeric(df$iron_concentration)
  df <- merge(df, meta, by ="ID") # Bind metadata information to df
  df$diet <- factor(df$diet, levels = c("50","500"))
  
  p = ironBoxplot(df, "iron_concentration", group = "diet", title = "Iron concentration in stools at day 35", y_axis_title = "µg Fe/g of stool", custom_colors = c("blue","red"))
  p <- p+
    scale_x_discrete(labels = c("50 ppm","500 ppm"))+
    labs(y = "µg Fe/g of stool", title = "")+
    ylim(0,NA)+
    guides(color = "none")+
    theme(panel.grid.major = element_blank(),
          axis.text.y = element_text(face = "bold", size = 12),
          axis.text.x = element_text(face = "bold", size = 10))
  p
  
  # Stats
  verifyStatsAssumptions(df = df, group = "diet", measure = "iron_concentration")
  wilcox.test(iron_concentration ~ diet, data = df)
  
  #TF ferrozine assay for stools
  df <- read_xlsx("tF_stool_ferrozine.xlsx")
  df <- df[,c(1,12)]
  colnames(df)[1:2] <- c("ID","iron_concentration")
  df <- df[-c(1),] 
  df <- df[!is.na(df$iron_concentration),]  # Remove empty rows
  df$iron_concentration <- as.numeric(df$iron_concentration)
  df <- merge(df, meta, by ="ID") # Bind metadata information to df
  df$diet <- factor(df$diet, levels = c("50","500"))
  df$treatment <- factor(df$treatment, levels = c("water","abx"))
  df$gg_group <- paste0(df$diet, ":", df$treatment)
  df$gg_group <- factor(df$gg_group, levels = c("50:water", "500:water", "50:abx", "500:abx"))
  
  p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Iron concentration in stools at final day", y_axis_title = "µg Fe/g of stool", custom_colors = c("blue","red", "deepskyblue", "brown1"))
  p <- p+
    scale_x_discrete(labels = c("50 ppm\nCtrl","500 ppm\nCtrl","50 ppm\nAbx","500 ppm\nAbx"))+
    labs(y = "µg Fe/g of stool", title = "")+
    ylim(0,NA)+
    guides(color = "none")+
    theme(panel.grid.major = element_blank(),
          axis.text.y = element_text(face = "bold", size = 12),
          axis.text.x = element_text(face = "bold", size = 10))
  p
  
  # Stats
  verifyStatsAssumptions(df = df, group = "gg_group", measure = "iron_concentration")
  pairwise.wilcox.test(df$iron_concentration, df$gg_group, p.adjust.method = "BH")
  
  # Tfinal ferrozine assay for liver
  df <- read_excel("liver_ferrozine.xlsx")
  df <- df[,c(3,14)]
  colnames(df)[1:2] <- c("ID","iron_concentration")
  df <- df[-c(1:3),] 
  df$iron_concentration <- as.numeric(df$iron_concentration)
  df <- merge(df, meta, by ="ID") # Bind metadata information to df
  df$diet <- factor(df$diet, levels = c("50","500"))
  df$treatment <- factor(df$treatment, levels = c("water", "abx"))
  df$gg_group <- factor(paste0(df$diet, ":", df$treatment), levels = c("50:water", "500:water", "50:abx", "500:abx"))
  
  p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Total iron in liver", y_axis_title = "µg Fe/g of liver", custom_colors = c("blue","red", "deepskyblue", "brown1"))
  p <- p+
    scale_x_discrete(labels = c("50 ppm\nCtrl","500 ppm\nCtrl","50 ppm\nAbx","500 ppm\nAbx"))+
    labs(y = "µg Fe/g dry weight", title = "")+
    ylim(0,NA)+
    guides(color = "none")+
    theme(panel.grid.major = element_blank(),
          axis.text.y = element_text(face = "bold", size = 12),
          axis.text.x = element_text(face = "bold", size = 10))
  p
  
  
  # Stats
  verifyStatsAssumptions(df = df, group = "gg_group", measure = "iron_concentration")
  pairwise.t.test(df$iron_concentration, df$gg_group, p.adjust.method = "BH", pool.sd = FALSE)
  
  # Spleen iron concentration at tfinal
  df <- read_excel("spleen_ferrozine.xlsx")
  df <- df[,c(3,14)]
  colnames(df)[1:2] <- c("ID","iron_concentration")
  df <- df[-c(1:3),] 
  df$iron_concentration <- as.numeric(df$iron_concentration)
  df <- merge(df, meta, by ="ID") # Bind metadata information to df
  df$diet <- factor(df$diet, levels = c("50","500"))
  df$treatment <- factor(df$treatment, levels = c("water", "abx"))
  df$gg_group <- factor(paste0(df$diet, ":", df$treatment), levels = c("50:water", "500:water", "50:abx", "500:abx"))
  
  p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Total iron in liver", y_axis_title = "µg Fe/g of liver", custom_colors = c("blue","red", "deepskyblue", "brown1"))
  p <- p+
    scale_x_discrete(labels = c("50 ppm\nCtrl","500 ppm\nCtrl","50 ppm\nAbx","500 ppm\nAbx"))+
    labs(y = "µg Fe/g dry weight", title = "")+
    ylim(0,NA)+
    guides(color = "none")+
    theme(panel.grid.major = element_blank(),
          axis.text.y = element_text(face = "bold", size = 12),
          axis.text.x = element_text(face = "bold", size = 10))
  p
  
  # Stats
  verifyStatsAssumptions(df = df, group = "gg_group", measure = "iron_concentration")
  # Assuming variances are unequal
  
  
  
  
  anova <- aov(iron_concentration ~ diet * treatment, data = df) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  
}

# PCR validation of F rodentium in young mice + Abx experiment
{
  # Load pcr df
  df <- read_excel("experiments/finished exp/young-abx-exp6/RT-PCR F rodentium.xlsx")
  colnames(df) <- df[3,]
  colnames(df)[c(1,4)] <- c("id","frod") # F. rodentium primer pair
  df <- df[-c(1:3),]
  
  # Load metadata and bind information to dataframe
  meta <- read.xlsx("experiments/finished exp/young-abx-exp6/dissection.xlsx")
  colnames(meta)[2] <- "id"
  meta$id <- substring(meta$id, first = 1, last = 5)
  df <- merge(df, meta, by = "id")
  df$diet <- factor(df$diet, levels = c("50", "500"))
  df$treatment <- factor(df$treatment, levels = c("water", "abx"))
  df$gg_group <- factor(paste0(df$diet, ":", df$treatment), levels = c("50:water", "500:water", "50:abx", "500:abx"), labels = c("50 ppm\nCtrl", "500 ppm\nCtrl", "50 ppm\nAbx", "500 ppm\nAbx"))
  df$frod <- as.numeric(df$frod)
  df <- df[-c(37,34,18),]# Remove outliers
  
  p <- ironBoxplot(df, "frod", group = "gg_group", title = "F. rodentium qPCR abundance", y_axis_title = "F. rodentium/16S", custom_colors = c("blue", "red", "deepskyblue", "brown1"), stats = TRUE,
                   test_results = c("n.s.","P=0.06","P=0.07","n.s."), text_sizes = c(3,3,3,3), upper_margin = 0.1)
  p
  ggsave("figures/youngAbx/frod_pcr_tfinal.png", bg = "white",height = 4, width =4, dpi = 300)
  
  pairwise_permuco(df, "frod", "diet", "treatment")
  
  p <- ironBoxplot(df[df$treatment == "abx",], "frod", group = "diet", title = "F. rodentium qPCR abundance", y_axis_title = "F. rodentium/16S", custom_colors = c("deepskyblue", "brown1"), stats = FALSE)
  p
  
  wilcox.test(frod ~ diet, data = df[df$treatment == "abx",])
  
  
  # Stats
  verifyStatsAssumptions(df = df, group = "gg_group", measure = "frod") 
  anova <- aov(frod ~ diet * treatment, data = df) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  pairwise.wilcox.test(df$frod, df$gg_group, p.adjust.method = "BH")
  
  
  verifyStatsAssumptions(df = df[df$treatment == "abx",], group = "diet", measure = "frod") 
  wilcox.test(frod ~ diet, data = df[df$treatment == "dss",], var.equal = TRUE)
  
  p <- ironBoxplot(df, "bpse", group = "gg_group", title = "B. pseudolongum qPCR abundance", y_axis_title = "B. pseudolongum /16S", custom_colors = c("blue","red","darkblue","darkred"))
  p+ylim(NA,NA)
  
  p <- ironBoxplot(df[df$treatment == "dss",], "bpse", group = "diet", title = "B. pseudolongum qPCR abundance", y_axis_title = "B. pseudolongum /16S", custom_colors = c("darkblue","darkred"))
  p
  
  # Removing two outliers
  df_sub <- df[-c(26,37),]
  
  p <- ironBoxplot(df_sub[df_sub$treatment == "dss",], "bpse", group = "diet", title = "B. pseudolongum qPCR abundance", y_axis_title = "B. pseudolongum /16S", custom_colors = c("darkblue","darkred"))
  p
  
  # Stats
  verifyStatsAssumptions(df = df, group = "gg_group", measure = "bpse")
  pairwise.wilcox.test(df$bpse, df$gg_group, p.adjust.method = "BH")
  
  verifyStatsAssumptions(df = df_sub, group = "gg_group", measure = "bpse")
  pairwise.wilcox.test(df_sub$bpse, df_sub$gg_group, p.adjust.method = "BH")
  
  
  verifyStatsAssumptions(df = df[df$treatment == "dss",], group = "diet", measure = "bpse") 
  wilcox.test(bpse ~ diet, data = df[df$treatment == "dss",])
  
  
  
}

# Adult mice + DSS experiment
{
  setwd("experiments/finished exp/adults-all-exp/")
  
  # Mice DAI evolution
  adult_dss_followup <- read.csv("combined_adult_dss_followup.csv", header = TRUE, sep = ";")
  adult_dss_followup <- dssFollowupManipulation(df = adult_dss_followup,groupInfoCols = 4,dateStart = "2024-05-29",nbrDays = 5, negativeOnly = TRUE)
  adult_dssflwup_plot <- dssDiseaseIndexPlot(adult_dss_followup)
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
  
  # Boxplot for body weight
  adult_dissec_bw <- dissecBoxplot(dissec_adult,"body") 
  adult_dissec_bw
  
  # Stats
  verifyStatsAssumptions(df = dissec_adult, group = "gg_group", measure = "body_weight")
  pairwise.wilcox.test(dissec_adult$body_weight, dissec_adult$gg_group, p.adjust.method = "BH")
  
  sub_df <- dissec_adult[!dissec_adult$ID== "36697B",] # Sub df removing the mouse for which they were no liver and spleen measures
  
  # Boxplot for std liver weight
  adult_dissec_lvr <- dissecBoxplot(sub_df,"liver") 
  adult_dissec_lvr
  
  # Stats
  verifyStatsAssumptions(df = dissec_adult, group = "gg_group", measure = "std_liver_weigth")
  pairwise.wilcox.test(sub_df$std_liver_weigth, sub_df$gg_group, p.adjust.method = "BH")
  
  # Boxplot for std spleen weight
  adult_dissec_spln <- dissecBoxplot(sub_df,"spleen") 
  adult_dissec_spln
  
  # Stats
  verifyStatsAssumptions(df = dissec_adult, group = "gg_group", measure = "std_spleen_weigth")
  pairwise.wilcox.test(sub_df$std_spleen_weigth, sub_df$gg_group, p.adjust.method = "BH")
  
  # Boxplot for colon length (non std)
  adult_dissec_cln <- dissecBoxplot(dissec_adult,"colon") 
  adult_dissec_cln
  
  # Stats
  verifyStatsAssumptions(df = dissec_adult, group = "gg_group", measure = "colon_length")
  pairwise.wilcox.test(dissec_adult$colon_length, dissec_adult$gg_group, p.adjust.method = "BH")
  
  
}

# Adult mice + Abx experiment
{
  setwd("experiments/finished exp/adults-all-exp/")
  
  # Dissection data
  dissec_adult <- read_excel("all-adults-dissection.xlsx")
  dissec_adult <- dissec_adult[dissec_adult$treatment != "dss",] # Remove dss mice
  dissec_adult <- dissec_adult[!is.na(dissec_adult$`body weight`),] # Remove dead mice (no ecorded values for final body weight)
  dissec_adult <- dissec_adult[!dissec_adult$batch =="C",] # Remove failed batch for abx experiment
  colnames(dissec_adult)[6:9] <- str_replace(colnames(dissec_adult)[6:9], pattern = " ", replacement = "_") # Replace spaces by underscore for variables colnames
  dissec_adult <- dissectionDataManipulation(dissec_adult, groupInfoCols = 4)
  dissec_adult$gg_group <- factor(dissec_adult$gg_group, levels = c("50 water","500 water","50 abx","500 abx"))
  dissec_adult$colon_length_nrm <- dissec_adult$colon_length/dissec_adult$body_weight
  
  # Boxplot for body weight
  adult_dissec_bw <- dissecBoxplot(dissec_adult,"body", abxExp = TRUE) 
  adult_dissec_bw
  
  # Stats
  verifyStatsAssumptions(df = dissec_adult, group = "gg_group", measure = "body_weight")
  anova <- aov(body_weight ~ diet * treatment, data = dissec_adult) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  
  # Boxplot for std liver weight
  adult_dissec_lvr <- dissecBoxplot(dissec_adult,"liver", abxExp = TRUE) 
  adult_dissec_lvr
  
  # Stats
  verifyStatsAssumptions(df = dissec_adult, group = "gg_group", measure = "std_liver_weigth")
  pairwise.wilcox.test(dissec_adult$std_liver_weigth, dissec_adult$gg_group, p.adjust.method = "BH")
  
  # Boxplot for std spleen weight
  adult_dissec_spln <- dissecBoxplot(dissec_adult,"spleen", abxExp = TRUE) 
  adult_dissec_spln
  
  # Stats
  verifyStatsAssumptions(df = dissec_adult, group = "gg_group", measure = "std_spleen_weigth")
  pairwise.wilcox.test(dissec_adult$std_spleen_weigth, dissec_adult$gg_group, p.adjust.method = "BH")
  
  # Boxplot for colon length (non std)
  adult_dissec_cln <- dissecBoxplot(dissec_adult,"colon", abxExp = TRUE) 
  adult_dissec_cln
  
  # Stats
  verifyStatsAssumptions(df = dissec_adult, group = "gg_group", measure = "colon_length")
  pairwise.wilcox.test(dissec_adult$colon_length, dissec_adult$gg_group, p.adjust.method = "BH")
  
  
}