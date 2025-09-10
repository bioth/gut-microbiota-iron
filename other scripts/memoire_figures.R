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
  library(ggbreak)
  library(lme4)
  library(lmerTest) # gives p-values for lmer
  library(emmeans)  # for post-hoc comparisons
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
  
  p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Stool iron concentration\nat end of experiment", y_axis_title = "yg of iron/g dry weight", custom_colors = c("#95BECF","#F2AA84","#325BAD","#B22222"),
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
    df$treatment <- factor(df$treatment, levels = c("water","dss"), labels = c("Ctrl","DSS"))
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
    df$gg_group <- factor(df$gg_group, levels = c("50:water","500:water","50:dss","500:dss"), labels = c(c("50:Ctrl","500:Ctrl","50:DSS","500:DSS")))
    df$iron_concentration <- as.numeric(df$iron_concentration)
    df$treatment <- factor(df$treatment, levels = c("water","dss"), labels = c("Ctrl","DSS"))
    dfTfinal <- df
    
    df_merged <- rbind(dfT35, dfTfinal)
    df_merged$gg_group 
    
    p = ironBoxplot(df_merged, "iron_concentration", group = "gg_group", title = "Stool iron", custom_colors = c("#95BECF","#F2AA84","#95BECF","#F2AA84","#325BAD","#B22222"),
                    stats = FALSE, shape = "treatment")
    p_combined <- p +
      guides(x = "axis_nested")+
      labs(y = "µg Fe/g dry weight")+
      theme(plot.title = element_text(size = 10),
            axis.title.y = element_text(size = 10))+
      geom_vline(xintercept = 2.5,
                 color      = "black",
                 linetype   = "dashed",
                 size       = 0.4)+
      ylim(0,7000)
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
    min(df_merged$iron_concentration)+ 4*min(df_merged$iron_concentration),
    min(df_merged$iron_concentration)+ 8*min(df_merged$iron_concentration),
    min(df_merged$iron_concentration)+ 12*min(df_merged$iron_concentration),
    min(df_merged$iron_concentration)+ 18*min(df_merged$iron_concentration))
  
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

  p_combined <- p_combined+
    scale_y_continuous(labels = scales::label_comma(big.mark = ","))

  # Tf ferrozine assay for liver
  df <- as.data.frame(sheets["Ferrozine Liver"])
  df <- df[1:51,c(3,6,8,14:18)]
  colnames(df) <- df[2,]
  colnames(df)[1:4] <- c("id","wet_weight","dry_weight","iron_concentration")
  df <- df[-c(1,2,27),]
  df$gg_group <- interaction(df$diet, df$treatment, sep=":")
  df$gg_group <- factor(df$gg_group, levels = c("50:water","500:water","50:dss","500:dss"), labels = c(c("50:Ctrl","500:Ctrl","50:DSS","500:DSS")))
  df$iron_concentration <- as.numeric(df$iron_concentration)
  df$liver_weight <- as.numeric(df$liver_weight)
  df$wet_weight <- as.numeric(df$wet_weight)
  df$dry_weight <- as.numeric(df$dry_weight)
  df$dry_to_wet_ratio <- df$dry_weight/df$wet_weight
  df <- df[-46,]
  df$treatment <- factor(df$treatment, levels = c("water","dss"), labels = c("Ctrl","DSS"))
  
  p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Liver iron", y_axis_title = "yg of iron", custom_colors = c("#95BECF","#F2AA84","#325BAD","#B22222"),
                  stats = TRUE, test_results = c("**","**","n.s.","n.s."), text_sizes = c(5,5,3,3), upper_margin = 25, vjustList = c(0.5,0.5,0,0), shape = "treatment")
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
  df$gg_group <- factor(df$gg_group, levels = c("50:water","500:water","50:dss","500:dss"), labels = c(c("50:Ctrl","500:Ctrl","50:DSS","500:DSS")))
  df$iron_concentration <- as.numeric(df$iron_concentration)
  df$spleen_weight <- as.numeric(df$spleen_weight)
  df$wet_weight <- as.numeric(df$wet_weight)
  df$dry_weight <- as.numeric(df$dry_weight)
  df$dry_to_wet_ratio <- df$wet_weight/df$dry_weight
  df <- df[-46,]
  df$treatment <- factor(df$treatment, levels = c("water","dss"), labels = c("Water","DSS"))
  
  p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Spleen iron", y_axis_title = "yg of iron", custom_colors = c("#95BECF","#F2AA84","#325BAD","#B22222"),
                  stats = TRUE, test_results = c("*","n.s.","*","*"), text_sizes = c(5,3,5,5), upper_margin = 400, vjustList = c(0.5,0,0.5,0.5), shape = "treatment")
  p4 <- p+
    guides(x = "axis_nested")+
    labs(y = "µg Fe/g dry weight")+
    theme(plot.title = element_text(size = 10),
          axis.title.y = element_text(size = 10))+
    scale_y_continuous(labels = scales::label_comma(big.mark = ","))+
    ylim(0,NA)
  p4

  # Stats 
  verifyStatsAssumptions(df, "gg_group" , "iron_concentration")
  bartlett.test(iron_concentration ~ gg_group, data = df)
  anova <- aov(iron_concentration ~ diet * treatment, data = df) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  
  kruskal.test(iron_concentration ~ gg_group, data = df)
  dunnTest(iron_concentration ~ gg_group, data = df, method = "bonferroni")
  
  df$diet <- factor(df$diet, levels = c("50","500"))
  pairwise_permuco(df, "iron_concentration", "diet", "treatment",np = 5000)
  pairwise_permuco2(df, "iron_concentration", "diet", "treatment",np = 5000)
  
  library(welchADF)
  
  wj <- welchADF.test(
    iron_concentration ~ diet * treatment,
    data     = df, 
    contrast = "omnibus"  # test all effects that appear in the formula
  )
  summary(wj)
  
  welchADF.test(iron_concentration ~ diet * treatment, data = df,
                contrast = "all.pairwise",
                effect   = c("diet"))
  
  welchADF.test(iron_concentration ~ diet * treatment, data = df,
                contrast = "all.pairwise",
                effect   = c("diet","treatment"))
  
  post_diet <- welchADF.test(
    iron_concentration ~ diet * treatment, data = df,
    contrast = "all.pairwise",
    effect = list(diet = levels(df$diet))
  )
  summary(post_diet)
  
  
  
  
  library(emmeans)
  library(sandwich)
  
  m  <- lm(iron_concentration ~ diet * treatment, data = df)
  
  # Robust covariance (HC3)
  V  <- vcovHC(m, type = "HC3")
  
  # Simple effects: diet within each treatment
  emm <- emmeans(m, ~ diet | treatment, vcov. = V)
  
  # Get 50 vs 500 within Water and within DSS, with BH adjustment
  pw  <- pairs(emm, adjust = "none") |> as.data.frame()
  pw$p_adj <- p.adjust(pw$p.value, method = "BH")
  pw[, c("treatment","contrast","estimate","SE","df","p.value","p_adj")]
  
  # planned comparisons
  comparisons <- list(
    "50Ctrl_vs_500Ctrl" = df$diet == "50" & df$treatment == "Ctrl",
    "50Ctrl_vs_500Ctrl_alt" = df$diet == "500" & df$treatment == "Ctrl",
    
    "50DSS_vs_500DSS"   = df$diet == "50" & df$treatment == "DSS",
    "50DSS_vs_500DSS_alt" = df$diet == "500" & df$treatment == "DSS",
    
    "Ctrl50_vs_DSS50"   = df$diet == "50" & df$treatment %in% c("Ctrl","DSS"),
    "Ctrl500_vs_DSS500" = df$diet == "500" & df$treatment %in% c("Ctrl","DSS")
  )
  
  # do Welch t-tests
  pvals <- c(
    t.test(iron_concentration ~ diet, data = subset(df, treatment=="Water"))$p.value,
    t.test(iron_concentration ~ diet, data = subset(df, treatment=="DSS"))$p.value,
    t.test(iron_concentration ~ treatment, data = subset(df, diet=="50"))$p.value,
    t.test(iron_concentration ~ treatment, data = subset(df, diet=="500"))$p.value
  )
  
  # adjust for multiple testing (FDR or Holm)
  p.adjust(pvals, method="fdr")   # or "BH"
  
  
  pairwise.wilcox.test(df$iron_concentration, df$gg_group, p.adjust.method = "fdr")
  
  
  # Combine plots
  fig1 <- wrap_elements(p1 + p2 +p3 +p4) +
    plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 20))
    
  fig1
  
  existingDirCheck("~/CHUM_git/figures/memoire/dss/")
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig1.png", plot = fig1, width = 5, height = 6, dpi = 500)
  
  # Another version with grouped stool iron
  fig1_v2 <- p_combined/(p3+p4)+
    plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 20))
  fig1_v2
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig1v2.png", plot = fig1_v2, width = 6, height = 5, dpi = 500)
  
  # Checking spleen weights
  {
    dissec_young <- read.csv("young48_dss_dissection.csv", sep = ";", header = TRUE)
    dissec_young <- dissectionDataManipulation(dissec_young, groupInfoCols = 4, numerical = FALSE)
    dissec_young$gg_group <- factor(dissec_young$gg_group, levels = c("50 water","500 water","50 dss","500 dss"))
    dissec_young$treatment <- factor(dissec_young$treatment, levels = c("water","dss"))
    
    verifyStatsAssumptions(dissec_young, "gg_group", "std_spleen_weigth")
    verifyStatsAssumptions(dissec_young, "gg_group", "spleen_weight")
    pairwise.wilcox.test(dissec_young$std_spleen_weigth, dissec_young$gg_group, p.adjust.method = "fdr")
    pairwise.wilcox.test(dissec_young$spleen_weight, dissec_young$gg_group, p.adjust.method = "fdr")
    
    # Boxplot for body weight of the spleen
    spleen = ironBoxplot(dissec_young, "spleen_weight", group = "gg_group", title = "At end of experiment", y_axis_title = "Spleen weight (g)", custom_colors = c("#95BECF","#F2AA84","#325BAD","#B22222"),
                     stats = TRUE, test_results = c("*","n.s.","*","n.s."), text_sizes = c(5,3,5,3),
                     upper_margin = 3, vjustList = c(0.5,0,0.5,0), tip_length = 0.1, shape = "treatment", scaleFactor = 30)
    spleen <- spleen +
      guides(x = "axis_nested")+
      theme(plot.title = element_text(size = 10),
            axis.title.y = element_text(size = 10))+
      ylim(0,NA)
    spleen
  }
  
}

# Figure 2 - dss induced colitis model
{
  # Set working directory
  setwd("experiments/finished exp/young-DSS-exp3")

  young_weight <- read.csv("young48_weight_cageChanges.csv", header = TRUE, sep = ";")
  young_weight <- young_weight[,-c(6:10)] # Remove weight measures before start of DSS
  young_weight <- weightDataManipulation(young_weight, groupInfoCols = 4, fromDay0 = TRUE)
  young_weight$gg_group <- factor(young_weight$gg_group, levels = c("50 water","50 DSS","500 water","500 DSS"),
                                  labels = c("50 ppm Ctrl","50 ppm DSS","500 ppm Ctrl","500 ppm DSS"))
  young_weight$timepoint <- factor(young_weight$time_numeric)
  
  young_weight50 <- young_weight[young_weight$diet == "50",]
  young_weight500 <- young_weight[young_weight$diet == "500",]
  
  # Scatter plot with 50ctrl vs 50dss - all timepoints
  young_weight50_plot <- weightPlot(young_weight50, group = "gg_group", percentage = FALSE, title = "Iron sufficient diet (50 ppm)", customColors = c("#95BECF","#325BAD"))
  young_weight50_plot <- young_weight50_plot+
    scale_x_continuous(
      breaks = seq(min(young_weight50$time_numeric), max(young_weight50$time_numeric), by = 7) # breaks every 7 days
    )+
    labs(x = "Days", y = "Body weight (g)", color = NULL)+
    theme(axis.text.x = element_text(size = 7.5),
          axis.title.x = element_text(size = 11))+
    annotate("rect", xmin = -Inf, xmax = 49, # 50 ppm
             fill = "#95BECF", alpha = 0.2,
             ymin = -Inf, ymax = Inf)+
    annotate("rect", xmin = 54, xmax = Inf, # 50 ppm
             fill = "#95BECF", alpha = 0.2,
             ymin = -Inf, ymax = Inf)+
    annotate("rect", xmin = 49, xmax = 54, # DSS
             fill = "gray56", alpha = 0.2,
             ymin = -Inf, ymax = Inf)
  
  young_weight50_plot <- putGgLastLayerBack(young_weight50_plot, nLayers = 3)
  young_weight50_plot <- young_weight50_plot +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
  young_weight50_plot <- young_weight50_plot+
    coord_cartesian(ylim = c(5, 25))+
    geom_signif(xmin = 49,
                xmax = 56,
                annotations = "***",
      y_position = 18.5,
      tip_length = -0.03,
      color = "black",
      size = 0.4,
      textsize = 5,
      hjust = 0.5,
      vjust = 1.6)
  young_weight50_plot
  
  # Stats between t49 and t56
  {
    # Fit the mixed model
    model <- lmer(weight ~ gg_group *timepoint + (1 | id), data = young_weight50)
    
    # Estimated marginal means
    emm <- emmeans(model, ~ timepoint | gg_group)
    
    # Pairwise comparisons between timepoints within each group
    pairs(emm, adjust = "fdr")
    }
  
  # Scatter plot with 500ctrl vs 500dss - all timepoints
  young_weight500_plot <- weightPlot(young_weight500, group = "gg_group", percentage = FALSE,
                                     title = "Iron excess diet (500 ppm)", customColors = c("#F2AA84","#B22222"), customShape = 22)
  young_weight500_plot <- young_weight500_plot+
    scale_x_continuous(
      breaks = seq(min(young_weight500$time_numeric), max(young_weight500$time_numeric), by = 7) # breaks every 7 days
    )+
    labs(x = "Days", y = "Body weigth (g)", color = NULL)+
    theme(axis.text.x = element_text(size = 7.5),
          axis.title.x = element_text(size = 11))+
    annotate("rect", xmin = -Inf, xmax = 35, # 500 ppm
             fill = "#F2AA84", alpha = 0.2,
             ymin = -Inf, ymax = Inf)+
    annotate("rect", xmin = 35, xmax = 49, # 50 ppm
             fill = "#95BECF", alpha = 0.2,
             ymin = -Inf, ymax = Inf)+
    annotate("rect", xmin = 54, xmax = Inf, # 50 ppm
             fill = "#95BECF", alpha = 0.2,
             ymin = -Inf, ymax = Inf)+
    annotate("rect", xmin = 49, xmax = 54, # DSS
             fill = "gray56", alpha = 0.2,
             ymin = -Inf, ymax = Inf)
  
  young_weight500_plot <- putGgLastLayerBack(young_weight500_plot, nLayers = 3)
  young_weight500_plot <- young_weight500_plot +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
  young_weight500_plot <- young_weight500_plot+
    coord_cartesian(ylim = c(5, 25))+
    geom_signif(xmin = 49,
                xmax = 56,
                annotations = "**",
                y_position = 18.5,
                tip_length = -0.03,
                color = "black",
                size = 0.4,
                textsize = 5,
                hjust = 0.5,
                vjust = 1.6)
  young_weight500_plot
  
  # Stats between t49 and t56 for dss group, and between ctrl and dss at last timepoint
  {
    # Fit the mixed model
    model <- lmer(weight ~ gg_group *timepoint + (1 | id), data = young_weight500)
    
    # Estimated marginal means
    emm <- emmeans(model, ~ timepoint | gg_group)
    
    # Pairwise comparisons between timepoints within each group
    pairs(emm, adjust = "fdr")
    
    
    # Estimated marginal means for group, within each timepoint
    emm_groups <- emmeans(model, ~ gg_group | timepoint)
    
    # Pairwise comparisons between groups at each timepoint
    pairs(emm_groups, adjust = "fdr")
    }

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
  dissec_young$gg_group <- factor(dissec_young$gg_group, levels = c("50:water","50:dss","500:water","500:dss"), labels = c("Ctrl:50","DSS:50","Ctrl:500","DSS:500"))
  dissec_young$treatment <- factor(dissec_young$treatment, levels = c("water","dss"))
  
  # Boxplot for body weight
  bw = ironBoxplot(dissec_young, "body_weight", group = "gg_group", title = "At end of experiment", y_axis_title = "Body weight (g)", custom_colors = c("#95BECF","#325BAD","#F2AA84","#B22222"),
                  stats = TRUE, test_results = c("n.s.","n.s.","*","n.s."), text_sizes = c(3,3,5,3),
                  upper_margin = 3, vjustList = c(0,0,0.5,0), tip_length = 0.1, shape = "treatment", scaleFactor = 30)
  bw <- bw+
    guides(x = "axis_nested")+
    theme(plot.title = element_text(size = 10),
          axis.title.y = element_text(size = 10))+
    ylim(15,NA)
  bw
  
  # Statistics
  verifyStatsAssumptions(df = dissec_young, group = "gg_group", measure = "body_weight")
  anova <- aov(body_weight ~ diet * treatment, data = dissec_young) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  
  # Combine plots
  fig2 <- ((young_weight50_plot /young_weight500_plot) | (bw / young_dssflwup_plot))+
    plot_layout(widths = c(2.5, 1))+
    plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 20))
  fig2
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig2.png", plot = fig2, width =10, height =5, dpi = 500)
  
  
  # Lcn-2 DSS mice day53
  df <- as.data.frame(read_excel("Lcn-2 (Thibault).xlsx"))
  colnames(df) <- df[1,]
  df <- df[-1,]
  df$diet <- factor(as.character(df$diet), levels = c("50","500"), labels = c("50 DSS", "500 DSS"))
  df <- df[-18,]
  df$treatment <- "dss"
  
  df$lcn2 <- as.numeric(df$`ng/mg feces`)

  lcn2 <- ironBoxplot(df, "lcn2", group = "diet", title = "Fecal LCN2", y_axis_title = "ng/mg feces", custom_colors = c("#325BAD","#B22222"),
              stats = FALSE, all.ns = TRUE, shape = "treatment")+
    scale_shape_manual(values = c(22))+
    scale_y_continuous(limits = c(0,200))+
    guides(x = "axis_nested")
  
  verifyStatsAssumptions(df, "diet","lcn2")
  wilcox.test(lcn2 ~ diet, data = df)
  
  t.test(lcn2 ~ diet, data = df, var.equal = FALSE)
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/lcn2.png", plot = lcn2, width =2, height =2, dpi = 500, bg = "white")
  
  
}

# Additionnal figure - dss induced colitis model but in the adults
{
  # Set working directory
  setwd("experiments/finished exp/adults-all-exp/")
  
  adult_weight <- as.data.frame(read_xlsx("adult-all/body_weight2.xlsx"))
  # adult_weight <- adult_weight[!(adult_weight$batch == "C" & adult_weight$treatment == "abx"), ] # Remove failed abx mice from batch C
  adult_weight <- adult_weight[!adult_weight$batch == "D" & !adult_weight$treatment == "abx",] # Keep only DSS data and batch A B C
  adult_weight <- adult_weight[-12,] # Remove dead mouse/mice
  adult_weight <- weightDataManipulation(adult_weight, groupInfoCols = 5, fromDay0 = TRUE, arrangeDateColsFormat = FALSE)
  adult_weight$gg_group <- factor(adult_weight$gg_group, levels = c("50 water","50 DSS","500 water","500 DSS"),
                                  labels = c("50 ppm Ctrl","50 ppm DSS","500 ppm Ctrl","500 ppm DSS"))
  adult_weight$timepoint <- factor(adult_weight$time_numeric)
  
  adult_weight50 <- adult_weight[adult_weight$diet == "50",]
  adult_weight500 <- adult_weight[adult_weight$diet == "500",]
  
  # Scatter plot with 50ctrl vs 50dss - all timepoints
  adult_weight50_plot <- weightPlot(adult_weight50, group = "gg_group", percentage = FALSE, title = "Iron sufficient diet (50 ppm)", customColors = c("#95BECF","#325BAD"))
  adult_weight50_plot <- adult_weight50_plot+
    scale_x_continuous(
      breaks = seq(min(adult_weight50$time_numeric), max(adult_weight50$time_numeric), by = 7) # breaks every 7 days
    )+
    labs(x = "Days", y = "Body weight (g)", color = NULL)+
    theme(axis.text.x = element_text(size = 7.5),
          axis.title.x = element_text(size = 11))+
    annotate("rect", xmin = -Inf, xmax = 49, # 50 ppm
             fill = "#95BECF", alpha = 0.2,
             ymin = -Inf, ymax = Inf)+
    annotate("rect", xmin = 54, xmax = Inf, # 50 ppm
             fill = "#95BECF", alpha = 0.2,
             ymin = -Inf, ymax = Inf)+
    annotate("rect", xmin = 49, xmax = 54, # DSS
             fill = "gray56", alpha = 0.2,
             ymin = -Inf, ymax = Inf)
  
  adult_weight50_plot <- putGgLastLayerBack(adult_weight50_plot, nLayers = 3)
  adult_weight50_plot <- adult_weight50_plot +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
  adult_weight50_plot <- adult_weight50_plot+
    coord_cartesian(ylim = c(18, 27.5))+
    geom_signif(xmin = 49,
                xmax = 56,
                annotations = "*",
                y_position = 21,
                tip_length = -0.03,
                color = "black",
                size = 0.4,
                textsize = 5,
                hjust = 0.5,
                vjust = 1.7)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))  # force ~5 ticks
  adult_weight50_plot
  
  # Stats between t49 and t56
  {
    # Fit the mixed model
    model <- lmer(weight ~ gg_group *timepoint + (1 | id), data = adult_weight50)
    
    # Estimated marginal means
    emm <- emmeans(model, ~ timepoint | gg_group)
    
    # Pairwise comparisons between timepoints within each group
    pairs(emm, adjust = "fdr")
  }
  
  # Scatter plot with 500ctrl vs 500dss - all timepoints
  adult_weight500_plot <- weightPlot(adult_weight500, group = "gg_group", percentage = FALSE,
                                     title = "Iron excess diet (500 ppm)", customColors = c("#F2AA84","#B22222"), customShape = 22)
  adult_weight500_plot <- adult_weight500_plot+
    scale_x_continuous(
      breaks = seq(min(adult_weight500$time_numeric), max(adult_weight500$time_numeric), by = 7) # breaks every 7 days
    )+
    labs(x = "Days", y = "Body weigth (g)", color = NULL)+
    theme(axis.text.x = element_text(size = 7.5),
          axis.title.x = element_text(size = 11))+
    annotate("rect", xmin = -Inf, xmax = 35, # 500 ppm
             fill = "#F2AA84", alpha = 0.2,
             ymin = -Inf, ymax = Inf)+
    annotate("rect", xmin = 35, xmax = 49, # 50 ppm
             fill = "#95BECF", alpha = 0.2,
             ymin = -Inf, ymax = Inf)+
    annotate("rect", xmin = 54, xmax = Inf, # 50 ppm
             fill = "#95BECF", alpha = 0.2,
             ymin = -Inf, ymax = Inf)+
    annotate("rect", xmin = 49, xmax = 54, # DSS
             fill = "gray56", alpha = 0.2,
             ymin = -Inf, ymax = Inf)
  
  adult_weight500_plot <- putGgLastLayerBack(adult_weight500_plot, nLayers = 3)
  adult_weight500_plot <- adult_weight500_plot +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
  adult_weight500_plot <- adult_weight500_plot+
    coord_cartesian(ylim = c(18, 27.5))+
    geom_signif(xmin = 49,
                xmax = 56,
                annotations = "**",
                y_position = 20.5,
                tip_length = -0.03,
                color = "black",
                size = 0.4,
                textsize = 5,
                hjust = 0.5,
                vjust = 1.8)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))  # force ~5 ticks
  adult_weight500_plot
  
  # Stats between t49 and t56 for dss group, and between ctrl and dss at last timepoint
  {
    # Fit the mixed model
    model <- lmer(weight ~ gg_group *timepoint + (1 | id), data = adult_weight500)
    
    # Estimated marginal means
    emm <- emmeans(model, ~ timepoint | gg_group)
    
    # Pairwise comparisons between timepoints within each group
    pairs(emm, adjust = "fdr")
    
    
    # Estimated marginal means for group, within each timepoint
    emm_groups <- emmeans(model, ~ gg_group | timepoint)
    
    # Pairwise comparisons between groups at each timepoint
    pairs(emm_groups, adjust = "fdr")
  }
  
  # Mice DAI evolution 
  adult_dss_followup <- read.csv("combined_adult_dss_followup.csv", header = TRUE, sep=";")
  adult_dss_followup <- dssFollowupManipulation(df = adult_dss_followup,groupInfoCols = 4,dateStart = "2024-05-29",nbrDays = 5, negativeOnly = TRUE)
  adult_dssflwup_plot <- dssDiseaseIndexPlot(adult_dss_followup, statBarLast = TRUE, signifAnnotation = "n.s.", signifPositionShift = 0.3, signifSize = 3)
  adult_dssflwup_plot <- adult_dssflwup_plot+
    theme(axis.ticks.x = element_line(),
          legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5))
  adult_dssflwup_plot
  
  # Statistics for last day of DAI
  verifyStatsAssumptions(df = adult_dss_followup[adult_dss_followup$time_numeric == "5",], group = "gg_group", measure = "index")
  wilcox.test(index ~ gg_group, data = adult_dss_followup[adult_dss_followup$time_numeric == "5",])
  
  # Mice dissection data
  dissec_adult <- read_excel("all-adults-dissection.xlsx")
  dissec_adult <- dissec_adult[!dissec_adult$treatment == "abx",] # Keep only DSS data and batch A B C
  dissec_adult <- dissec_adult[-12,] # Remove dead mouse/mice
  dissec_adult$gg_group <- interaction(dissec_adult$diet, dissec_adult$treatment, sep = ":")
  dissec_adult$gg_group <- factor(dissec_adult$gg_group, levels = c("50:water","50:dss","500:water","500:dss"), labels = c("Ctrl:50","DSS:50","Ctrl:500","DSS:500"))
  dissec_adult$treatment <- factor(dissec_adult$treatment, levels = c("water","dss"))
  dissec_adult$diet <- factor(dissec_adult$diet, levels = c("50","500"))
  colnames(dissec_adult)[6] <- "body_weight"
  
  # Boxplot for body weight
  bw = ironBoxplot(dissec_adult, "body_weight", group = "gg_group", title = "At end of experiment", y_axis_title = "Body weight (g)", custom_colors = c("#95BECF","#325BAD","#F2AA84","#B22222"),
                   stats = FALSE, test_results = c("n.s.","n.s.","n.s.","n.s."), text_sizes = c(3,3,5,3),
                   upper_margin = 3, vjustList = c(0,0,0.5,0), tip_length = 0.1, shape = "treatment", scaleFactor = 30)
  bw <- bw+
    guides(x = "axis_nested")+
    theme(plot.title = element_text(size = 10),
          axis.title.y = element_text(size = 10))+
    ylim(15,NA)
  bw
  
  # Statistics
  verifyStatsAssumptions(df = dissec_adult, group = "gg_group", measure = "body_weight")
  anova <- aov(body_weight ~ diet * treatment, data = dissec_adult) #Fit a ANOVA model
  summary(anova)
  results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
  results
  
  pairwise.wilcox.test(dissec_adult$body_weight, dissec_adult$gg_group, p.adjust.method = "fdr")
  
  # Combine plots
  fig2_adults <- ((adult_weight50_plot /adult_weight500_plot) | (bw / adult_dssflwup_plot))+
    plot_layout(widths = c(2.5, 1))+
    plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 20))
  fig2_adults
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fig2_adults.png", plot = fig2_adults, width =10, height =5, dpi = 500)
}

# PCR validations
{
  # Set working directory
  setwd("experiments/finished exp/young-DSS-exp3")
  
  # PCR for L. murinus and M. intestinale and F. rodentium
  df <- read_excel("pcr_validation/F rodentium RT-PCR YoungDSS Exp.xlsx")
  colnames(df) <- df[3,]
  df <- df[!is.na(df$`Muribaculum intestinale`),]
  df <- df[-1,]

  # Metadata handling
  meta <- read.csv("young48_dss_dissection.csv", sep = ";")
  meta$id <- substring(meta$id, 1, 5)
  
  # Bind metadata information to dataframe
  df <- merge(df, meta, by.x = "Mouse ID", by.y = "id")
  # df$lm <- log10(as.numeric(df$`Lm/16S`))
  # df$mi <- log10(as.numeric(df$`Mi/16S`))
  
  df$lm <- as.numeric(df$`Lm/16S`)
  df$mi <- as.numeric(df$`Mi/16S`)
  df$fr <- as.numeric(df$`Fr/16S`)
  df$diet <- factor(as.character(df$diet), levels = c("50","500"), labels = c("50 DSS", "500 DSS"))
  df <- df[,-c(5,10,15)]
  
  lm <- ironBoxplot(df, "lm", group = "diet", title = "L. murinus", y_axis_title = "L.murinus/16S", custom_colors = c("#325BAD","#B22222"),
              stats = FALSE, all.ns = FALSE, shape = "treatment")+
    scale_shape_manual(values = c(22))+
    scale_y_continuous(limits = c(0,0.06), n.breaks = 5)+
    guides(x = "axis_nested")
  lm
  
  verifyStatsAssumptions(df, "diet","lm")
  t.test(lm ~ diet, data = df, var.equal = TRUE)
  wilcox.test(lm ~ diet, data = df)
  
  mi <- ironBoxplot(df, "mi", group = "diet", title = "M. intestinale", y_axis_title = "M. intestinale/16S", custom_colors = c("#325BAD","#B22222"),
                      stats = FALSE, all.ns = TRUE, shape = "treatment")+
    scale_shape_manual(values = c(22))+
    scale_y_continuous(limits = c(0,0.6), n.breaks = 4)+
    guides(x = "axis_nested")
  mi
  
  
  verifyStatsAssumptions(df, "diet","mi")
  t.test(mi ~ diet, data = df, var.equal = TRUE)
  t.test(mi ~ diet, data = df, var.equal = FALSE)
  wilcox.test(mi ~ diet, data = df)
  
  df$log_mi <- log10(df$mi)
  verifyStatsAssumptions(df, "diet","log_mi")
  t.test(log_mi ~ diet, data = df, var.equal = FALSE)
  wilcox.test(log_mi ~ diet, data = df)
  
  fr <- ironBoxplot(df, "fr", group = "diet", title = "F. rodentium", y_axis_title = "F. rodentium/16S", custom_colors = c("#325BAD","#B22222"),
                    stats = FALSE, all.ns = TRUE, shape = "treatment")+
    scale_shape_manual(values = c(22))+
    scale_y_continuous(limits = c(0,NA), n.breaks = 4)+
    guides(x = "axis_nested")
  fr

  
  verifyStatsAssumptions(df, "diet","fr")
  t.test(fr ~ diet, data = df, var.equal = TRUE)
  wilcox.test(fr ~ diet, data = df)
  
  
  pcrs <- (fr | mi | lm)+
    plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 20))
  pcrs
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/pcrs.png", plot = pcrs, width =6, height =2.5, dpi = 500, bg = "white")
  
  
  # Correlation between the species qpcrs abundance
  df_cor <- df[,c(1,18,24,26)]
  row.names(df_cor) <- df_cor[,1] 
  df_cor <- df_cor[,-1]
  df_cor$diet <- factor(df_cor$diet, labels = c("50 ppm DSS","500 ppm DSS"))
  
  cor_results <- rcorr(as.matrix(df_cor[,-1]), type = "spearman")
  cor_matrix <- cor_results$r  # Extract correlation coefficients
  p_values <- cor_results$P    # Extract p-values
  r = -0.6116601
  
  cor <- ggplot(df_cor, mapping = aes(x = lm, y = fr, fill = diet))+
    geom_smooth(aes(x = lm, y = fr, group = 1),
                method = "lm", se = TRUE, color = "black", inherit.aes = FALSE) +
    geom_point(shape = 22, color = "black")+
    labs(x= "L. murinus/16S", y = "F. rodentium/16S", fill = NULL)+
    scale_fill_manual(values = c("#325BAD","#B22222"))+
    scale_x_continuous(n.breaks = 4, limits = c(0,0.06))+
    scale_y_continuous(n.breaks = 4)+
    my_theme()
  
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/fr_lm_cor.png", plot = cor, width =4, height =3, dpi = 500, bg = "white")
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
source("~/CHUM_git/gut-microbiota-iron/other scripts/dataManipFunctions.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/alpha_diversity_graphs_and_stats.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/beta_diversity_graphs_and_stats.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/correlation_graphs_and_stats.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/relab_analysis_graphs_and_stats.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/taxa_distrib_graphs_and_stats.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/plot_microbiota_extension.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/deseq2_log2fold_change_analysis.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/chronobiome.R")
source("~/CHUM_git/gut-microbiota-iron/pipeline/microbiota_analysis/picrust2_graphs.R")

# For microbiota 18 - dss experiment
{
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

# Figure 3 - alpha diversity
{
  
  # Alpha diveristy timeline
  {
    # Week as numeric
    sample_data(ps)$timepoint  <- as.numeric(as.character(sample_data(ps)$timepoint))
    sample_data(ps)$gg_group2 
    sample_data(ps)$gg_group2 <- factor(sample_data(ps)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
    sample_data(ps)$treatment
    
    #Estinate richness measures for dataset
    richness_data <- estimate_richness(ps, measures = c("Chao1", "Shannon", "InvSimpson"))
    alpha_d <- cbind(as.data.frame(sample_data(ps)), richness_data)
    
    custom_colors <- c("#95BECF","#F2AA84","#325BAD","#B22222")
    graphs <- alphaDiversityTimeline(ps, time = "timepoint", group = "gg_group2", shape = "treatment", custom_colors, semRibbons = TRUE)
    
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
      ylim(0,NA)
    p1
    
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
        breaks = seq(min(sample_data(ps)$timepoint), max(sample_data(ps)$timepoint), by = 7) # breaks every 7 days                    # "T" before each label
      ) +
      labs(y = "Index", x = "Days", color = NULL, fill = "", title = 'Shannon')+
      guides(shape = "none", fill = "none")+
      theme(axis.ticks.x = element_line(),
            legend.background = element_rect(color = "black", fill  = "white", linewidth = 0.5),
            axis.text.x = element_text(size = 6)
      )+
      ylim(0,NA)
    p2
    
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
      ylim(0,NA)
    p3
    
    alpha_d_plot <- p1 / p2 /p3 +plot_layout(guides = "collect")
    ggsave(filename = "~/CHUM_git/figures/memoire/dss/alphad.png", plot = alpha_d_plot, width = 5, height = 5, dpi = 500)
  }
  
  # Alpha diveristy timeline v2 (ctrl vs dss for each diet) + boxplots at last timepoint
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
    alpha_d_v2 <- (alpha_d50 | alpha_d500 | tfinal_stack)+
      plot_layout(widths = c(2, 2, 1))
    alpha_d_v2
    
    ggsave(filename = "~/CHUM_git/figures/memoire/dss/alphad_v2.png", plot = alpha_d_v2, width = 12, height = 6, dpi = 500)
  }
  
}

# Verifying that there is no issue regarding alpha diversity calculation
{
  richness_data <- estimate_richness(ps_t54, measures = c("Chao1", "Shannon", "InvSimpson"))
  rownames(richness_data) <- substring(rownames(richness_data), 2 , 10)
  df <- merge(richness_data, sample_data(ps_t54), by = "row.names")
  mean(df$Shannon[df$gg_group2 == "50:dss"])
  mean(df$Shannon[df$gg_group2 == "50:water"])
  
  mean(df$Shannon[df$gg_group2 == "500:dss"])
  mean(df$Shannon[df$gg_group2 == "500:water"])
  
  
  library(vegan)
}

# Figure 4 - Beta diversity
{
  # Beta diversity
  {
    set.seed(100)
    
    phy_tree(ps_flt_diet) <- midpoint(tree) # Add rooted tree to phyloseq object
    
    # For diet timepoints
    sample_data(ps_flt_diet)$diet <- factor(sample_data(ps_flt_diet)$diet, labels = c("50 ppm","500 ppm"))
    
    # Weighted unifrac
    diet_beta_d_plot <- betaDiversityTimepoint2Factors(ps_flt_diet, sample_id = "sample_id", timeVariable = "timepoint",
                                                       varToCompare =  "diet", shape = "treatment", distMethod ="wunifrac",
                                                       transform = "rel_ab", customColors = c("#95BECF","#F2AA84"),
                                                       font = "Arial", path = "../figures/Thibault_dss/beta_diversity/diet/", 
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
                                                           font = "Arial", path = "../figures/Thibault_dss/beta_diversity/diet_dss/",
                                                           additionnalAes = my_theme()+theme(plot.title = element_text(size = 10)), dim = c(4,8), combineGraphs = TRUE, returnFig = TRUE,
                                                           customTitles = c("End of DSS treatment","End of experiment"), pairwiseAdonis = FALSE, displayPValue = FALSE, hideLegend = TRUE)
    
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
    
    # Compute pvalues of pairwise comparisons with FDR correction
    {
      # For diet + treatment at t54 and tfinal
      phy_tree(ps_flt_dss) <- midpoint(tree) # Add rooted tree to phyloseq object
      sample_data(ps_flt_dss)$gg_group2 <- factor(sample_data(ps_flt_dss)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl", "50 ppm DSS", "500 ppm DSS"))
      # Weighted unifrac => not very informative
      betaDiversityTimepoint2Factors(ps_flt_dss, sample_id = "sample_id", timeVariable = "timepoint",
                                     varToCompare =  "gg_group2", distMethod ="wunifrac", shape = "treatment",
                                     customColors = c("#95BECF","#F2AA84","#325BAD","#B22222"),
                                     font = "Arial", path = "../figures/Thibault_dss/beta_diversity/diet_dss/",
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
                                                              font = "Arial", path = "../figures/Thibault_dss/beta_diversity/dss_only/",
                                                              additionnalAes = theme_beta_d(), dim = c(2.5,4), returnFig = TRUE, displayPValue = FALSE,
                                                              customTitles = c("50 ppm Ctrl\nvs 500 ppm Ctrl"), combineGraphs = FALSE, hideLegend = TRUE, pvalueToDisplay = 0.268, subGraph = TRUE,
                                                              size = 1, stroke = 0.15)
    
    phy_tree(ps_tfinal_flt) <- midpoint(tree) 
    Ctrl50vs500_beta_d_plot <- betaDiversityTimepoint2Factors(ps_tfinal_flt , sample_id = "sample_id", timeVariable = "timepoint",
                                                              varToCompare =  "gg_group2", distMethod ="wunifrac", shape = "treatment",
                                                              customColors = c("#95BECF","#F2AA84"),
                                                              font = "Arial", path = "../figures/Thibault_dss/beta_diversity/dss_only/",
                                                              additionnalAes = theme_beta_d(), dim = c(2.5,4), returnFig = TRUE, displayPValue = FALSE,
                                                              customTitles = c("50 ppm Ctrl\nvs 500 ppm Ctrl"), combineGraphs = FALSE, hideLegend = TRUE, pvalueToDisplay = 0.268, subGraph = TRUE,
                                                              size = 1, stroke = 0.15, selectedPairwise = c("50:water","500:water"))
    
    
    Ctrl50vs500_beta_d_plot
    
    
    # Weighted unifrac for dss groups - last timepoint
    ps_sub <- prune_samples(sample_data(ps_tfinal_flt)$treatment == "dss", ps_tfinal_flt)
    sample_data(ps_sub)$diet <- factor(sample_data(ps_sub)$diet, labels = c("50 ppm DSS","500 ppm DSS"))
    phy_tree(ps_sub) <- midpoint(tree) # Add rooted tree to phyloseq object
    DSS_50v500_beta_d_plot <- betaDiversityTimepoint2Factors(ps_sub , sample_id = "sample_id", timeVariable = "timepoint",
                                                             varToCompare =  "diet", distMethod ="wunifrac", shape = "treatment",
                                                             customColors = c("#325BAD","#B22222"),
                                                             font = "Arial", path = "../figures/Thibault_dss/beta_diversity/dss_only/",
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
                                                               font = "Arial", path = "../figures/Thibault_dss/beta_diversity/dss_only/",
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
                                                                 font = "Arial", path = "../figures/Thibault_dss/beta_diversity/dss_only/",
                                                                 additionnalAes =theme_beta_d(), dim = c(2.5,4), returnFig = TRUE, displayPValue = FALSE,
                                                                 customTitles = c("500 ppm Ctrl\nvs 500 ppm DSS"), combineGraphs = FALSE, hideLegend = TRUE, pvalueToDisplay = 0.048, subGraph = TRUE, positionPvalue = "right",
                                                                 size = 1, stroke = 0.15)
    Ctrl500vDSS500_beta_d_plot
    
    
  }
  
  # Combine all graphs
  stack <- Ctrl50vs500_beta_d_plot+DSS_50v500_beta_d_plot + Ctrl50vDSS50_beta_d_plot+Ctrl500vDSS500_beta_d_plot &
    theme(plot.margin = margin(3, 8, 3, 8))
  main <- (diet_beta_d_plot | diet_dss_beta_d_plot)+
    plot_layout(widths = c(1.5,1))
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/beta_main.png", plot = main, width = 12, height = 2.5, dpi = 500)
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/beta_sub.png", plot = stack, width = 2.6, height = 2.5, dpi = 500)
}

# Stackbar extended at last day of DSS - before figure 5?
{
  # First, timepoints and groups must be ordered properly and as factors
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
  
  library(ggh4x)
  
  facet_scales <- list(
    scale_x_discrete(labels = as.character(1:12)),
    scale_x_discrete(labels = as.character(13:24)),
    scale_x_discrete(labels = as.character(25:35)),
    scale_x_discrete(labels = as.character(36:45))
  )
  
  # Custom the plot
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
  
  # Saving the plot and the associated stats
  existingDirCheck("../figures/Thibault_dss/stackbar")
  writeStackbarExtendedSigTable(main_table = diet_dss_phyla_fam$significant_table_main, includeSubTable = TRUE, sub_table = diet_dss_phyla_fam$significant_table_sub, filepath = "../figures/Thibault_dss/stackbar/diet_dss_stackbar_stats.xlsx")
  
  # pvalues heatmap for the main lvl stats
  pvalHmapPhyla <- pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_dss/stackbar/diet_dss_stackbar_stats.xlsx")),
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
  pvalHmapFamily <- pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_dss/stackbar/diet_dss_stackbar_stats.xlsx")),
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
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/stackbar.png", plot = stackbar_t54, width = 7, height = 7, dpi = 500)
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/stackbar_hmap.png", plot = statHmap, width = 3.3, height = 9, dpi = 500)
  
  
}

# Figure 5 - Chronobiome
{
  
  # Chronobiome 
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
      labs(x = "Days")
    
    
    chronobiome <- chronobiome+
      coord_cartesian(xlim = c(0, 112), ylim = c(0, 100), expand = FALSE)+
      theme(panel.spacing = unit(0.75, "cm"))
  }

  ggsave(filename = "~/CHUM_git/figures/memoire/dss/chronobiome.png", plot = chronobiome, width = 15, height = 12, dpi = 500)
  
  # Stats at last timepoint - Family
  ps_taxa <- tax_glom(ps_tfinal_flt, taxrank = "Family", NArm = FALSE)
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
  sig_df_list_family <- sig_df_list
  
  # Stats at last timepoint - Phylum
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
  
  # Saving the plot and the associated stats
  existingDirCheck("../figures/Thibault_dss/stackbar")
  writeStackbarExtendedSigTable(main_table =  sig_df_list_phylum, includeSubTable = TRUE, sub_table =  sig_df_list_family, filepath = "../figures/Thibault_dss/stackbar/final_stackbar_stats.xlsx")
  
  # pvalues heatmap for the main lvl stats
  pvalHmapPhyla <- pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_dss/stackbar/final_stackbar_stats.xlsx")),
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
  pvalHmapFamily <- pvaluesHmap(stats = as.data.frame(readxl::read_excel("../figures/Thibault_dss/stackbar/final_stackbar_stats.xlsx")),
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
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/final_timepoint_stats_hmap.png", plot = statHmap, width = 3.3, height = 8.7, dpi = 500)
  
  
  # Stats timepoint t0 t35 and t49 between 50 ppm and 500 ppm diets, for family level 
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
  
  
  # Chronobiome diet
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
        )
    }
    
    sample_data(ps_flt_diet)$diet <- factor(sample_data(ps_flt_diet)$diet, labels = c("50 ppm","500 ppm"))
    
    diet_chronobiome <- plot_timeline_2_groups(
      ps_object = ps_flt_diet,
      exp_group =  "diet", # must be as factor
      time_group = "timepoint", # must be as factor
      sample_name = "sample_id",
      main_level = 'Phylum',
      sub_level = 'Family',
      average_relab_per_group = TRUE,
      smoothing = FALSE,
      n_phy = 4,
      hues = c("Blues", "Greens", "Purples", "Oranges"),
      color_bias = 2,
      custom_theme = theme_chronobiome(),
      stats = TRUE
    )
    
    diet_chronobiome$main_fig
    diet_chronobiome$sub_stat_plot
    diet_chronobiome$main_stat_plot
    
    main_fig <-  diet_chronobiome$main_fig+
      facet_wrap2(~ diet, 
                  scales  = "free_x", nrow = 2, ncol = 1,
                  strip = strip_themed(background_x = elem_list_rect(fill = c("#95BECF","#F2AA84"))))+
      scale_x_continuous(breaks = seq(min(as.numeric(levels(sample_data(ps_flt_diet)$timepoint))), max(as.numeric(levels(sample_data(ps_flt_diet)$timepoint))), by = 7))+
      labs(x = "Days")
    
    main_fig
    
    existingDirCheck("~/CHUM_git/figures/memoire/dss/newVisuals")
    ggsave(filename = "~/CHUM_git/figures/memoire/dss/newVisuals/chronobiome.png", plot = main_fig, width = 8, height = 8, dpi = 500)
    ggsave(filename = "~/CHUM_git/figures/memoire/dss/newVisuals/main_stats.png", plot =  diet_chronobiome$main_stat_plot, width = 4, height = 2, dpi = 500, bg = "white")
    ggsave(filename = "~/CHUM_git/figures/memoire/dss/newVisuals/sub_stats.png", plot =  diet_chronobiome$sub_stat_plot, width = 4, height = 5.5, dpi = 500, bg = "white")
  }
  
  # Function that produces a stat table for comparison between 2 groups with DESeq2 for different timepoints
  # Timepoints and group variables must be ordered as factors
  timelineStatTable <- function(ps, groupVar, timeVar, taxonLevel = "Family",
                                top_n = 30, p_col = "padj", p_sig = 0.05,
                                title = NULL) {
    
    if(taxonLevel != "Species"){
      ps <- tax_glom(ps, taxrank = taxonLevel)
    }
    
    res_all <- data.frame() # to store all results
    
    for(timePoint in levels(sample_data(ps)[[timeVar]])){
      
      # Subset to this timepoint
      ps_subset <- prune_samples(sample_data(ps)[[timeVar]] == timePoint, ps)
      ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)
      
      # DESeq2 analysis
      deseq <- phyloseq_to_deseq2(ps_subset, as.formula(paste("~", groupVar)))
      deseq <- DESeq(deseq, test = "Wald", fitType = "parametric", quiet = TRUE)
      res <- results(deseq)
      
      # Convert to data frame and add taxonomy
      sigtab <- cbind(as.data.frame(res),
                      as(tax_table(ps_subset)[rownames(res), ], "matrix"))
      sigtab$padj[is.na(sigtab$padj)] <- 1
      sigtab[[timeVar]] <- timePoint
      
      # rbind to the main table
      res_all <- rbind(res_all, sigtab)
    }
    
    dat <- within(res_all, {
      taxon <- ifelse(is.na(get(taxonLevel)) | get(taxonLevel) == "",
                      "Unclassified", get(taxonLevel))
      tp    <- get(timeVar)
      coef  <- as.numeric(log2FoldChange)
      pval  <- as.numeric(get(p_col))
      pval[is.na(pval)] <- 1
      sig   <- ifelse(pval < p_sig, "*", "")
      sizevar <- abs(coef)
      p_for_col <- pmin(pmax(pval, 1e-16), 1)
    })
    
    View(dat)
    return(NULL)
    
    # choose top taxa to keep the y-axis readable
    keep <- aggregate(sizevar ~ taxon, dat, function(x) mean(x, na.rm = TRUE))
    keep <- head(keep[order(-keep$sizevar), "taxon"], top_n)
    dat  <- dat[dat$taxon %in% keep, , drop = FALSE]
    
    # order axes
    ord  <- aggregate(sizevar ~ taxon, dat, mean, na.rm = TRUE)
    dat$taxon <- factor(dat$taxon, levels = rev(ord[order(ord$sizevar), "taxon"]))
    if (is.numeric(dat$tp)) dat$tp <- factor(dat$tp, levels = sort(unique(dat$tp)))
    
    p <- ggplot(dat, aes(x = tp, y = taxon)) +
      geom_point(aes(size = sizevar, fill = p_for_col),
                 shape = 21, stroke = 0.2, color = "grey20", alpha = 0.9) +
      geom_text(aes(label = sig), size = 4, color = "white") +
      scale_size_continuous(name = "Coefficient", range = c(1.5, 10)) +
      scale_fill_gradientn(
        name = "P.Value",
        colours = c("#b2182b","#fddbc7","#f7f7f7","#d1e5f0","#2166ac"),
        values = rescale(c(0.001, 0.01, 0.05, 0.25, 0.75)),
        limits = c(0, 0.75), oob = squish
      ) +
      labs(x = NULL, y = NULL, title = title) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5))
    
    return(list(table = res_all, plot = p))
  }
  
  out <- timelineStatTable(ps_flt_diet, "diet", "timepoint", "Family",
                           top_n = 30, p_col = "padj",
                           title = "Diet effect over time")
  
  sample_data(ps_flt_diet)$timepoint
  sample_data(ps_flt_diet)$diet
  
  timelineStatTable(ps_flt_diet, "diet","timepoint")
  
}

# Old figure 5 - bacteria relative abundance and overall composition of the gut microbiota community
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
                  strip = strip_themed(background_x = elem_list_rect(fill = c("#95BECF","#325BAD","#F2AA84","#B22222"))), 
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
                  strip = strip_themed(background_x = elem_list_rect(fill = c("#95BECF","#F2AA84","#325BAD","#B22222"))))+
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

# Figure 6 - differentially abundant bacterial species at last timepoint
{
  # Volcano plot - 50 DSS VS 500 DSS last timepoint
  {
      deseq <- phyloseq_to_deseq2(ps_tfinal_flt, ~ treatment+diet+diet:treatment)
      deseq <- DESeq(deseq, test="Wald", fitType = "parametric")
      print(resultsNames(deseq))
      
      # Manually add species names that we identified with BLASTn (ASV 8,38,19,6 had no species assignments) 
      # To retrieve sequences to compare theme to blast, run 
      # "as.character(c(refseq(ps_tfinal_flt)["ASV8"],
      # refseq(ps_tfinal_flt)["ASV38"],
      # refseq(ps_tfinal_flt)["ASV19"],
      # refseq(ps_tfinal_flt)["ASV6"]))"
      
      tax_table(ps_tfinal_flt)["ASV8","Species"] <- "murinus (ASV8)"
      tax_table(ps_tfinal_flt)["ASV38","Species"] <- "murinus (ASV38)"
      tax_table(ps_tfinal_flt)["ASV19","Species"] <- "intestinale"
      tax_table(ps_tfinal_flt)["ASV6","Species"] <- "rodentium (ASV6)"
      tax_table(ps_tfinal_flt)["ASV1","Species"] <- "rodentium (ASV1)"

      #For a given taxononical levels, creates graph for each timepoint, displaying a volcano plot for taxononical level of interest
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
          )
      volcanoPlot50vs500dss
      
      volcanoPlot50vs500dss <- volcanoPlot50vs500dss+labs(color = "Significance")+
        scale_color_manual(
          breaks = c("n.s.", "Unknown species", "Up", "Down"),
          values = c('grey','black','#B22222','#325BAD'))
    
  }
  
  # Individual graphs for each species displayed as full timeline
  {
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
      axis.title.y   = element_text(size = 8, margin = margin(r = 10)),
      plot.title     = element_text(size = 10, margin = margin(b = 10)),
      plot.margin    = margin(5, 5, 5, 5),    # optional: equalize outer margins
      legend.position = "none"
      # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
    
    # Species M. intestinale ASV19
    asv19 <- asvRelAbDistributionTimeline(ps_flt_all, "ASV19", taxon = "Species",
                                 "gg_group2", "timepoint", c("#95BECF","#F2AA84","#325BAD","#B22222"),
                                 displayASVNumber = FALSE)+
      labs(color = NULL, x = "Days")+
      scale_x_continuous(
        breaks = seq(min(sample_data(ps_flt_all)$timepoint), max(sample_data(ps_flt_all)$timepoint), by = 7)                   # "T" before each label
      )+
      timeline_theme+
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
    
    asv19 <- putGgLastLayerBack(asv19, nLayers = 4)
    asv19

    
    # Species ligilactobacillus murinus ASV 8
    asv8 <- asvRelAbDistributionTimeline(ps_flt_all, "ASV8", taxon = "Species",
                                         "gg_group2", "timepoint", c("#95BECF","#F2AA84","#325BAD","#B22222"),
                                         displayASVNumber = FALSE)+
      labs(color = NULL, x = "Days")+
      scale_x_continuous(
        breaks = seq(min(sample_data(ps_flt_all)$timepoint), max(sample_data(ps_flt_all)$timepoint), by = 7)                   # "T" before each label
      )+
      timeline_theme+
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
    
    asv8 <- putGgLastLayerBack(asv8, nLayers = 4)
    asv8
    
    # Species ligilactobacillus murinus ASV 38
    asv38 <- asvRelAbDistributionTimeline(ps_flt_all, "ASV38", taxon = "Species",
                                 "gg_group2", "timepoint", c("#95BECF","#F2AA84","#325BAD","#B22222"),
                                 displayASVNumber = FALSE)+
      labs(color = NULL, x = "Days")+
      scale_x_continuous(
        breaks = seq(min(sample_data(ps_flt_all)$timepoint), max(sample_data(ps_flt_all)$timepoint), by = 7)                   # "T" before each label
      )+
      timeline_theme+
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
    
    asv38 <- putGgLastLayerBack(asv38, nLayers = 4)
    asv38
    
    # F rodentium - ASV1
    asv1 <- asvRelAbDistributionTimeline(ps_flt_all, "ASV1", taxon = "Species",
                                        "gg_group2", "timepoint", c("#95BECF","#F2AA84","#325BAD","#B22222"),
                                        displayASVNumber = FALSE)+
      labs(color = NULL, x = "Days")+
      scale_x_continuous(
        breaks = seq(min(sample_data(ps_flt_all)$timepoint), max(sample_data(ps_flt_all)$timepoint), by = 7)                   # "T" before each label
      )+timeline_theme+
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
    
    asv1 <- putGgLastLayerBack(asv1, nLayers = 4)
    asv1
    
    # F rodentium - ASV6
    asv6 <- asvRelAbDistributionTimeline(ps_flt_all, "ASV6", taxon = "Species",
                                         "gg_group2", "timepoint", c("#95BECF","#F2AA84","#325BAD","#B22222"),
                                         displayASVNumber = FALSE)+
      labs(color = NULL, x = "Days")+
      scale_x_continuous(
        breaks = seq(min(sample_data(ps_flt_all)$timepoint), max(sample_data(ps_flt_all)$timepoint), by = 7)                   # "T" before each label
      )+timeline_theme+
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
    
    asv6 <- putGgLastLayerBack(asv6, nLayers = 4)
    asv6
    
    
    # Combine individual graphs into one with merged legend
    species_relab_timeline <- (asv1 / asv6 / asv19 / asv8 / asv38)+
      plot_layout(guides = "collect")
    species_relab_timeline
      
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
    
    timeline_theme <- theme(
      axis.text.x    = element_text(size = 7),
      axis.title.x   = element_text(size = 10, margin = margin(t = 10)),
      axis.title.y   = element_text(size = 8, margin = margin(r = 10)),
      plot.title     = element_text(size = 10, margin = margin(b = 10)),
      plot.margin    = margin(5, 5, 5, 5),    # optional: equalize outer margins
      legend.position = "none"
      # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
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
  
  # Same as above but with breaks within the graph
  {
    sample_data(ps_flt_all)$gg_group2 <- factor(sample_data(ps_flt_all)$gg_group2, labels = c("50 ppm Ctrl","500 ppm Ctrl","50 ppm DSS","500 ppm DSS"))
    sample_data(ps_flt_all)$timepoint <- as.numeric(as.character(sample_data(ps_flt_all)$timepoint))
    sample_data(ps_flt_all)$period <- ifelse(sample_data(ps_flt_all)$timepoint < 35, "Early", "Late")

    # Manually add species names that we identified with BLASTn
    tax_table(ps_flt_all)["ASV8","Species"] <- "murinus (ASV8)"
    tax_table(ps_flt_all)["ASV38","Species"] <- "murinus (ASV38)"
    tax_table(ps_flt_all)["ASV19","Species"] <- "intestinale"
    tax_table(ps_flt_all)["ASV6","Species"] <- "rodentium (ASV6)"
    tax_table(ps_flt_all)["ASV1","Species"] <- "rodentium (ASV1)"
    
    timeline_theme <- theme(
      axis.text.x    = element_text(size = 7),
      axis.title.x   = element_text(size = 10, margin = margin(t = 10)),
      axis.title.y   = element_text(size = 8, margin = margin(r = 10)),
      plot.title     = element_text(size = 10, margin = margin(b = 10)),
      plot.margin    = margin(5, 5, 5, 5),    # optional: equalize outer margins
      legend.position = "none"
      # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
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
                                                   diet = "50")+
      facet_wrap(~period, scales = "free_y", ncol = 1)+
      theme(strip.text = element_blank())+
      scale_y_continuous(n.breaks = 3)
    asv19
    
    # Species ligilactobacillus murinus ASV 8
    asv8 <- asvRelAbDistributionTimelineExtension(ps_sub50, "ASV8", taxon = "Species",
                                                  "gg_group2", "timepoint", c("#95BECF","#325BAD"),
                                                  diet = "50")+
      facet_wrap(~period, scales = "free_y", ncol = 1)+
      theme(strip.text = element_blank())+
      scale_y_continuous(n.breaks = 3)
    asv8
    
    # Species ligilactobacillus murinus ASV 38
    asv38 <- asvRelAbDistributionTimelineExtension(ps_sub50, "ASV38", taxon = "Species",
                                                   "gg_group2", "timepoint", c("#95BECF","#325BAD"),
                                                   diet = "50")+
      facet_wrap(~period, scales = "free_y", ncol = 1)+
      theme(strip.text = element_blank())+
      scale_y_continuous(n.breaks = 3)
    asv38
    
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
                                                   diet = "500")+
      facet_wrap(~period, scales = "free_y", ncol = 1)+
      theme(strip.text = element_blank())+
      scale_y_continuous(n.breaks = 3)
    asv19
    
    # Species ligilactobacillus murinus ASV 8
    asv8 <- asvRelAbDistributionTimelineExtension(ps_sub500, "ASV8", taxon = "Species",
                                                  "gg_group2", "timepoint", c("#F2AA84","#B22222"),
                                                  diet = "500")+
      facet_wrap(~period, scales = "free_y", ncol = 1)+
      theme(strip.text = element_blank())+
      scale_y_continuous(n.breaks = 3)
    
    # Species ligilactobacillus murinus ASV 38
    asv38 <- asvRelAbDistributionTimelineExtension(ps_sub500, "ASV38", taxon = "Species",
                                                   "gg_group2", "timepoint", c("#F2AA84","#B22222"),
                                                   diet = "500")+
      facet_wrap(~period, scales = "free_y", ncol = 1)+
      theme(strip.text = element_blank())+
      scale_y_continuous(n.breaks = 3)
    
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
  
  # Create figure 5
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/volcano.png", plot = volcanoPlot50vs500dss, width = 6, height = 5, dpi = 500)

  rel_ab <- (wrap_elements(full = species_relab_timeline50)+ wrap_elements(full = species_relab_timeline500) + wrap_elements(full = species_relab_time_final))+
    plot_layout(widths = c(2,2,1))
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/species_relab.png", plot = rel_ab, width = 10.5, height = 8, dpi = 500)
  
}

# Identifying differentially abundant taxa at the different timepoints
{
  
  # Relative abundance analysis: finding differential abundant bugs at the species level, for diet groups only
  {
    # Path where to save graphs
    pathToSave <- "~/CHUM_git/figures/Thibault_dss/relative_abundance_diet/"
    existingDirCheck(pathToSave)
    
    #customColors for graph display
    customColors = c("#95BECF","#F2AA84")
    customPhylaColors = c("#e6550d","#31a354", "#583093","skyblue3")
    
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
      relabSingleTimepoint(ps_subset, deseq_subset, measure = "log2fold", "diet", timePoint = timePoint, taxa = "Species", threshold = 0.05, LDA = FALSE, FDR = TRUE, customColors = customColors, path = newPath, includeUnknownSpecies = TRUE)
      
      # log2fold change graph based on deseq2 results
      log2foldChangeGraphSingleTimepoint(ps_subset, deseq_subset, timePoint = timePoint, taxa = "Species", threshold = 0.05, customColors = customColors, customPhylaColors = customPhylaColors, path = newPath, dim =c(6,10), includeUnknownSpecies = TRUE)
      
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
        ps_taxa <- tax_glom(ps_subset, taxrank = txnLevel, NArm = FALSE)
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
    pathToSave <- "~/CHUM_git/figures/Thibault_dss/relative_abundance_dss_diet_all_groups/"
    existingDirCheck(pathToSave)
    
    #customColors for graph display
    customColors = c("#95BECF","#F2AA84","#325BAD","#B22222")
    
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
                  dim = c(4,4.5), displayPvalue = FALSE, displaySignificance = TRUE, includeUnknownSpecies = TRUE, additionnalAes =
                    list(scale_x_discrete(labels = c("50 ppm\nCtrl","500 ppm\nCtrl","50 ppm\nDSS","500 ppm\nDSS")),
                         my_theme(),
                         labs(color = "", x=""))) # Include axis lines  # Include axis bar)
    }
    
    #customColors for graph display
    customColors = c("#95BECF","#F2AA84","#325BAD","#B22222")
    
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
        ps_taxa <- tax_glom(ps_subset, taxrank = txnLevel, NArm = FALSE)
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

}

# Manual annotation of ASVs species using BLAST
{
  as.character(refseq(ps_flt_diet)["ASV213"]) # Gordonibacter urolithinfaciens 
  as.character(refseq(ps_flt_diet)["ASV166"]) # Lactobacillus gasseri
  as.character(refseq(ps_flt_diet)["ASV200"]) # Lactobacillus gasseri
  as.character(refseq(ps_flt_diet)["ASV8"])  # Ligilactobacillus murinus
  as.character(refseq(ps_flt_diet)["ASV38"])  # Ligilactobacillus murinus
}

# Correlations
{
  
  # With differentially abundant bugs at end timepoint
  {
  df <- as.data.frame(read_excel("~/CHUM_git/gut-microbiota-iron/experiments/data_both_exp.xlsx"))
  df$id <- substring(df$id, first = 1, last = 5)
  df <- df[df$exp == "dss",]
  df <- df[-46,]
  df <- df[df$treatment == "dss",]
  df$diet <- factor(df$diet, levels = c(50,500))
  rownames(df) <- df$id
  rownames(df) <- paste0(rownames(df), "_Tfinal")
  df <- df[,-c(1:7,9,12:14,16)]
  
  colnames(df) <- c("Stool iron at\nend of exposure","Spleen iron","Liver iron","DAI")
  
  {
    ps_counts = transformCounts(tax_glom(ps_t35_flt, taxrank = "Family", NArm = FALSE), transformation = "rel_ab")
    counts <- t(otu_table(ps_counts))[,c("ASV2","ASV3")]
    row.names(counts) <- str_replace(row.names(counts), pattern = "_T35", replacement = "_Tfinal")
    df <- merge(df, counts, by = "row.names")
    rownames(df) <- df$Row.names
    df <- df[,-1]
  }
  
  ps_sub <- prune_samples(sample_data(ps_tfinal_flt)$treatment %in% c("dss"), ps_tfinal_flt)
  tax_table(ps_sub)["ASV6","Species"] <- "Rodentium"
  tax_table(ps_sub)["ASV8","Species"] <- "Murinus"
  tax_table(ps_sub)["ASV38","Species"] <- "Murinus"
    
  correlation2VarSelectASVs(ps_sub, measure = "log2fold", asvList = c("ASV38","ASV8","ASV6","ASV1"),
                            taxa = "Species", threshold = 0.05, FDR = TRUE, displayPvalue = TRUE,
                            path = "~/CHUM_git/figures/Thibault_dss/cor_hmaps/", df = df, showIndivCor = FALSE, transformation = "rel_ab",
                            displayOnlySig = FALSE, returnMainFig = TRUE, displaySpeciesASVNumber = TRUE, colorsHmap = c("#d1762d","#639381"))+
    my_theme()+
    theme(axis.line = element_blank(),
          legend.background = element_blank())
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/cor_hmap_speciest112.png", width = 7, height = 5, dpi = 500, bg = "white")
  
  df <- as.data.frame(read_excel("~/CHUM_git/gut-microbiota-iron/experiments/data_both_exp.xlsx"))
  df$id <- substring(df$id, first = 1, last = 5)
  df <- df[df$exp == "dss",]
  df <- df[-46,]
  df <- df[df$treatment == "dss",]
  df$diet <- factor(df$diet, levels = c(50,500))
  rownames(df) <- df$id
  rownames(df) <- paste0(rownames(df), "_Tfinal")
  df <- df[!is.na(df$lcn2),]
  df <- df[,-c(1:7,9,12:13,16)]
  colnames(df) <- c("Stool iron at\nend of exposure","Spleen iron","Liver iron","Fecal\nLCN2","DAI")
  
  ps_sub <- prune_samples(sample_data(ps_tfinal_flt)$treatment %in% c("dss"), ps_tfinal_flt)
  tax_table(ps_sub)["ASV6","Species"] <- "Rodentium"
  tax_table(ps_sub)["ASV8","Species"] <- "Murinus"
  tax_table(ps_sub)["ASV38","Species"] <- "Murinus"
  
  correlation2VarSelectASVs(ps_sub, measure = "log2fold", asvList = c("ASV1","ASV6","ASV8","ASV38"),
                            taxa = "Species", threshold = 0.05, FDR = TRUE, displayPvalue = TRUE,
                            path = "~/CHUM_git/figures/Thibault_dss/cor_hmaps/", df = df, showIndivCor = FALSE, transformation = "rel_ab",
                            displayOnlySig = FALSE, returnMainFig = TRUE, displaySpeciesASVNumber = TRUE, colorsHmap = c("#d1762d","#639381"))+
    my_theme()+
    theme(axis.line = element_blank(),
          legend.background = element_blank())
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/cor_hmap_speciest112_LCN2.png", width = 8.75, height = 5, dpi = 500, bg = "white")
  }
  
  # With some differentially abundant ASVs at t49 - B. acidifaciens asv70
  {
    df <- as.data.frame(read_excel("~/CHUM_git/gut-microbiota-iron/experiments/data_both_exp.xlsx"))
    df$id <- substring(df$id, first = 1, last = 5)
    df <- df[df$exp == "dss",]
    df <- df[-46,]
    df <- df[df$treatment == "dss",]
    df$diet <- factor(df$diet, levels = c(50,500))
    rownames(df) <- df$id
    rownames(df) <- paste0(rownames(df), "_T49")
    df <- df[,-c(1:7,9,12:14,16)]
    colnames(df) <- c("Stool iron at\nend of exposure","Spleen iron","Liver iron","DAI")
    
    ps_sub <- prune_samples(sample_data(ps_t49_flt)$treatment %in% c("dss"), ps_t49_flt)
    
    correlation2VarSelectASVs(ps_sub, measure = "log2fold", asvList = c("ASV70","ASV64"),
                              taxa = "Species", threshold = 0.05, FDR = TRUE, displayPvalue = TRUE,
                              path = "~/CHUM_git/figures/Thibault_dss/cor_hmaps/", df = df, showIndivCor = FALSE, transformation = "rel_ab",
                              displayOnlySig = FALSE, returnMainFig = TRUE, displaySpeciesASVNumber = TRUE, colorsHmap = c("#d1762d","#639381"))+
      my_theme()+
      theme(axis.line = element_blank(),
            legend.background = element_blank())
  }
  
  # With Lacto family and other families that were differentially abundant at t35
  {
    df <- as.data.frame(read_excel("~/CHUM_git/gut-microbiota-iron/experiments/data_both_exp.xlsx"))
    df$id <- substring(df$id, first = 1, last = 5)
    df <- df[df$exp == "dss",]
    df <- df[-46,]
    df <- df[df$treatment == "dss",]
    df$diet <- factor(df$diet, levels = c(50,500))
    rownames(df) <- df$id
    rownames(df) <- paste0(rownames(df), "_T35")
    # df <- df[!is.na(df$lcn2),]
    df <- df[,-c(1:7,9,12:14,16)]
    colnames(df) <- c("Stool iron at\nend of exposure","Spleen iron","Liver iron","DAI")
    
    ps_sub <- prune_samples(sample_data(ps_t35_flt)$treatment %in% c("dss"), ps_t35_flt)
    ps_sub <- tax_glom(ps_sub, taxrank = "Family", NArm = FALSE)
    
    correlation2VarSelectASVs(ps_sub, measure = "log2fold", asvList = c("ASV2","ASV3","ASV9","ASV80","ASV7","ASV194"),
                              taxa = "Family", threshold = 0.05, FDR = TRUE, displayPvalue = TRUE,
                              path = "~/CHUM_git/figures/Thibault_dss/cor_hmaps/", df = df, showIndivCor = FALSE, transformation = "rel_ab",
                              displayOnlySig = FALSE, returnMainFig = TRUE, displaySpeciesASVNumber = TRUE, colorsHmap = c("#d1762d","#639381"))+
      my_theme()+
      theme(axis.line = element_blank(),
            legend.background = element_blank())
    
    
    ps_sub <- prune_samples(sample_data(ps_t35_flt)$treatment %in% c("dss"), ps_t35_flt)
    ps_sub <- tax_glom(ps_sub, taxrank = "Genus", NArm = FALSE)
    
    correlation2VarSelectASVs(ps_sub, measure = "log2fold", asvList = c("ASV8","ASV9","ASV252"),
                              taxa = "Genus", threshold = 0.05, FDR = TRUE, displayPvalue = TRUE,
                              path = "~/CHUM_git/figures/Thibault_dss/cor_hmaps/", df = df, showIndivCor = FALSE, transformation = "rel_ab",
                              displayOnlySig = FALSE, returnMainFig = TRUE, displaySpeciesASVNumber = TRUE, colorsHmap = c("#d1762d","#639381"))+
      my_theme()+
      theme(axis.line = element_blank(),
            legend.background = element_blank())
    
    
    
    
    
  }
  
}

# Correlation but with qPCR results
{
  # All variables
  df_pcr <- as.data.frame(read_excel("~/CHUM_git/gut-microbiota-iron/experiments/finished exp/young-DSS-exp3/pcr_validation/F rodentium RT-PCR YoungDSS Exp.xlsx"))
  colnames(df_pcr) <- df_pcr[3,]
  df_pcr <- df_pcr[!is.na(df_pcr$`Muribaculum intestinale`),]
  df_pcr <- df_pcr[-1,]
  df_pcr$lm <- as.numeric(df_pcr$`Lm/16S`)
  df_pcr$mi <- as.numeric(df_pcr$`Mi/16S`)
  df_pcr$fr <- as.numeric(df_pcr$`Fr/16S`)
  df_pcr <- df_pcr[,-c(5,10,15)]
  row.names(df_pcr) <- df_pcr$`Mouse ID`
  df_pcr <- df_pcr[,c(17:19)]
  
  df_var <- as.data.frame(read_excel("~/CHUM_git/gut-microbiota-iron/experiments/data_both_exp.xlsx"))
  df_var$id <- substring(df_var$id, first = 1, last = 5)
  df_var <- df_var[df_var$exp == "dss",]
  df_var <- df_var[-46,]
  df_var <- df_var[df_var$treatment == "dss",]
  df_var$diet <- factor(df_var$diet, levels = c(50,500))
  rownames(df_var) <- df_var$id
  df_var <- df_var[,-c(1:7,9,12:14,16)]
  colnames(df_var) <- c("Stool iron at\nend of exposure","Spleen iron","Liver iron","DAI")
  
  {
    ps_counts = transformCounts(tax_glom(ps_t35_flt, taxrank = "Family", NArm = FALSE), transformation = "rel_ab")
    counts <- t(otu_table(ps_counts))[,c("ASV2","ASV3","ASV9")]
    row.names(counts) <- str_replace(row.names(counts), pattern = "_T35", replacement = "")
    df_var<- merge(df_var, counts, by = "row.names")
    rownames(df_var) <- df_var$Row.names
    df_var <- df_var[,-1]
  }
  
  df_cor <- merge(df_pcr, df_var, by = "row.names")
  row.names(df_cor) <- df_cor[,1] 
  df_cor <- df_cor[,-1]
  
  cor_results <- rcorr(as.matrix(df_cor), type = "spearman")
  cor_matrix <- cor_results$r  # Extract correlation coefficients
  cor_matrix <- cor_matrix[4:10,1:3]
  p_values <- cor_results$P    # Extract p-values
  p_values <- p_values[4:10,1:3]
  
  # Specific to LCN2
  df_var <- as.data.frame(read_excel("~/CHUM_git/gut-microbiota-iron/experiments/data_both_exp.xlsx"))
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
  
  ggplot(cor_melt, aes(x = Variables, y = Species, fill = value)) + 
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
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/cor_hmap.png", width = 9.5, height = 3, dpi = 500, bg = "white")
}

# Additionnal figures
{
  # B. acidifaciens abundance at t49 and t54 (ASV70)
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
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/baT54.png", plot = baT54, width = 5, height = 3, dpi = 500)
}

# Supplementary chronobiome figures
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
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/chronobiome_suppl/b.png", plot = b, width = 9, height = 5, dpi = 500)
  
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
  
  ggsave(filename = "~/CHUM_git/figures/memoire/dss/chronobiome_suppl/e.png", plot = e, width = 14, height = 6, dpi = 500)
  
  
}

# Additionnal correlations
{
  # Load dataframes of variable
  df_pcr <- as.data.frame(read_excel("~/CHUM_git/gut-microbiota-iron/experiments/finished exp/young-DSS-exp3/pcr_validation/F rodentium RT-PCR YoungDSS Exp.xlsx"))
  colnames(df_pcr) <- df_pcr[3,]
  df_pcr <- df_pcr[!is.na(df_pcr$`Muribaculum intestinale`),]
  df_pcr <- df_pcr[-1,]
  df_pcr$lm <- as.numeric(df_pcr$`Lm/16S`)
  df_pcr$mi <- as.numeric(df_pcr$`Mi/16S`)
  df_pcr$fr <- as.numeric(df_pcr$`Fr/16S`)
  df_pcr <- df_pcr[,-c(5,10,15)]
  row.names(df_pcr) <- df_pcr$`Mouse ID`
  df_pcr <- df_pcr[,c(17:19)]
  
  df_var <- as.data.frame(read_excel("~/CHUM_git/gut-microbiota-iron/experiments/data_both_exp.xlsx"))
  df_var$id <- substring(df_var$id, first = 1, last = 5)
  df_var <- df_var[df_var$exp == "dss",]
  df_var <- df_var[-46,]
  rownames(df_var) <- df_var$id
  df_var <- df_var[,-c(1:3,16)]
  colnames(df_var) <- c("Body weight\n(T112)","Colon length\n(T112)","Spleen weight\n(T112)",
                        "Liver weight\n(T112)","Stool iron\n(T35)","Stool iron\n(T112)",
                        "Spleen iron\n(T112)","Liver iron\n(T112)","Body weight\n(T35)",
                        "Body weight\n(T49)","Fecal LCN-2\n(T54)","DAI\n(T54)")
  
  # T35 correlations
  {
    
    # Family level
    df <- df_var[,c(5,9)]
    ps_counts = transformCounts(tax_glom(ps_t35_flt, taxrank = "Family", NArm = FALSE), transformation = "rel_ab")
    
    p <- correlationSelectedASVs(df_var = df, ps = ps_counts, selected_ASVs = c("ASV2","ASV3","ASV9","ASV80","ASV194"),
                            taxrank = "Family", timepoint = "T35")
    # Species level
    tax_table(ps_t35_flt)["ASV213","Species"] <- "urolithinfaciens"
    
    tax_table(ps_t35_flt)["ASV166","Species"] <- "gasseri"
    tax_table(ps_t35_flt)["ASV200","Species"] <- "gasseri"
    tax_table(ps_t35_flt)["ASV166","Genus"] <- "Lactobacillus"
    tax_table(ps_t35_flt)["ASV200","Genus"] <- "Lactobacillus"
    
    tax_table(ps_t35_flt)["ASV8","Species"]  <- "murinus"
    tax_table(ps_t35_flt)["ASV38","Species"]  <- "murinus"
    tax_table(ps_t35_flt)["ASV49","Species"]  <- "acidifaciens"
    
    tax_table(ps_t35_flt)[,"Species"] <- ifelse(is.na(tax_table(ps_t35_flt)[,"Species"]),
                                                paste0(tax_table(ps_t35_flt)[,"Genus"]," unkwown", " (",row.names(tax_table(ps_t35_flt)),")"),
                                                paste0(tax_table(ps_t35_flt)[,"Genus"], " ",tax_table(ps_t35_flt)[,"Species"], " (", row.names(tax_table(ps_t35_flt)),")"))
    
    p <- correlationSelectedASVs(df_var = df, ps = ps_t35_flt, selected_ASVs = c("ASV213","ASV53","ASV70","ASV49","ASV166","ASV200","ASV8","ASV38"),
                                 taxrank = "Species", timepoint = "T35")
    
  }
  
  # T49 correlations
  {
    
    # Species level
    df <- df_var[,c(10:12)]
    
    tax_table(ps_t49_flt)["ASV213","Species"] <- "urolithinfaciens"
    tax_table(ps_t49_flt)[,"Species"] <- ifelse(is.na(tax_table(ps_t49_flt)[,"Species"]),
                                                paste0(tax_table(ps_t49_flt)[,"Genus"]," unkwown", " (",row.names(tax_table(ps_t49_flt)),")"),
                                                paste0(tax_table(ps_t49_flt)[,"Genus"], " ",tax_table(ps_t49_flt)[,"Species"], " (", row.names(tax_table(ps_t49_flt)),")"))
    
    p <- correlationSelectedASVs(df_var = df, ps = ps_t49_flt, selected_ASVs = c("ASV213","ASV70","ASV9"),
                                 taxrank = "Species", timepoint = "T49")
    # Species level - DSS only
    df <- df_var[!is.na(df_var[,11]),]
    df <- df[,c(10:12)]
    
    df <- df_var[!is.na(df_var[,12]),]
    df <- df[,c(10:12)]
    
    p <- correlationSelectedASVs(df_var = df, ps = ps_t49_flt, selected_ASVs = c("ASV213","ASV70","ASV9"),
                                 taxrank = "Species", timepoint = "T49")
    p

    
  }
  
  # T54 correlations
  {
    
    # Species level - DSS only
    df <- df_var[!is.na(df_var[,12]),]
    df <- df[,c(11:12)]
    
    tax_table(ps_t54_flt)["ASV8","Species"] <- "murinus"
    tax_table(ps_t54_flt)["ASV38","Species"] <- "murinus"
    tax_table(ps_t54_flt)["ASV19","Species"] <- "intestinale"
    tax_table(ps_t54_flt)["ASV6","Species"] <- "rodentium"
    tax_table(ps_t54_flt)["ASV1","Species"] <- "rodentium"
    
    tax_table(ps_t54_flt)[,"Species"] <- ifelse(is.na(tax_table(ps_t54_flt)[,"Species"]),
                                                paste0(tax_table(ps_t54_flt)[,"Genus"]," unkwown", " (",row.names(tax_table(ps_t54_flt)),")"),
                                                paste0(tax_table(ps_t54_flt)[,"Genus"], " ",tax_table(ps_t54_flt)[,"Species"], " (", row.names(tax_table(ps_t54_flt)),")"))
    
    p <- correlationSelectedASVs(df_var = df, ps = ps_t54_flt, selected_ASVs = c("ASV1","ASV6","ASV8","ASV38","ASV19"),
                                 taxrank = "Species", timepoint = "T54")
    p

    # Family level
    ps_counts = transformCounts(tax_glom(ps_t54_flt, taxrank = "Family", NArm = FALSE), transformation = "rel_ab")
    View(tax_table(ps_counts))
    
    p <- correlationSelectedASVs(df_var = df, ps = ps_counts, selected_ASVs = c("ASV1","ASV2","ASV4","ASV5",
                                                                                "ASV7","ASV22","ASV57","ASV80",
                                                                                "ASV87","ASV124","ASV164","ASV241",
                                                                                "ASV242","ASV246","ASV390"),
                                 taxrank = "Family", timepoint = "T54")
    p
    
    
  }
  
  # Tfinal correlations
  {
    
    # Species level - DSS only
    df <- df_var[is.na(df_var[,12]),]
    df <- df[,c(1:4,6:8)]
    
    tax_table(ps_tfinal_flt)["ASV8","Species"] <- "murinus"
    tax_table(ps_tfinal_flt)["ASV38","Species"] <- "murinus"
    tax_table(ps_tfinal_flt)["ASV19","Species"] <- "intestinale"
    tax_table(ps_tfinal_flt)["ASV6","Species"] <- "rodentium"
    tax_table(ps_tfinal_flt)["ASV1","Species"] <- "rodentium"
    
    tax_table(ps_tfinal_flt)[,"Species"] <- ifelse(is.na(tax_table(ps_tfinal_flt)[,"Species"]),
                                                paste0(tax_table(ps_tfinal_flt)[,"Genus"]," unkwown", " (",row.names(tax_table(ps_tfinal_flt)),")"),
                                                paste0(tax_table(ps_tfinal_flt)[,"Genus"], " ",tax_table(ps_tfinal_flt)[,"Species"], " (", row.names(tax_table(ps_tfinal_flt)),")"))
    
    p <- correlationSelectedASVs(df_var = df, ps = ps_tfinal_flt, selected_ASVs = c("ASV1","ASV6","ASV8","ASV38","ASV19"),
                                 taxrank = "Species", timepoint = "Tfinal")
    p
    
    # Family level
    ps_counts = transformCounts(tax_glom(ps_t54_flt, taxrank = "Family", NArm = FALSE), transformation = "rel_ab")
    View(tax_table(ps_counts))
    
    p <- correlationSelectedASVs(df_var = df, ps = ps_counts, selected_ASVs = c("ASV1","ASV2","ASV4","ASV5",
                                                                                "ASV7","ASV22","ASV57","ASV80",
                                                                                "ASV87","ASV124","ASV164","ASV241",
                                                                                "ASV242","ASV246","ASV390"),
                                 taxrank = "Family", timepoint = "T54")
    p
    
    
  }
  


}
}





