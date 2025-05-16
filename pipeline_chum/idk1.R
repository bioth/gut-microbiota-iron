#loading libraries
{
  library("tidyverse") #loading bunch of packages
  library("ggplot2") #come on, everyone knows what it is used for
  library("dplyr") #arranging and manipulating data easily
  library("geepack") #library for loading GEE tests
  library("lme4") #library for loading ANOVA
  library("car") #for anova too
  library("ggsignif") #adding significance bars to ggplots
  library("Hmisc") #idk lol
  library("esquisse")
}

#If working from CHUM pc
setwd("I:/Chercheurs/Santos_Manuela/Thibault M/gut-microbiota-iron/")
source("other scripts/dataManipFunctions.R")


#repeated young mice experiment
###LOADING YOUNG MICE DATA
young_dss_followup <- read.csv("experiments/finished exp/young-DSS-exp3/young48_dss_followup.csv", header = TRUE, sep = ";")
young_dss_followup <- dssFollowupManipulation(df = young_dss_followup,groupInfoCols = 4,dateStart = "2024-05-29",nbrDays = 5, negativeOnly = TRUE)


#creating scatter plot with the four different treatments (diet combined with dss or control)
#this graph has a disease index score in the y column
young_dssflwup_plot <- dssDiseaseIndexPlot(young_dss_followup)
young_dssflwup_plot

#saving figure
ggsave("figures/young48/young48_dIndex_evolution.png", width = 5, height = 5, dpi = 300, bg = "white")

#DSI measures at final day of DSS
young_final_DSI_plot <- dssDsiFinalDay(young_dss_followup)
young_final_DSI_plot
ggsave("young48_final_dIndex.png", width = 5, height = 5, dpi = 300, bg = "white")


#adult mice dss experiments
###LOADING DATA
adult_dss_followup <- read.csv("experiments/ongoing exp/combined exp/combined_adult_dss_followup.csv", header = TRUE, sep = ";")
adult_dss_followup <- dssFollowupManipulation(df = adult_dss_followup,groupInfoCols = 4,dateStart = "2024-05-29",nbrDays = 5, negativeOnly = TRUE)


#creating scatter plot with the four different treatments (diet combined with dss or control)
#this graph has a disease index score in the y column
adult_dssflwup_plot <- dssDiseaseIndexPlot(adult_dss_followup)
adult_dssflwup_plot

#saving figure
ggsave("figures/adult_batch/adult_dIndex_evolution.png", width = 5, height = 5, dpi = 300, bg = "white")

#DSI measures at final day of DSS
adult_final_DSI_plot <- dssDsiFinalDay(adult_dss_followup)
adult_final_DSI_plot
ggsave("adult_final_dIndex.png", width = 5, height = 5, dpi = 300, bg = "white")