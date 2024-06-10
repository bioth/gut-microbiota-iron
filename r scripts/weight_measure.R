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



#If working from huawei pc
setwd("D:/CHUM_git/gut-microbiota-iron/")
#If working from CHUM pc
setwd("I:/Chercheurs/Santos_Manuela/Thibault M/gut-microbiota-iron/")
#If working from la bête
setwd("C:/Users/Thibault/Documents/CHUM_git/gut-microbiota-iron/")

{
  #loading functions for data manipulation
  setwd("r scripts/")
  source("dataManipFunctions.R")
}
  
#Loading weight measure file
setwd("../adult-DSS-exp/")
adult_weight <- read.csv("adult36dss_weight_measurement.csv", header = TRUE, sep = ";")

#Manipulating weight measures data
adult_weight <- weightDataManipulation(adult_weight,4)

#replacing abnormal values
adult_weight$weight[4] <- 24.2
adult_weight$weight[28] <- 23.1

#using this enables to verify if date variable is of type "date" and weight of type numeric
str(adult_weight)

#creating scatter plot with the four different treatments (diet combined with dss or control)
adult_weight_plot <- weightPlot(adult_weight, percentage = TRUE)
adult_weight_plot

#saving scatter plot
setwd("../figures/adult36")
ggsave("adults_weight.png", width = 6, height = 6, dpi = 300, bg = "white")
  
#YOUNG MICE
#Loading weight measure file for young mice
setwd("../../young-DSS-exp2/")
young_weight <- read.csv("young32dss_weight.csv", header = TRUE, sep = ";")

#data manipulation
young_weight <- weightDataManipulation(young_weight,4)

#creating scatter plot with the four different treatments (diet combined with dss or control)
young_weight_plot <- weightPlot(young_weight, percentage = TRUE)
young_weight_plot

#saving scatter plot
setwd("../figures/young32/")
ggsave("young32_weight.png", width = 10, height = 6, dpi = 300, bg = "white")


#Young 48 mice
setwd("../young-DSS-exp3")
young_weight <- read.csv("young48_weight.csv", header = TRUE, sep = ";")

#data manipulation
young_weight <- weightDataManipulation(young_weight,4)

#creating scatter plot with the four different treatments (diet combined with dss or control)
young_weight_plot <- weightPlot(young_weight, percentage = FALSE, diet_only = FALSE)
young_weight_plot


#Loading weight measure file
setwd("../young-DSS-exp3/")
young_weight <- read.csv("young48_weight_cageChanges.csv", header = TRUE, sep = ";")

#Manipulating weight measures data
young_weight <- weightDataManipulation(young_weight,4)

#creating scatter plot with the four different treatments (diet combined with dss or control)
young_weight_plot <- weightPlot(young_weight, percentage = FALSE, diet_only = FALSE)
young_weight_plot

#saving scatter plot
setwd("../figures/young48")
ggsave("young_weight_percent.png", width = 11, height = 6, dpi = 300, bg = "white")








###DSS FOLLOW UP SHEET DATA
#loading the dss followup sheet data
setwd("../adult-DSS-exp/")
adult_dss_followup <- read.csv("adult36dss_followup.csv", header = TRUE, sep = ";")
setwd("../young-DSS-exp2/")
young_dss_followup <- read.csv("young32dss_followup.csv", header = TRUE, sep = ";")


#loading the functions
setwd("../r scripts/")
source("dataManipFunctions.R")

###LOADING ADULT MICE DATA
adult_dss_followup[,21][adult_dss_followup[,21] == "Yes"] <- "G" #Checking if changing this changes anything
adult_dss_followup <- dssFollowupManipulation(df = adult_dss_followup,groupInfoCols = 4,dateStart = "2023-12-04",nbrDays = 5, negativeOnly = FALSE) #negative only FALSE if absolute differences in weight are taken into account 


###LOADING YOUNG MICE DATA
young_dss_followup <- dssFollowupManipulation(df = young_dss_followup,groupInfoCols = 4,dateStart = "2024-02-14",nbrDays = 5, negativeOnly = FALSE)

#creating scatter plot with the four different treatments (diet combined with dss or control)
#this graph has a disease index score in the y column
adult_dssflwup_plot <- dssDiseaseIndexPlot(adult_dss_followup)
adult_dssflwup_plot

#DSI measures at final day of DSS
adult_final_DSI_plot <- dssDsiFinalDay(adult_dss_followup)
adult_final_DSI_plot

#saving figure
setwd("../figures/adult36/")
ggsave("adult36_dIndex.png", width = 9, height = 5, dpi = 300, bg = "white")

young_dssflwup_plot <- dssDiseaseIndexPlot(young_dss_followup)
young_dssflwup_plot

young_final_DSI_plot <- dssDsiFinalDay(young_dss_followup)
young_final_DSI_plot

#saving figure
setwd("../young32/")
ggsave("young32_dIndex.png", width = 9, height = 5, dpi = 300, bg = "white")






#repeated young mice experiment
###LOADING YOUNG MICE DATA
setwd("../young-DSS-exp3/")
young_dss_followup <- read.csv("young48_dss_followup.csv", header = TRUE, sep = ";")

young_dss_followup <- dssFollowupManipulation(df = young_dss_followup,groupInfoCols = 4,dateStart = "2024-05-29",nbrDays = 5, negativeOnly = TRUE)


#creating scatter plot with the four different treatments (diet combined with dss or control)
#this graph has a disease index score in the y column
young_dssflwup_plot <- dssDiseaseIndexPlot(young_dss_followup)
young_dssflwup_plot

#DSI measures at final day of DSS
young_final_DSI_plot <- dssDsiFinalDay(young_dss_followup)
young_final_DSI_plot




###Final dissection data###
#loading the data
setwd("../adult-DSS-exp/")
dissec_adult <- read.csv("adult36dss_dissection.csv", sep = ";", header = TRUE)
setwd("../r scripts/")
source("dataManipFunctions.R")

dissec_adult <- dissectionDataManipulation(dissec_adult,4)

#boxplot for body weight
adult_dissec_bw <- dissecBoxplot(dissec_adult,"body",display_significance_bars = FALSE) 
adult_dissec_bw

#boxplot for std liver weight
adult_dissec_lvr <- dissecBoxplot(dissec_adult,"liver") 
adult_dissec_lvr

#boxplot for std spleen weight
adult_dissec_spln <- dissecBoxplot(dissec_adult,"spleen") 
adult_dissec_spln

#boxplot for colon length (non std)
adult_dissec_cln <- dissecBoxplot(dissec_adult,"colon") 
adult_dissec_cln


#saving figures
setwd("../figures/adult36")
ggsave(plot = adult_dissec_spln,"spleen_weight.png", width = 8, height = 5, dpi = 300, bg = "white")
ggsave(plot = adult_dissec_lvr,"liver_weight.png", width = 9, height = 5, dpi = 300, bg = "white")
ggsave(plot = adult_dissec_cln,"colon_length.png", width = 7, height = 5, dpi = 300, bg = "white")
ggsave(plot = adult_dissec_bw,"body_wweight.png", width = 9, height = 5, dpi = 300, bg = "white")


#loading the data
setwd("../../young-DSS-exp2/")
dissec_young <- read.csv("young32dss_dissection.csv", sep = ";", header = TRUE)
setwd("../r scripts/")
source("dataManipFunctions.R")

dissec_young <- dissectionDataManipulation(dissec_young,4)

#boxplot for body weight
young_dissec_bw <- dissecBoxplot(dissec_young,"body",display_significance_bars = FALSE) 
young_dissec_bw

#boxplot for std liver weight
young_dissec_lvr <- dissecBoxplot(dissec_young,"liver") 
young_dissec_lvr

#boxplot for std spleen weight
young_dissec_spln <- dissecBoxplot(dissec_young,"spleen") 
young_dissec_spln

#boxplot for colon length (non std)
young_dissec_cln <- dissecBoxplot(dissec_young,"colon") 
young_dissec_cln


#saving figures
setwd("../figures/young32")
ggsave(plot = young_dissec_spln,"spleen_weight.png", width = 8, height = 5, dpi = 300, bg = "white")
ggsave(plot = young_dissec_lvr,"liver_weight.png", width = 9, height = 5, dpi = 300, bg = "white")
ggsave(plot = young_dissec_cln,"colon_length.png", width = 7, height = 5, dpi = 300, bg = "white")
ggsave(plot = young_dissec_bw,"body_weight.png", width = 9, height = 5, dpi = 300, bg = "white")
