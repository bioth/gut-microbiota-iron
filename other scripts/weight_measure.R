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
  library("readxl")
  library("esquisse")
}



#If working from huawei pc
setwd("D:/CHUM_git/gut-microbiota-iron/")
#If working from CHUM pc
setwd("I:/Chercheurs/Santos_Manuela/Thibault M/gut-microbiota-iron/")
#If working from la bête
setwd("C:/Users/Thibault/Documents/CHUM_git/gut-microbiota-iron/")

#loading functions for data manipulation
source("other scripts/dataManipFunctions.R")
  
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
setwd("../gut-microbiota-iron/young-DSS-exp3")
young_weight <- read.csv("young48_weight.csv", header = TRUE, sep = ";")

#data manipulation
young_weight <- weightDataManipulation(young_weight,4, fromDay0 = TRUE)

#creating scatter plot with the four different treatments (diet combined with dss or control)
young_weight_plot <- weightPlot(young_weight, percentage = TRUE, diet_only = FALSE)
young_weight_plot


#Loading weight measure file
setwd("../young-DSS-exp3/")
young_weight <- read.csv("young48_weight_cageChanges.csv", header = TRUE, sep = ";")

#removing first week body weight measurements
young_weight <- young_weight[,-c(6:10)]

#Manipulating weight measures data
young_weight <- weightDataManipulation(young_weight,4, fromDay0 = TRUE)

#creating scatter plot with the four different treatments (diet combined with dss or control)
young_weight_plot <- weightPlot(young_weight, percentage = TRUE, diet_only = FALSE)
young_weight_plot

#saving scatter plot
setwd("../figures/young48")
ggsave("young_weight_percentage.png", width = 11, height = 6, dpi = 300, bg = "white")








#Loading weight measure file
setwd("../../adult-abx-exp5/")
young_weight <- read.csv("weighabmeeting.csv", header = TRUE, sep = ";")

#Manipulating weight measures data
young_weight <- weightDataManipulation(young_weight,4, fromDay0 = TRUE)

#creating scatter plot with the four different treatments (diet combined with dss or control)
young_weight_plot <- weightPlot(young_weight, percentage = TRUE, diet_only = TRUE)
young_weight_plot

#saving scatter plot
setwd("../figures/young48")
ggsave("yapaletemps.png", width = 11, height = 6, dpi = 300, bg = "white")


















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

#saving figure
setwd("../figures/young48/")
ggsave("young48_dIndex_evolution.png", width = 5, height = 5, dpi = 300, bg = "white")

#DSI measures at final day of DSS
young_final_DSI_plot <- dssDsiFinalDay(young_dss_followup)
young_final_DSI_plot
ggsave("young48_final_dIndex.png", width = 5, height = 5, dpi = 300, bg = "white")





#adult mice dss experiments
###LOADING DATA
setwd("../experiments/ongoing exp/combined exp/")
young_dss_followup <- read.csv("combined_adult_dss_followup.csv", header = TRUE, sep = ";")

young_dss_followup <- dssFollowupManipulation(df = young_dss_followup,groupInfoCols = 4,dateStart = "2024-05-29",nbrDays = 5, negativeOnly = TRUE)


#creating scatter plot with the four different treatments (diet combined with dss or control)
#this graph has a disease index score in the y column
young_dssflwup_plot <- dssDiseaseIndexPlot(young_dss_followup)
young_dssflwup_plot

#saving figure
setwd("../../figures/adult_batch/")
ggsave("adult_dIndex_evolution.png", width = 5, height = 5, dpi = 300, bg = "white")

#DSI measures at final day of DSS
young_final_DSI_plot <- dssDsiFinalDay(young_dss_followup)
young_final_DSI_plot
ggsave("adult_final_dIndex.png", width = 5, height = 5, dpi = 300, bg = "white")









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











###Final dissection data young48###
#loading the data
setwd("../young-DSS-exp3/")
dissec_adult <- read.csv("young48_dss_dissection.csv", sep = ";", header = TRUE)
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
setwd("../figures/meetingurgent")
ggsave(plot = adult_dissec_spln,"spleen_weight.png", width = 6, height = 5, dpi = 300, bg = "white")
ggsave(plot = adult_dissec_lvr,"liver_weight.png", width = 6, height = 5, dpi = 300, bg = "white")
ggsave(plot = adult_dissec_cln,"colon_length.png", width = 6, height = 5, dpi = 300, bg = "white")
ggsave(plot = adult_dissec_bw,"body_wweight.png", width = 6, height = 5, dpi = 300, bg = "white")


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






































#Iron content measurements
#For young48 dss experiment
setwd("experiments/finished exp/young-DSS-exp3/")


#T35 ferrozine assay for stools
df <- read_xlsx("young48_dss_ferrozine_t35.xlsx")
df <- df[,c(2,14:16)]
colnames(df) <- df[2,]
colnames(df)[1:2] <- c("id","iron_concentration")
df <- df[-c(1,2,27),]
df$gg_group <- paste(df$treatment, "+", df$diet, sep = "")
df$gg_group <- factor(df$gg_group, levels = c("water+50","dss+50","water+500","dss+500"))
df$iron_concentration <- as.numeric(df$iron_concentration)

ironBoxplot(df, "iron_concentration", display_significance_bars = F, title = "Iron concentration in stools at day 35", y_axis_title = "yg of iron per g of stools", custom_colors = c("blue","red"), path = "D:/CHUM_git/figures/iron_measures")

#Fit a ANOVA model for the current time point
anova <- aov(iron_concentration ~ diet * treatment, data = df)
summary(anova)
# Perform Tukey's HSD test and store the results in the list
results <- TukeyHSD(anova)
print(results)

#Tf ferrozine assay for liver
#They are multiple sheets so we use a custom function
sheets <- read_excel_allsheets("young48_dss_ferrozine_tF.xlsx")

df <- as.data.frame(sheets["Ferrozine Liver"])
df <- df[1:51,c(3,6,8,14:18)]
colnames(df) <- df[2,]
colnames(df)[1:4] <- c("id","wet_weight","dry_weight","iron_concentration")
df <- df[-c(1,2,27),]
df$gg_group <- paste(df$treatment, "+", df$diet, sep = "")
df$gg_group <- factor(df$gg_group, levels = c("water+50","dss+50","water+500","dss+500"))
# df$gg_group <- factor(df$gg_group, levels = c("water+50","water+500","dss+50","dss+500")) #For Claire comparaison with BalbC mice
df$iron_concentration <- as.numeric(df$iron_concentration)
df$liver_weight <- as.numeric(df$liver_weight)
df$wet_weight <- as.numeric(df$wet_weight)
df$dry_weight <- as.numeric(df$dry_weight)
df$dry_to_wet_ratio <- df$dry_weight/df$wet_weight


#To calculate total iron in organ = iron concentration per g of dry weight*total organ weight
#need to take into account the wet to dry ratio!
df$total_iron <- df$iron_concentration*df$liver_weight*df$dry_to_wet_ratio


ironBoxplot(df, "total_iron", display_significance_bars = F, title = "Total liver iron at final day", y_axis_title = "yg of iron", custom_colors = c("blue","red"), path = "D:/CHUM_git/figures/iron_measures")

#Fit a ANOVA model for the current time point
anova <- aov(total_iron ~ diet * treatment, data = df)
summary(anova)
# Perform Tukey's HSD test and store the results in the list
results <- TukeyHSD(anova)
print(results)



#Tf ferrozine assay for spleen
df <- as.data.frame(sheets["Ferrozine Spleen"])
df <- df[1:51,c(3,6,8,14:18)]
colnames(df) <- df[2,]
colnames(df)[1:4] <- c("id","wet_weight","dry_weight","iron_concentration")
df <- df[-c(1,2,27),]
df$gg_group <- paste(df$treatment, "+", df$diet, sep = "")
df$gg_group <- factor(df$gg_group, levels = c("water+50","dss+50","water+500","dss+500"))
df$iron_concentration <- as.numeric(df$iron_concentration)
df$spleen_weight <- as.numeric(df$spleen_weight)
df$wet_weight <- as.numeric(df$wet_weight)
df$dry_weight <- as.numeric(df$dry_weight)
df$dry_to_wet_ratio <- df$wet_weight/df$dry_weight


#To calculate total iron in organ = iron concentration per g of dry weight*total organ weight
#need to take into account the wet to dry ratio!
df$total_iron <- df$iron_concentration*df$spleen_weight*df$dry_to_wet_ratio

ironBoxplot(df, "total_iron", display_significance_bars = F, title = "Total spleen iron at final day", y_axis_title = "yg of iron", custom_colors = c("blue","red"), path = "D:/CHUM_git/figures/iron_measures")

#Fit a ANOVA model for the current time point
anova <- aov(iron_concentration ~ diet * treatment, data = df)
summary(anova)
# Perform Tukey's HSD test and store the results in the list
results <- TukeyHSD(anova)
print(results)


#Tf ferrozine assay for stools
df <- as.data.frame(sheets["Ferrozine Stool Tfinal"])
df <- df[4:53,c(3,14:16)]
colnames(df) <- df[1,]
colnames(df)[1:2] <- c("id","iron_concentration")
df <- na.omit(df)
df <- df[-1,]
df$gg_group <- paste(df$treatment, "+", df$diet, sep = "")
df$gg_group <- factor(df$gg_group, levels = c("water+50","dss+50","water+500","dss+500"))
df$iron_concentration <- as.numeric(df$iron_concentration)

ironBoxplot(df, "iron_concentration", display_significance_bars = F, title = "Iron concentration in stools at final day", y_axis_title = "yg of iron per g of stools", custom_colors = c("blue","red"), path = "D:/CHUM_git/figures/iron_measures")

#Fit a ANOVA model for the current time point
anova <- aov(iron_concentration ~ diet * treatment, data = df)
summary(anova)
# Perform Tukey's HSD test and store the results in the list
results <- TukeyHSD(anova)
print(results)

