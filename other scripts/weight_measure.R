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
setwd("D:/wow these are nice graphs/gut-microbiota-iron/")
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
setwd("adult-DSS-exp/")
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
setwd("experiments/finished exp/young-DSS-exp3/")
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




# All graphs regarding young mice colitis experiment 
setwd("experiments/finished exp/young-DSS-exp3/")

# Mice weight evolution
young_weight <- read.csv("young48_weight_cageChanges.csv", header = TRUE, sep = ";")
young_weight <- young_weight[,-c(6:10)] #removing first week body weight measurements
young_weight <- weightDataManipulation(young_weight,4, fromDay0 = TRUE) #Manipulating weight measures data

# creating scatter plot with the four different treatments (diet combined with dss or control)
young_weight_plot <- weightPlot(young_weight, percentage = FALSE, diet_only = FALSE)
young_weight_plot

# Mice DAI evolution
young_dss_followup <- read.csv("young48_dss_followup.csv", header = TRUE, sep = ";")
young_dss_followup <- dssFollowupManipulation(df = young_dss_followup,groupInfoCols = 4,dateStart = "2024-05-29",nbrDays = 5, negativeOnly = TRUE)

# creating scatter plot with the four different treatments (diet combined with dss or control)
young_dssflwup_plot <- dssDiseaseIndexPlot(young_dss_followup)
ggsave(plot = young_dssflwup_plot, filename = "~/Documents/CHUM_git/figures/thibault_new/icm_seminar/young_dai.png", width = 5, height = 5, dpi = 300, bg = "white")

#DAI measures at final day of DSS
young_final_DSI_plot <- dssDsiFinalDay(young_dss_followup)
young_final_DSI_plot

# Mice dissection data
dissec_young <- read.csv("young48_dss_dissection.csv", sep = ";", header = TRUE)
dissec_young <- dissectionDataManipulation(dissec_young,4)
dissec_young$gg_group <- factor(dissec_young$gg_group, levels = c("50 water","500 water","50 dss","500 dss"))
dissec_young$colon_length_nrm <- dissec_young$colon_length/dissec_young$body_weight

#boxplot for body weight
young_dissec_bw <- dissecBoxplot(dissec_young,"body") 
young_dissec_bw

#boxplot for std liver weight
young_dissec_lvr <- dissecBoxplot(dissec_young,"liver") 
young_dissec_lvr

# boxplot for std spleen weight
young_dissec_spln <- dissecBoxplot(dissec_young,"spleen") 
young_dissec_spln
# Statistics for spleen measurements
df <- young_dissec_spln$data 
dss50 <- df[df$gg_group == "50 dss",]
dss500 <- df[df$gg_group == "500 dss",]
water50 <- df[df$gg_group == "50 water",]
water500 <- df[df$gg_group == "500 water",]

# Shapiro-Wilk test for normality
print(shapiro.test(dss50[["std_spleen_weigth"]]))
print(shapiro.test(dss500[["std_spleen_weigth"]]))
print(shapiro.test(water50[["std_spleen_weigth"]]))
print(shapiro.test(water500[["std_spleen_weigth"]])) # Conclusion => most of the data is not normally distributed

df$log_spleen <- log(df$std_spleen_weigth)
print(shapiro.test(dss50[["log_spleen"]]))
hist(dss50[["log_spleen"]], breaks = 15)
dss500 <- dss500[-9,] # Remove super high spleen weight outlier
print(shapiro.test(dss500[["log_spleen"]]))
hist(dss500[["log_spleen"]], breaks = 15)
print(shapiro.test(water50[["log_spleen"]]))
print(shapiro.test(water500[["log_spleen"]]))

df <- df[df$id != "10966B", ]

# Levene's test for homogeneity of variance
print(leveneTest(df[["log_spleen"]] ~ gg_group, data = df))

# Perform anova
result <- aov(df[["log_spleen"]] ~ treatment * diet, data = df)
model <- lm(log_spleen ~ treatment * diet, data = df)
result <- anova(model)
print(result)
library(emmeans)
# Compute estimated marginal means (EMMs)
emm <- emmeans(model, pairwise ~ treatment * diet)
# Display results
summary(emm)
print(TukeyHSD(result))

result <- oneway.test(log_spleen ~ gg_group, data = df, var.equal = FALSE)
print(result)
library(rstatix)
post_hoc <- df %>% games_howell_test(log_spleen ~ gg_group)
print(post_hoc)

#boxplot for colon length (non std)
young_dissec_cln <- dissecBoxplot(dissec_young,"colon") 
young_dissec_cln

# boxplot for colon length normalized for body weight
young_dissec_cln <- dissecBoxplot(dissec_young,"colon_nrm") 
young_dissec_cln



# All graphs regarding adult mice colitis experiment 
setwd("experiments/finished exp/adults-all-exp/")

# Mice DAI evolution
adult_dss_followup <- read.csv("combined_adult_dss_followup.csv", header = TRUE, sep = ";")
adult_dss_followup <- dssFollowupManipulation(df = adult_dss_followup,groupInfoCols = 4,dateStart = "2024-05-29",nbrDays = 5, negativeOnly = TRUE)

# creating scatter plot with the four different treatments (diet combined with dss or control)
adult_dssflwup_plot <- dssDiseaseIndexPlot(adult_dss_followup)
adult_dssflwup_plot
ggsave(plot = adult_dssflwup_plot, filename = "~/Documents/CHUM_git/figures/thibault_new/icm_seminar/adult_dai.png", width = 5, height = 5, dpi = 300, bg = "white")

#DAI measures at final day of DSS
adult_final_DSI_plot <- dssDsiFinalDay(adult_dss_followup)
adult_final_DSI_plot

# Mice dissection data
dissec_adult <- read_excel("dss-groups-dissection.xlsx")
dissec_adult <- dissec_adult[dissec_adult$treatment != "abx",] # Remove abx mice
dissec_adult <- dissec_adult[-12,] # Remove dead mouse
dissec_adult[dissec_adult$ID == "36697B", 'liver weight'] <- 1 # For 36697B we are missing liver and spleen measures, replace by 1 and then exlcude if from analysis for spleen and liver measures
dissec_adult[dissec_adult$ID == "36697B", 'spleen weight'] <- 1 
colnames(dissec_adult)[6:9] <- str_replace(colnames(dissec_adult)[6:9], pattern = " ", replacement = "_") # Replace spaces by underscore for variables colnames
dissec_adult <- dissectionDataManipulation(dissec_adult, groupInfoCols = 4)
dissec_adult$gg_group <- factor(dissec_adult$gg_group, levels = c("50 water","500 water","50 dss","500 dss"))
dissec_adult$colon_length_nrm <- dissec_adult$colon_length/dissec_adult$body_weight

#boxplot for body weight
adult_dissec_bw <- dissecBoxplot(dissec_adult,"body") 
adult_dissec_bw

sub_df <- dissec_adult[!dissec_adult$ID== "36697B",] # Sub df removing the mouse for which they were no liver and spleen measures

#boxplot for std liver weight
adult_dissec_lvr <- dissecBoxplot(sub_df,"liver") 
adult_dissec_lvr

# boxplot for std spleen weight
adult_dissec_spln <- dissecBoxplot(sub_df,"spleen") 
adult_dissec_spln

# Statistics for spleen measurements
df <- adult_dissec_spln$data 
dss50 <- df[df$gg_group == "50 dss",]
dss500 <- df[df$gg_group == "500 dss",]
water50 <- df[df$gg_group == "50 water",]
water500 <- df[df$gg_group == "500 water",]

# Shapiro-Wilk test for normality
print(shapiro.test(dss50[["std_spleen_weigth"]]))
print(shapiro.test(dss500[["std_spleen_weigth"]]))
print(shapiro.test(water50[["std_spleen_weigth"]]))
print(shapiro.test(water500[["std_spleen_weigth"]])) # Conclusion => most of the data is not normally distributed

df$log_spleen <- log(df$std_spleen_weigth)
print(shapiro.test(dss50[["log_spleen"]]))
hist(dss50[["log_spleen"]], breaks = 15)
dss500 <- dss500[-9,] # Remove super high spleen weight outlier
print(shapiro.test(dss500[["log_spleen"]]))
hist(dss500[["log_spleen"]], breaks = 15)
print(shapiro.test(water50[["log_spleen"]]))
print(shapiro.test(water500[["log_spleen"]]))

df <- df[df$id != "10966B", ]

# Levene's test for homogeneity of variance
print(leveneTest(df[["log_spleen"]] ~ gg_group, data = df))

# Perform anova
result <- aov(df[["log_spleen"]] ~ treatment * diet, data = df)
model <- lm(log_spleen ~ treatment * diet, data = df)
result <- anova(model)
print(result)
library(emmeans)
# Compute estimated marginal means (EMMs)
emm <- emmeans(model, pairwise ~ treatment * diet)
# Display results
summary(emm)
print(TukeyHSD(result))

result <- oneway.test(log_spleen ~ gg_group, data = df, var.equal = FALSE)
print(result)
library(rstatix)
post_hoc <- df %>% games_howell_test(log_spleen ~ gg_group)
print(post_hoc)

#boxplot for colon length (non std)
adult_dissec_cln <- dissecBoxplot(dissec_adult,"colon") 
adult_dissec_cln

# boxplot for colon length normalized for body weight
adult_dissec_cln <- dissecBoxplot(dissec_adult,"colon_nrm") 
adult_dissec_cln








# All graphs regarding young mice abx experiment 
setwd("experiments/ongoing exp/young-abx-exp6/")

# Mice dissection data
dissec_young <- read_excel("dissection.xlsx")
dissec_young <- dissec_young[-17,] # Remove empty rows
colnames(dissec_young)[5:8] <- str_replace(colnames(dissec_young)[5:8], pattern = " ", replacement = "_") # Replace spaces by underscore for variables colnames
dissec_young <- dissectionDataManipulation(dissec_young,groupInfoCols = 3)
dissec_young$gg_group <- factor(dissec_young$gg_group, levels = c("50 water","500 water","50 abx","500 abx"))
dissec_young$colon_length_nrm <- dissec_young$colon_length/dissec_young$body_weight

#boxplot for body weight
young_dissec_bw <- dissecBoxplot(dissec_young,"body") 
young_dissec_bw

#boxplot for std liver weight
young_dissec_lvr <- dissecBoxplot(dissec_young,"liver") 
young_dissec_lvr

# boxplot for std spleen weight
young_dissec_spln <- dissecBoxplot(dissec_young,"spleen") 
young_dissec_spln
# Statistics for spleen measurements
df <- young_dissec_spln$data 
dss50 <- df[df$gg_group == "50 dss",]
dss500 <- df[df$gg_group == "500 dss",]
water50 <- df[df$gg_group == "50 water",]
water500 <- df[df$gg_group == "500 water",]

# Shapiro-Wilk test for normality
print(shapiro.test(dss50[["std_spleen_weigth"]]))
print(shapiro.test(dss500[["std_spleen_weigth"]]))
print(shapiro.test(water50[["std_spleen_weigth"]]))
print(shapiro.test(water500[["std_spleen_weigth"]])) # Conclusion => most of the data is not normally distributed

df$log_spleen <- log(df$std_spleen_weigth)
print(shapiro.test(dss50[["log_spleen"]]))
hist(dss50[["log_spleen"]], breaks = 15)
dss500 <- dss500[-9,] # Remove super high spleen weight outlier
print(shapiro.test(dss500[["log_spleen"]]))
hist(dss500[["log_spleen"]], breaks = 15)
print(shapiro.test(water50[["log_spleen"]]))
print(shapiro.test(water500[["log_spleen"]]))

df <- df[df$id != "10966B", ]

# Levene's test for homogeneity of variance
print(leveneTest(df[["log_spleen"]] ~ gg_group, data = df))

# Perform anova
result <- aov(df[["log_spleen"]] ~ treatment * diet, data = df)
model <- lm(log_spleen ~ treatment * diet, data = df)
result <- anova(model)
print(result)
library(emmeans)
# Compute estimated marginal means (EMMs)
emm <- emmeans(model, pairwise ~ treatment * diet)
# Display results
summary(emm)
print(TukeyHSD(result))

result <- oneway.test(log_spleen ~ gg_group, data = df, var.equal = FALSE)
print(result)
library(rstatix)
post_hoc <- df %>% games_howell_test(log_spleen ~ gg_group)
print(post_hoc)

#boxplot for colon length (non std)
young_dissec_cln <- dissecBoxplot(dissec_young,"colon") 
young_dissec_cln

# boxplot for colon length normalized for body weight
young_dissec_cln <- dissecBoxplot(dissec_young,"colon_nrm") 
young_dissec_cln

#Iron content measurements
#For young48 dss experiment
{
setwd("experiments/finished exp/young-DSS-exp3/")


#T35 ferrozine assay for stools
df <- read_xlsx("young48_dss_ferrozine_t35.xlsx")
df <- df[,c(2,14:16)]
colnames(df) <- df[2,]
colnames(df)[1:2] <- c("id","iron_concentration")
df <- df[-c(1,2,27),]
df$gg_group <- paste(df$treatment, "+", df$diet, sep = "")
df$gg_group <- factor(df$gg_group, levels = c("water+50","dss+50","water+500","dss+500"))
df$gg_group <- df$diet
df$iron_concentration <- as.numeric(df$iron_concentration)

p = ironBoxplot(df, "iron_concentration", display_significance_bars = F, title = "Iron concentration in stools at day 35", y_axis_title = "yg of iron per g of stools", custom_colors = c("blue","red"), path = "D:/CHUM_git/figures/iron_measures")
p <- p+
  scale_x_discrete(labels = c("50 ppm","500 ppm"))+
  labs(y = "µg iron per g of stool", title = "")+
  guides(color = "none")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 10))
ggsave(filename = "D:/icm_poster_2025/iron_stool_8w.png", height =3, width =2.3 , dpi = 300, bg = "white")


#Fit a ANOVA model for the current time point
anova <- aov(iron_concentration ~ diet * treatment, data = df)
summary(anova)
# Perform Tukey's HSD test and store the results in the list
results <- TukeyHSD(anova)
print(results)

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

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


p = ironBoxplot(df, "total_iron", display_significance_bars = F, title = "Total liver iron at final day", y_axis_title = "yg of iron", custom_colors = c("blue","red"), path = "D:/CHUM_git/figures/iron_measures")
p <- p+
  scale_x_discrete(labels = c("50 ppm\nCtrl","50 ppm\nDSS","500 ppm\nCtrl","500 ppm\nDSS"))+
  labs(y = "µg iron per g of liver", title = "")+
  guides(color = "none")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 9))
ggsave(filename = "D:/icm_poster_2025/iron_liver.png", height =3, width =3 , dpi = 300, bg = "white")

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

p = ironBoxplot(df, "total_iron", display_significance_bars = F, title = "Total spleen iron at final day", y_axis_title = "yg of iron", custom_colors = c("blue","red"), path = "D:/CHUM_git/figures/iron_measures")
p+
  scale_x_discrete(labels = c("50 ppm\nCtrl","50 ppm\nDSS","500 ppm\nCtrl","500 ppm\nDSS"))+
  labs(y = "µg iron per g of spleen", title = "")+
  guides(color = "none")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 9))
ggsave(filename = "D:/icm_poster_2025/iron_spleen.png", height =3, width =3 , dpi = 300, bg = "white")


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
df <- df[-c(1,26,48),]
df$gg_group <- paste(df$treatment, "+", df$diet, sep = "")
df$gg_group <- factor(df$gg_group, levels = c("water+50","dss+50","water+500","dss+500"))
df$iron_concentration <- as.numeric(df$iron_concentration)

a = readxl::read_excel("")
setwd("../")
df




p = ironBoxplot(df, "iron_concentration", display_significance_bars = F, title = "Iron concentration in stools at final day", y_axis_title = "yg of iron per g of stools", custom_colors = c("blue","red"), path = "D:/CHUM_git/figures/iron_measures")
p
p+
  scale_x_discrete(labels = c("50 ppm\nCtrl","50 ppm\nDSS","500 ppm\nCtrl","500 ppm\nDSS"))+
  labs(y = "µg iron per g of liver", title = "")+
  guides(color = "none")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 9))
ggsave(filename = "D:/icm_poster_2025/iron_liver.png", height =3, width =3 , dpi = 300, bg = "white")
#Fit a ANOVA model for the current time point
anova <- aov(iron_concentration ~ diet * treatment, data = df)
summary(anova)
# Perform Tukey's HSD test and store the results in the list
results <- TukeyHSD(anova)
print(results)
}

#Iron content measurements
#For young48 abx experiment
{
setwd("../abx48_ferrozine/")

# Load metadata to retrieve group information
meta <- read_xlsx("abx48_dissection.xlsx")
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
df$gg_group <- df$diet

p = ironBoxplot(df, "iron_concentration", display_significance_bars = F, title = "Iron concentration in stools at day 35", y_axis_title = "µg Fe/g of stool", custom_colors = c("blue","red"), path = "D:/CHUM_git/figures/iron_measures")
p <- p+
  scale_x_discrete(labels = c("50 ppm","500 ppm"))+
  labs(y = "µg Fe/g of stool", title = "")+
  ylim(0,NA)+
  guides(color = "none")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 10))
p
ggsave(filename = "abx48_iron_stool_8w.png", height =3, width =2.3 , dpi = 300, bg = "white")

# Test normality and homoscedasticity
# Levene's Test for homogeneity of variance
leveneTest(iron_concentration ~ gg_group, data = df)

# Shapiro test per group for normality assumptions
by(df$iron_concentration, df$gg_group, shapiro.test)

# Wilcoxon rank-sum test (Mann–Whitney U) (equivalent of t test for non normal and heteroscedastic data)
wilcox.test(df$iron_concentration[df$gg_group == "50"], df$iron_concentration[df$gg_group == "500"], paired = FALSE)

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

p = ironBoxplot(df, "iron_concentration", display_significance_bars = F, title = "Iron concentration in stools at final day", y_axis_title = "µg Fe/g of stool", custom_colors = c("blue","red", "darkblue", "darkred"))
p <- p+
  scale_x_discrete(labels = c("50 ppm\nctrl","500 ppm\nctrl","50 ppm\nAbx","500 ppm\nAbx"))+
  labs(y = "µg Fe/g of stool", title = "")+
  ylim(0,NA)+
  guides(color = "none")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 10))
p
ggsave(filename = "abx48_iron_stool_tf.png", height =3, width =3.5 , dpi = 300, bg = "white")

# Test normality and homoscedasticity
# Levene's Test for homogeneity of variance
leveneTest(iron_concentration ~ gg_group, data = df)

# Shapiro test per group for normality assumptions
by(df$iron_concentration, df$gg_group, shapiro.test)

df$log_iron_c <- log(df$iron_concentration) # Log transform because one group is not normally distributed 

# Levene's Test for homogeneity of variance
leveneTest(log_iron_c ~ gg_group, data = df)

# Shapiro test per group for normality assumptions
by(df$log_iron_c, df$gg_group, shapiro.test)

#Fit a ANOVA model for the current time point
anova <- aov(log_iron_c ~ diet * treatment, data = df)
summary(anova)
# Perform Tukey's HSD test and store the results in the list
results <- TukeyHSD(anova)
print(results)

p = ironBoxplot(df, "log_iron_c", display_significance_bars = F, title = "Iron concentration in stools at final day", y_axis_title = "log(µg Fe/g of stool)", custom_colors = c("blue","red", "darkblue", "darkred"))
p <- p+
  scale_x_discrete(labels = c("50 ppm\nctrl","500 ppm\nctrl","50 ppm\nAbx","500 ppm\nAbx"))+
  labs(y = "log(µg Fe/g of stool)", title = "")+
  ylim(0,NA)+
  guides(color = "none")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 10))
p
ggsave(filename = "abx48_iron_stool_tf_log_transformed.png", height =3, width =3.5 , dpi = 300, bg = "white")
}

