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

adult_weight <- weightDataManipulation(adult_weight,4)

#replacing abnormal values
adult_weight$weight[4] <- 24.2
adult_weight$weight[28] <- 23.1

#using this enables to verify if date variable is of type "date" and weight of type numeric
str(adult_weight)

#creating scatter plot with the four different treatments (diet combined with dss or control)
adult_weight_plot <- weightPlot(adult_weight)
adult_weight_plot

#saving scatter plot
setwd("../figures/adult36")
ggsave("adults_weight.png", width = 9, height = 6, dpi = 300, bg = "white")
  
#statistics
#here we will attest for differences in weight, first need to assess normality and homoscedasticity
#testing for homoscedasticity

leveneTest(weight ~ treatment * diet * time_numeric, data = weight_measure_file)


shapiro.test(residuals(model))

qqnorm(residuals(model))
qqline(residuals(model))


#performing a repeated ANOVA

anova_result <- aov(weight ~ treatment * diet * time_numeric + Error(ID/time_numeric),
                    data = weight_measure_file)



summary(anova_result)



#YOUNG MICE
#Loading weight measure file for young mice
setwd("../young-DSS-exp2/")
young_weight <- read.csv("young32dss_weight.csv", header = TRUE, sep = ";")

#data manipulation
young_weight <- weightDataManipulation(young_weight,4)

#creating scatter plot with the four different treatments (diet combined with dss or control)
young_weight_plot <- weightPlot(young_weight)
young_weight_plot

#saving scatter plot
setwd("../figures/young32/")
ggsave("young32_weight.png", width = 10, height = 6, dpi = 300, bg = "white")

#performing a repeated ANOVA
anova_result <- aov(weight ~ treatment * diet * time_numeric + Error(id/time_numeric),
                    data = young_weight)
summary(anova_result)



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
setwd("../figures/young32/")
ggsave("young32_dIndex.png", width = 9, height = 5, dpi = 300, bg = "white")

###going for the statistical measurements
model <- geeglm(index ~ diet * treatment + time_numeric, id = id, data = combined_dfa)
summary(model) #NOT WORKING? Need to figure that out


#performing a repeated ANOVA
anova_result <- aov(index ~ treatment * diet * time_numeric + Error(id/time_numeric),
                    data = young_dss_followup)
summary(anova_result)



#defining a 3-way repeated ANOVA (3 independant variables including date, treatment and diet)
model <- Anova(aov(weight ~ Diet * Treatment * Date + Error(ID/Date), data = combined_dfa))
summary(model)

#Though you gotta verify that the residuals are normally distributed and are homoscedastic
#To do that you need to test the residuals separately (NEED TO BE VERIFIED THOUGH)

# Assuming 'model' is your mixed-effects model
residuals_ID <- residuals(model$ID)
residuals_ID_date <- residuals(model$`ID:Date`)

# Q-Q plot for normality - ID residuals
qqnorm(residuals_ID)
qqline(residuals_ID)
shapiro.test(residuals_ID)
#assumed to be normally distributed

# Q-Q plot for normality - ID:date residuals
qqnorm(residuals_ID_date)
qqline(residuals_ID_date)
shapiro.test(residuals_ID_date)
#assumed to be normally distributed

# Levene's test for homogeneity of variances - ID residuals
leveneTest(residuals_ID ~ Diet * Treatment * Date, data = combined_dfa)
#NOT WORKING

# Levene's test for homogeneity of variances - ID:date residuals
leveneTest(residuals_ID_date ~ Diet * Treatment * date, data = combined_dfa)
#NOT WORKING








###Final dissection data###
#loading the data
setwd("../../adult-DSS-exp/")
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
ggsave(plot = adult_dissec_cln,"colon_length.png", width = 9, height = 5, dpi = 300, bg = "white")
ggsave(plot = adult_dissec_bw,"body_wweight.png", width = 9, height = 5, dpi = 300, bg = "white")


#loading the data
setwd("../young-DSS-exp2//")
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
ggsave(plot = young_dissec_cln,"colon_length.png", width = 9, height = 5, dpi = 300, bg = "white")
ggsave(plot = young_dissec_bw,"body_weight.png", width = 9, height = 5, dpi = 300, bg = "white")













#trying stuff with significance bars
# Extract p-value for the Group factor from ANOVA result
group_p_value <- summary(result)[[1]]$`Pr(>F)`[1]

# Add significance bar if p-value is below the threshold (e.g., 0.05)
if(group_p_value < 0.05) {
  liver_box_plot + 
    geom_signif(comparisons = list(c("50 ppm FeSO4 + DSS", "500 ppm FeSO4 + water")),
                annotations = c(pval = round(group_p_value, 3)), 
                map_signif_level = TRUE, textsize = 5, vjust = -0.5)
} else {
  liver_box_plot
}


dissec <- dissec_young

#messing up with stats
# Subsetting the dataframe to include only the DSS groups
dss_data <- dissec[dissec$gg_group %in% c("50 ppm FeSO4 + DSS", "500 ppm FeSO4 + DSS"), ]


# Shapiro-Wilk test for normality
shapiro.test(dss_data)
shapiro.test(group2)

# Levene's test for homogeneity of variance
library(car)
leveneTest(group1, group2)

# Perform t-test
result_colon <- t.test(colon_length ~ gg_group, data = dss_data)
result_spleen <- t.test(std_spleen_weigth ~ gg_group, data = dss_data)
result_liver <- t.test(std_liver_weigth ~ gg_group, data = dss_data)
result_nrm_colon <- t.test(std_colon_len ~ gg_group, data = dss_data)


# Print the result
print(result_colon)
print(result_spleen)
print(result_liver)
print(result_nrm_colon)


#this one works I think
anova_result <- aov(index ~ treatment * diet * time_numeric + Error(id/time_numeric),
                    data = combined_dfa)

# Perform two-way ANOVA
result <- aov(colon.length ~ Group + Diet + Group:Diet, data = dissec)

# Print the summary of the ANOVA
summary(result)


# Perform two-way ANOVA
result <- aov(std_colon_len ~ Group + Diet + Group:Diet, data = dissec)

# Print the summary of the ANOVA
summary(result)

# Perform two-way ANOVA
result_liver <- aov(std_liver_weigth ~ Group + Diet + Group:Diet, data = dissec)

# Print the summary of the ANOVA
summary(result_liver)

# Perform two-way ANOVA
result <- aov(std_spleen_weigth ~ Group + Diet + Group:Diet, data = dissec)

# Print the summary of the ANOVA
summary(result)

# Perform two-way ANOVA
result <- aov(body_weight ~ Group + Diet + Group:Diet, data = dissec)

# Print the summary of the ANOVA
summary(result)


#ancova test (taking into account age at start of exp)
# Fit ANCOVA model
model <- lm(colon.length~ Diet + Group + Diet:Group + Age.start..d. + body_weight, data = dissec)

# Print the summary of the model
summary(model)
