#loading libraries
{
  library("tidyverse") #loading bunch of packages
  library("ggplot2") #come on, everyone knows what it is used for
  library("dplyr") #arranging and manipulating data easily
  library("geepack") #library for loading GEE tests
  library("lme4") #library for loading ANOVA
  library("car") #for anova too
  library("ggsignif") #adding significance bars to ggplots
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
adult_weight <- read.csv("adult_dss_weight_measurement.csv", header = TRUE, sep = ";")

#removing useless cols
adult_weight <- adult_weight[,-c(1,3:5,7)]

#modifying wrong colnames (colnames with dates managed by function changeDateCols)
colnames(adult_weight)[1:4] <- c("id","diet","cage","treatment")

adult_weight <- weightDataManipulation(adult_weight,4)

#replacing abnormal values
weight_measure_file$weight[4] <- 24.2
weight_measure_file$weight[28] <- 23.1

#using this enables to verify if date variable is of type "date" and weight of type numeric
str(weight_measure_file)

#creating scatter plot with the four different treatments (diet combined with dss or control)
adult_weight_plot <- weightPlot(adult_weight)
adult_weight_plot

#saving scatter plot
setwd("D:/CHUM_git/figures/")
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
setwd("../young-DSS-exp2//")
young_weight <- read.csv("weight_measurement.csv", header = TRUE, sep = ";")

#removing useless cols
young_weight <- young_weight[,-c(5,6)]

#modifying wrong colnames (colnames with dates managed by function changeDateCols)
colnames(young_weight)[1:4] <- c("diet","treatment","cage","id")

#data manipulation
young_weight <- weightDataManipulation(young_weight,4)

#creating scatter plot with the four different treatments (diet combined with dss or control)
young_weight_plot <- weightPlot(young_weight)
young_weight_plot






###DSS FOLLOW UP SHEET DATA
#loading the dss followup sheet data
setwd("../adult-DSS-exp/")
adult_dss_followup <- read.csv("adult_mice_DSS1_followup.csv", header = TRUE, sep = "\t")
setwd("../young-DSS-exp2/")
young_dss_followup <- read.csv("dss_followup_young_32.csv", header = TRUE, sep = ";")

#loading the functions
setwd("../r scripts/")
source("dataManipFunctions.R")


###LOADING ADULT MICE DATA
#replacing colnames
colnames(adult_dss_followup)[0:4] <- c("cage","diet","treatment","id")

#replacing "bleeding" for "hemoccult" for one of the columns
adult_dss_followup$Day.0.1[1] <- "hemoccult"

adult_dss_followup <- dssFollowupManipulation(df = adult_dss_followup,groupInfoCols = 4,dateStart = "2023-12-04",nbrDays = 5)


###LOADING YOUNG MICE DATA
#getting rid of useless col
young_dss_followup <- young_dss_followup[,-3]

#replacing colnames
colnames(young_dss_followup)[0:4] <- c("cage","id","diet","treatment")

young_dss_followup <- dssFollowupManipulation(df = young_dss_followup,groupInfoCols = 4,dateStart = "2024-02-14",nbrDays = 5)

#creating scatter plot with the four different treatments (diet combined with dss or control)
#this graph has a disease index score in the y column
adult_dssflwup_plot <- dssDiseaseIndexPlot(adult_dss_followup)
adult_dssflwup_plot

young_dssflwup_plot <- dssDiseaseIndexPlot(young_dss_followup)
young_dssflwup_plot

#saving figure
setwd("D:/CHUM_git/figures/")
ggsave("disease_index.png", width = 7, height = 4, dpi = 300, bg = "white")

###going for the statistical measurements
# Convert 'date' to Date format if not already done
combined_dfa$date <- as.Date(combined_dfa$date)

# Choose a reference date
reference_date <- min(combined_dfa$date)

# Convert dates to numeric values representing the number of days elapsed since the reference date
combined_dfa$time_numeric <- as.numeric(combined_dfa$date - reference_date)

model <- geeglm(index ~ diet * treatment + time_numeric, id = id, data = combined_dfa)
summary(model) #NOT WORKING? Need to figure that out


#anova test
anova_result <- aov(index ~ treatment * diet * time_numeric + Error(id/time_numeric),
                          data = combined_dfa)



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
setwd("../adult-DSS-exp/")
dissec_adult <- read.csv("exp1_DSS_dissection.csv", sep = ";", header = TRUE)
setwd("../r scripts/")
source("dataManipFunctions.R")

#getting rid of useless cols
dissec_adult <- dissec_adult[,-(9:10)]

#changing colnames
colnames(dissec_adult)[1:4] <- c("cage","diet","treatment","id")

dissec_adult <- dissectionDataManipulation(dissec_adult,4)

#boxplot for body weight
adult_dissec_bw <- dissecBoxplot(dissec_adult,"body") 
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
setwd("D:/CHUM_git/figures/")
ggsave(plot = adult_dissec_spln,"spleen_weight.png", width = 8, height = 4, dpi = 300, bg = "white")
ggsave(plot = adult_dissec_lvr,"liver_weight.png", width = 9, height = 5, dpi = 300, bg = "white")
ggsave(plot = adult_dissec_cln,"colon_length.png", width = 9, height = 5, dpi = 300, bg = "white")
ggsave(plot = adult_dissec_bw,"body_wweight.png", width = 9, height = 5, dpi = 300, bg = "white")





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
result_colon <- t.test(colon.length ~ gg_group, data = dss_data)
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
