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
weight_measure_file <- read.csv("adult_dss_weight_measurement.csv", header = TRUE, sep = ";")


#check file for checking cols to remove and colnames to change
View(weight_measure_file)

#removing useless cols
weight_measure_file <- weight_measure_file[,-c(1,3:5,7)]

#modifying wrong colnames (colnames with dates managed by function changeDateCols)
colnames(weight_measure_file)[1:4] <- c("id","diet","cage","treatment")

weight_measure_file <- weightDataManipulation(weight_measure_file,4)

#replacing abnormal values
weight_measure_file$weight[4] <- 24.2
weight_measure_file$weight[28] <- 23.1

#using this enables to verify if date variable is of type "date" and weight of type numeric
str(weight_measure_file)

#creating scatter plot with the four different treatments (diet combined with dss or control)
data <- as.data.frame(weight_measure_file)
data %>%
  ggplot(aes(x = date, y = weight, color = diet)) +
  stat_summary(aes(group = group, shape = treatment), fun = "mean", geom = "point", size = 3) +
  stat_summary(fun = "mean", geom = "line", aes(group = group, linetype = ifelse(grepl("DSS", group), "DSS", "Water")), size = 1) +
  labs(title = "Adult mice exposed to iron diets and later DSS, body weight evolution",
       x = "Date",
       y = "Weight (g)",
       color = "Diet") +
  scale_linetype_manual(name = "Treatment", 
                        values = c("DSS" = "dashed", "Water" = "solid")) +
  guides(shape = 'none')+
  theme_minimal() +
  ylim(15, 25)+
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
    axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 12),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
    legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", size = 1)  # Include axis lines  # Include axis bars
  )

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

#check file for checking cols to remove and colnames to change
View(young_weight)

#removing useless cols
young_weight <- young_weight[,-c(5,6)]

#modifying wrong colnames (colnames with dates managed by function changeDateCols)
colnames(young_weight)[1:4] <- c("diet","treatment","cage","id")

young_weight <- weightDataManipulation(young_weight)








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
testPlot <- dssDiseaseIndexPlot(adult_dss_followup)
testPlot

testPlot2 <- dssDiseaseIndexPlot(young_dss_followup)
testPlot2

setwd("D:/CHUM_git/figures/")
ggsave("disease_index.png", width = 7, height = 4, dpi = 300, bg = "white")



#this graph has the weight measurement during the 6 days of DSS for the y column = doesnt work anymore lol
data <- as.data.frame(combined_dfa)
data %>%
  ggplot(aes(x = date...5, y = weight, color = group))+
  stat_summary(aes(group = group),fun = "mean", geom = "point", shape = 18, size =3)+
  stat_summary(
    fun = "mean",
    geom = "line",
    aes(group = group),
    size = 1
  ) +
  labs(title = "Weight measurements for the DSS week",
       x = "Day",
       y = "Weight (g)")+
  theme_minimal()





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
setwd("D:/CHUM_git/gut-microbiota-iron/adult-DSS-exp/")
dissec <- read.csv("exp1_DSS_dissection.csv", sep = ";", header = TRUE)

#removing null lines (corresponding to dead mice)
#empty is an array with the lines numbers where there is nothing
empty = NULL
for (i in 1:length(dissec$Cage)){
  if(dissec$ID[i] == ""){
    empty <- c(empty,i)
  }
}
dissec <- dissec[-c(empty),]



#####box plots for colon, spleen, and liver
#transforming the chr into numeric (measures, because of the ","), creating a function
char_into_num <- function(df,start,stop) {
  for(i in start:stop){
    for(k in 1:nrow(df)){
      df[k,i] <- gsub("\\,", ".", df[k,i])
      
    }
  }
  df[, start:stop] <- lapply(df[, start:stop], as.numeric)
  return(df)
}

dissec <- char_into_num(dissec,5,8)

#dividing these measures by final weight of the mice (normalization)
dissec$std_spleen_weigth <- dissec$spleen.weight/dissec$body_weight
dissec$std_liver_weigth <- dissec$liver.weight/dissec$body_weight
dissec$std_colon_len <- dissec$colon.length/dissec$body_weight

#creating 4 groups for easier graphic interpretations
for(i in 1:nrow(dissec)){
  if(any(dissec$Cage[i] %in% c(1,3,5))){
    dissec$gg_group[i] <- "50 ppm FeSO4 + DSS"
  }
  if(any(dissec$Cage[i] %in% c(2,4,6))){
    dissec$gg_group[i] <- "50 ppm FeSO4 + water"
  }
  if(any(dissec$Cage[i] %in% c(7,9,11))){
    dissec$gg_group[i] <- "500 ppm FeSO4 + DSS"
  }
  if(any(dissec$Cage[i] %in% c(8,10,12))){
    dissec$gg_group[i] <- "500 ppm FeSO4 + water"
  }
}


dissec <- dissec %>%
  mutate(grouping = ifelse(grepl("DSS", gg_group), "DSS", "Water"))

#box plot for spleen
spleen_box_plot <- dissec %>%
  ggplot(aes(x = factor(gg_group, levels = c("50 ppm FeSO4 + DSS", "500 ppm FeSO4 + DSS", "50 ppm FeSO4 + water", "500 ppm FeSO4 + water")), 
             y = std_spleen_weigth, color = gg_group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +  # Customize box width and hide outliers
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Normalized spleen weight measures at final day of the experiment",
       x = "Treatment",
       y = "Normalized weight(spleen weight/body weight)") +
  ylim(0, 0.01) +
  theme_minimal()
spleen_box_plot


spleen_box_plot <- dissec %>%
  ggplot(aes(x = factor(grouping), y = std_spleen_weigth, color = as.character(Diet)))+
  geom_boxplot(width = 0.5, outlier.shape = NA, aes(linetype = ifelse(grepl("DSS", grouping), "DSS", "Water"))) +  # Customize box width and hide outliers
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), alpha = 0.5, aes(shape = grouping)) +
  labs(title = "Normalized spleen weight measures at final day of the experiment",
       x = "Treatment",
       y = "Normalized spleen weight (liver weight/body weight)",
       color = "Diet")+ 
  scale_linetype_manual(name = "Treatment", 
                        values = c("DSS" = "dashed", "Water" = "solid")) +
  guides(shape = 'none')+
  theme_minimal() +  # Use minimal theme
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
    axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 12),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
    legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", size = 1)  # Include axis lines  # Include axis bars
  )+
  ylim(0, 0.01)




#box plot for liver
liver_box_plot <- dissec %>%
  ggplot(aes(x = factor(grouping), y = std_liver_weigth, color = as.character(Diet)))+
  geom_boxplot(width = 0.5, outlier.shape = NA, aes(linetype = ifelse(grepl("DSS", grouping), "DSS", "Water"))) +  # Customize box width and hide outliers
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), alpha = 0.5, aes(shape = grouping)) +
  labs(title = "Normalized liver weight measures at final day of the experiment",
       x = "Treatment",
       y = "Normalized liver weight(liver weight/body weight)",
       color = "Diet")+ 
  scale_linetype_manual(name = "Treatment", 
                        values = c("DSS" = "dashed", "Water" = "solid")) +
  guides(shape = 'none')+
  theme_minimal() +  # Use minimal theme
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
    axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 12),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
    legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", size = 1)  # Include axis lines  # Include axis bars
  )

liver_box_plot

#box plot of colon length
colon_box_plot <- dissec %>%
  ggplot(aes(x = factor(gg_group, levels = c("50 ppm FeSO4 + DSS", "500 ppm FeSO4 + DSS", "50 ppm FeSO4 + water", "500 ppm FeSO4 + water")), y = colon.length, color = gg_group))+
  geom_boxplot(width = 0.5, outlier.shape = NA) +  # Customize box width and hide outliers
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Colon length measures at final day of the experiment",
       x = "Treatment",
       y = "Colon length (cm)")+
  theme_minimal()

colon_box_plot <- dissec %>%
  ggplot(aes(x = factor(grouping), y = colon.length, color = as.character(Diet)))+
  geom_boxplot(width = 0.5, outlier.shape = NA, aes(linetype = ifelse(grepl("DSS", grouping), "DSS", "Water"))) +  # Customize box width and hide outliers
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), alpha = 0.5, aes(shape = grouping)) +
  labs(title = "Colon length measures at final day of the experiment",
       x = "Treatment",
       y = "Colon length (cm)",
       color = "Diet")+ 
  scale_linetype_manual(name = "Treatment", 
                        values = c("DSS" = "dashed", "Water" = "solid")) +
  guides(shape = 'none')+
  theme_minimal() +  # Use minimal theme
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
    axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 12),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
    legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", size = 1)  # Include axis lines  # Include axis bars
  )

#box plot for standardized colon length
nrm_colon_box_plot <- dissec %>%
  ggplot(aes(x = factor(gg_group, levels = c("50 ppm FeSO4 + DSS", "500 ppm FeSO4 + DSS", "50 ppm FeSO4 + water", "500 ppm FeSO4 + water")), y = std_colon_len, color = gg_group))+
  geom_boxplot(width = 0.5, outlier.shape = NA) +  # Customize box width and hide outliers
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Normalized colon length measures at final day of the experiment",
       x = "Treatment",
       y = "Normalized colon length (divided by body weight)")+
  theme_minimal()

#box plot for body weight
body_weight_plot <- dissec %>%
  ggplot(aes(x = factor(gg_group, levels = c("50 ppm FeSO4 + DSS", "500 ppm FeSO4 + DSS", "50 ppm FeSO4 + water", "500 ppm FeSO4 + water")), y = body_weight, color = gg_group))+
  geom_boxplot(width = 0.5, outlier.shape = NA) +  # Customize box width and hide outliers
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Body weight measures at final day of the experiment",
       x = "Treatment",
       y = "Body weight (g)")+
  theme_minimal()





#saving figures
setwd("D:/CHUM_git/figures/")
ggsave(plot = spleen_box_plot,"spleen_weight.png", width = 8, height = 4, dpi = 300, bg = "white")
ggsave(plot = liver_box_plot,"liver_weight.png", width = 9, height = 5, dpi = 300, bg = "white")
ggsave(plot = colon_box_plot,"colon_length.png", width = 9, height = 5, dpi = 300, bg = "white")
ggsave(plot = nrm_colon_box_plot,"normalized_colon_length.png", width = 9, height = 5, dpi = 300, bg = "white")
ggsave(plot = body_weight_plot,"body_wweight.png", width = 9, height = 5, dpi = 300, bg = "white")



min(dissec$colon.length)











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
