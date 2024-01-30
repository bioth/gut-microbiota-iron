#loading libraries


library("tidyverse")
library("ggplot2")
library("dplyr")


###WEIGHT MEASURES FILE DATA
#loading the weight_measure file from github (using function in fetchGHdata.R)
weight_measure_file <- gh_data

#if the function does not work
setwd("D:/CHUM_git/gut-microbiota-iron/")
weight_measure_file <- read.csv("adult_dss_weight_measurement.csv", header = TRUE, sep = "\t")

#removing the 3 dead mice
weight_measure_file = weight_measure_file[-c(16,27,36,37),]

#using pivot longer from dplyr package to reorganize the data
weight_measure_file = pivot_longer(weight_measure_file, c(10:length(weight_measure_file)), cols_vary = "slowest", names_to = "date", values_to = "weight")

#modifying the date format so that it is recognizable by R
for(i in c(1:length(weight_measure_file$date))){
  weight_measure_file$date[i] <- gsub("\\.", "-", weight_measure_file$date[i])
  weight_measure_file$date[i] <- substring(weight_measure_file$date[i], first = 3)
  weight_measure_file$date[i] <- substr( weight_measure_file$date[i], start = 1, stop = nchar( weight_measure_file$date[i]) - 1)
  
}

#tranforming in date format
weight_measure_file$date <- as.Date(weight_measure_file$date)

#setting the diet column data into a string format so that it can be used into ggplot
weight_measure_file$Diet..Fe.ppm. = as.character(weight_measure_file$Diet..Fe.ppm.)

#setting the weight column data into numeric values
for(i in c(1:length(weight_measure_file$weight))){
  weight_measure_file$weight[i] <- gsub("\\,", ".", weight_measure_file$weight[i])
  
}
weight_measure_file$weight <- as.numeric(weight_measure_file$weight)

#using this enables to verify if date variable is of type "date" and weight of type numeric
str(weight_measure_file)

#replacing abnormal values
weight_measure_file$weight[5] = 24.2
weight_measure_file$weight[29] = 23.1


#creating 4 groups for easier graphic interpretations
for(i in c(1:length(weight_measure_file$arbitrary.number))){
  if(any(weight_measure_file$arbitrary.number[i] %in% c(1,3,5))){
    weight_measure_file$group[i] <- "50 ppm FeSO4 + DSS"
  }
  if(any(weight_measure_file$arbitrary.number[i] %in% c(2,4,6))){
    weight_measure_file$group[i] <- "50 ppm FeSO4 + water"
  }
  if(any(weight_measure_file$arbitrary.number[i] %in% c(7,9,11))){
    weight_measure_file$group[i] <- "500 ppm FeSO4 + DSS"
  }
  if(any(weight_measure_file$arbitrary.number[i] %in% c(8,10,12))){
    weight_measure_file$group[i] <- "500 ppm FeSO4 + water"
  }
}

#creating scatter plot with the four different treatments (diet combined with dss or control)
data <- as.data.frame(weight_measure_file)
data %>%
  ggplot(aes(x = date, y = weight, color = group))+
  geom_point(aes(color = group), size = 2)+
  stat_summary(fun ="mean")+
  labs(title = "Body weight through time",
       x = "Day",
       y = "Weight (g)")+
  theme_minimal()+
  ylim(1,30)

  


###DSS FOLLOW UP SHEET DATA
#loading the dss followup sheet from github
github_repo <- "gut-microbiota-iron"
file <- "adult_mice_DSS1_followup.csv"
dss_followup <- read.table(text = fetchGHdata(github_repo,file), header = TRUE, sep = "\t")

#if the function is not working
setwd("D:/CHUM_git/gut-microbiota-iron/")
dss_followup <- read.csv("adult_mice_DSS1_followup.csv", header = TRUE, sep = "\t")

#replacing "bleeding" for "hemoccult" for one of the columns
dss_followup$Day.0.1[1] <- "hemoccult"


#creating 3 different dataframes for each metric (weight, hemoccult, stool consistency)

#defining a function that takes a dataframe, the string we are looking for in the first row
#and the columns we wanna ignore and add at the end
colStringSelect <- function(df,string,cols_ignored){
  selected_cols <- grep(string, df[1,(cols_ignored+1):ncol(df)])
  selected_cols <- selected_cols + cols_ignored
  df_select <- df[, selected_cols, drop = FALSE]
  result_df <- cbind(df[,1:cols_ignored],df_select)
  return(as.data.frame(result_df))
}

dss_weight <- colStringSelect(dss_followup,"weight",4)

dss_hemo <- colStringSelect(dss_followup,"hemoccult",4)

dss_consistency <- colStringSelect(dss_followup,"consistency",4)


#replacing the colnames by the dates for the DSS
day_colnames <- as.Date(c("2023-12-04","2023-12-05","2023-12-06","2023-12-07","2023-12-08","2023-12-09"))
colnames(dss_weight)[5:10] <- format(day_colnames, "%Y-%m-%d")
colnames(dss_consistency)[5:10] <- format(day_colnames, "%Y-%m-%d")
colnames(dss_hemo)[5:10] <- format(day_colnames, "%Y-%m-%d")

#removing the first row and the mouse that died before the DSS
colnames(dss_weight)[1:4] <- dss_weight[1,1:4]
dss_weight <- dss_weight[-c(1,28),]

colnames(dss_consistency)[1:4] <- dss_consistency[1,1:4]
dss_consistency <- dss_consistency[-c(1,28),]

colnames(dss_hemo)[1:4] <- dss_hemo[1,1:4]
dss_hemo <- dss_hemo[-c(1,28),]

#using pivot_longer so that the TEMPORAL data can be displayed by ggplot
dss_weight_gg <- pivot_longer(dss_weight, c(5:length(dss_weight)), cols_vary = "slowest", names_to = "date", values_to = "weight")
dss_consistency_gg <- pivot_longer(dss_consistency, c(5:length(dss_consistency)), cols_vary = "slowest", names_to = "date", values_to = "consistency")
dss_hemo_gg <- pivot_longer(dss_hemo, c(5:length(dss_hemo)), cols_vary = "slowest", names_to = "date", values_to = "hemo")

###finally merging the data...###
combined_df <- bind_cols(dss_weight_gg, dss_consistency_gg, dss_hemo_gg)

# Identify duplicated columns
duplicated_cols <- duplicated(names(combined_df)) | duplicated(t(combined_df))

# Remove duplicated columns
combined_df <- combined_df[, !duplicated_cols, drop = FALSE]

#replace NL (normal-loose) by just loose, for the consistency, I was just misinterpreting
for(i in 1:length(combined_df$consistency)){
  if(combined_df$consistency[i]=="NL"){
    combined_df$consistency[i] <- "N"
  }
}

###plotting the data during the DSS
#creating a disease severity index, can go from 0 (no bleeding normal consistency) to 2 (bleeding and loose)
combined_df$index <- 0
combined_dfa <- combined_df
for(i in 1:length(combined_dfa$consistency)){
  if(combined_dfa$consistency[i]=="N"){
    combined_dfa$index[i] <- combined_dfa$index[i]+0
  }else{
    combined_dfa$index[i] <- combined_dfa$index[i]+1
  }
  if(combined_dfa$hemo[i]=="No"){
    combined_dfa$index[i] <- combined_dfa$index[i]+0
  }else{
    combined_dfa$index[i] <- combined_dfa$index[i]+1
  }
}

#creating 4 groups for easier graphic interpretations
for(i in c(1:length(combined_dfa$Cage...1))){
  if(any(combined_dfa$Cage...1[i] %in% c(1,3,5))){
    combined_dfa$group[i] <- "50 ppm FeSO4 + DSS"
  }
  if(any(combined_dfa$Cage...1[i] %in% c(2,4,6))){
    combined_dfa$group[i] <- "50 ppm FeSO4 + water"
  }
  if(any(combined_dfa$Cage...1[i] %in% c(7,9,11))){
    combined_dfa$group[i] <- "500 ppm FeSO4 + DSS"
  }
  if(any(combined_dfa$Cage...1[i] %in% c(8,10,12))){
    combined_dfa$group[i] <- "500 ppm FeSO4 + water"
  }
}

#creating scatter plot with the four different treatments (diet combined with dss or control)
data <- as.data.frame(combined_dfa)
data %>%
  ggplot(aes(x = date...5, y = index, color = group))+
  stat_summary(fun ="mean", geom = "point", shape = 18, size = 3)+
  stat_summary(
    fun = "mean",
    geom = "line",
    aes(group = group),
    size = 1
  ) +
  labs(title = "Disease severity index through time",
       x = "Day",
       y = "Severity")+
  theme_minimal()

#setting the weight column data into numeric values
for(i in c(1:length(combined_dfa$weight))){
  combined_dfa$weight[i] <- gsub("\\,", ".", combined_dfa$weight[i])
  
}
combined_dfa$weight <- as.numeric(combined_dfa$weight)

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


#going for the statistical measurements
library(geepack)

model <- geeglm(weight ~ Diet...2 * Treatment...3, id = ID...4, data = combined_dfa)
helsummary(model)


library("lme4")
combined_dfa$date <- as.factor(combined_dfa$date)
colnames(combined_dfa)[1:5] <- c("Cage","Diet","Treatment","ID","Date")

model <- aov(weight ~ Diet * Treatment * Date + Error(ID/Date), data = combined_dfa)

summary(model)

# Residuals from your model
 
residuals(model)

# Q-Q plot
qqnorm(residuals)
qqline(residuals)

# Shapiro-Wilk test for normality
shapiro.test(residuals)

# Levene's test for homogeneity of variances
library(car)
leveneTest(residuals ~ Diet * Treatment * Date, data = combined_dfa)


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

# Levene's test for homogeneity of variances - ID:date residuals
leveneTest(residuals_ID_date ~ Diet * Treatment * date, data = combined_dfa)


model <- aov(weight ~ Diet * Treatment * Date + Error(ID/Date), data = combined_dfa)

summary(model)
