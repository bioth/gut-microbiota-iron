#loading libraries
library("tidyverse") #loading bunch of packages
library("ggplot2") #come on, everyone knows what it is used for
library("dplyr") #arranging and manipulating data easily


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
  stat_summary(aes(group = group),fun = "mean", geom = "point", shape = 18, size =3)+
  stat_summary(fun = "mean",geom = "line",aes(group = group, linetype = ifelse(grepl("DSS", group), "Solid", "Dashed")),size = 1)+
  stat_summary(fun ="mean")+
  labs(title = "Adult mice exposed to iron diets and later DSS, body weight evolution",
       x = "Date",
       y = "Weight (g)")+
  #scale_color_manual(values = c("50 ppm FeSO4 + DSS" = "#96E1E3", "500 ppm FeSO4 + DSS" = "#FE918B",
                                #"50 ppm FeSO4 + water" = "#96E1E3", "500 ppm FeSO4 + water" = "#FE918B")) +
  theme_minimal()+
  ylim(15,25)+
  guides(linetype = "none")

ggsave("your_plot.png", width = 7, height = 4, dpi = 300)
  


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
#this graph has a disease index score in the y column
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


#this graph has the weight measurement during the 6 days of DSS for the y column
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
library(geepack) #library for loading GEE tests
library("lme4") #library for loading ANOVA

model <- geeglm(weight ~ Diet...2 * Treatment...3, id = ID...4, data = combined_dfa)
helsummary(model) #NOT WORKING? Need to figure that out

#Renaming the cols names that had been modified prior to that so that it's more intuitive
combined_dfa$date <- as.factor(combined_dfa$date)
colnames(combined_dfa)[1:5] <- c("Cage","Diet","Treatment","ID","Date")

#defining a 3-way repeated ANOVA (3 independant variables including date, treatment and diet)
model <- aov(weight ~ Diet * Treatment * Date + Error(ID/Date), data = combined_dfa)
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

#box plot for spleen
spleen_box_plot <- dissec %>%
  ggplot(aes(x = gg_group, y = std_spleen_weigth, color = gg_group))+
  geom_boxplot(width = 0.5, outlier.shape = NA) +  # Customize box width and hide outliers
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Normalized spleen weight measures at final day of the experiment",
       x = "Treatment",
       y = "Normalized weight(spleen weight/body weight)")+
  ylim(0,0.01)+
  theme_minimal()

#box plot for liver
liver_box_plot <- dissec %>%
  ggplot(aes(x = gg_group, y = std_liver_weigth, color = gg_group))+
  geom_boxplot(width = 0.5, outlier.shape = NA) +  # Customize box width and hide outliers
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Normalized liver weight measures at final day of the experiment",
       x = "Treatment",
       y = "Normalized weight(liver weight/body weight)")+
  theme_minimal()

#box plot of colon length
colon_box_plot <- dissec %>%
  ggplot(aes(x = gg_group, y = colon.length, color = gg_group))+
  geom_boxplot(width = 0.5, outlier.shape = NA) +  # Customize box width and hide outliers
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Colon length measures at final day of the experiment",
       x = "Treatment",
       y = "Colon length (cm)")+
  theme_minimal()

#box plot for standardized colon length
nrm_colon_box_plot <- dissec %>%
  ggplot(aes(x = gg_group, y = std_colon_len, color = gg_group))+
  geom_boxplot(width = 0.5, outlier.shape = NA) +  # Customize box width and hide outliers
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Normalized colon length measures at final day of the experiment",
       x = "Treatment",
       y = "Normalized colon length (divided by body weight)")+
  theme_minimal()

#box plot for body weight
body_weight_plot <- dissec %>%
  ggplot(aes(x = gg_group, y = body_weight, color = gg_group))+
  geom_boxplot(width = 0.5, outlier.shape = NA) +  # Customize box width and hide outliers
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Body weight measures at final day of the experiment",
       x = "Treatment",
       y = "Body weight (g)")+
  theme_minimal()



#saving figures
setwd("D:/CHUM_git/figures/")
ggsave(plot = spleen_box_plot,"spleen_weight.png", width = 7, height = 4, dpi = 300, bg = "white")
ggsave(plot = liver_box_plot,"liver_weight.png", width = 7, height = 4, dpi = 300, bg = "white")
ggsave(plot = colon_box_plot,"colon_length.png", width = 7, height = 4, dpi = 300, bg = "white")
ggsave(plot = nrm_colon_box_plot,"normalized_colon_length.png", width = 7, height = 4, dpi = 300, bg = "white")
ggsave(plot = colon_box_plot,"body_wweight.png", width = 7, height = 4, dpi = 300, bg = "white")


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



# Perform two-way ANOVA
result <- aov(colon.length ~ Group + Diet + Group:Diet, data = dissec)

# Print the summary of the ANOVA
summary(result)


# Perform two-way ANOVA
result <- aov(std_colon_len ~ Group + Diet + Group:Diet, data = dissec)

# Print the summary of the ANOVA
summary(result)

# Perform two-way ANOVA
result <- aov(std_liver_weigth ~ Group + Diet + Group:Diet, data = dissec)

# Print the summary of the ANOVA
summary(result)

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
