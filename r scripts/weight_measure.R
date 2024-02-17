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





#Loading weight measure file
setwd("adult-DSS-exp/")
weight_measure_file <- read.csv("adult_dss_weight_measurement.csv", header = TRUE, sep = ";")


#Function removing rows with empty values
emptyRow <- function(df){
  rows_to_remove <- c()  # Initialize an empty vector to store row indices to remove
  for(i in 1:nrow(df)){
    for(k in 1:ncol(df)){
      if(is.na(df[i,k]) || (is.character(df[i,k]) && df[i,k] == "")){ #check if they are empty cells or empty chars
        rows_to_remove <- c(rows_to_remove, i)  # Add the index of the row to remove
        break  # Exit the inner loop once an empty cell is found in the row
      }
    }
  }
  df <- df[-rows_to_remove, ]  # Remove the rows with empty cells
  return(df)  # Return the modified dataframe
}

#gets rid of empty rows (dead mice)
weight_measure_file <- emptyRow(weight_measure_file)

#removing useless cols
print(colnames(weight_measure_file))
weight_measure_file <- weight_measure_file[,-c(1:5,7)]

#changing colnames with dates
changeDateCols <- function(df,startDateCol){ #startDateCol corresponds to the position where cols with dates as names start
  for(i in startDateCol:ncol(df)){
   colnames(df)[i] <- gsub("\\.", "-",colnames(df)[i]) #replacing . with -
   colnames(df)[i]<- substring(colnames(df)[i], first = 2) #getting rid of X at the start
   #colnames(weight_measure_file)[i] <- substr(colnames(weight_measure_file)[i], start = 1, stop = nchar(colnames(weight_measure_file)[i]) - 1)
    
  }
  return(df)
}


print(colnames(weight_measure_file))
weight_measure_file <- changeDateCols(weight_measure_file,4)
colnames(weight_measure_file)[1:3] <- c("diet","cage","treatment")

#function that does pivot_longer from cols corresponding to weight measurement on a particular date
#long tidy format
pivotLongerWeightDate <- function(df){
  startDate <- grep("202",colnames(df))[1] #finding at which position the weight measurements date cols start appearing
  return(as.data.frame(pivot_longer(df, c(startDate:ncol(df)), cols_vary = "slowest", names_to = "date", values_to = "weight")))
}

weight_measure_file <- pivotLongerWeightDate(weight_measure_file)

#putting "date" col as a date variable type (not character)
weight_measure_file$date <- as.Date(weight_measure_file$date)

#setting the diet column data into a string format so that it can be used into ggplot
weight_measure_file$diet <- as.character(weight_measure_file$diet)

#setting the weight column data into numeric values
charToNum <- function(df){
  for(i in 1:ncol(df)){
      for(k in 1:nrow(weight_measure_file)){
        isMeasureCol <- grepl(",",df[k,i])
        if(isMeasureCol){
          df[k,i] <- as.numeric(gsub("\\,", ".", df[k,i]))}#replaces every "," by "." and makes it a num variable
      }

  
  }
  return(df)
}

weight_measure_file <- charToNum(weight_measure_file)

#using this enables to verify if date variable is of type "date" and weight of type numeric
str(weight_measure_file)

#replacing abnormal values
weight_measure_file$weight[4] = 24.2
weight_measure_file$weight[28] = 23.1

#Transforming dates into numeric format for statistical measurements
dateToNum<- function(df){
  for(i in 1:ncol(df)){
    if(class(df[,i])=="Date"){ #check if dataframe col is of type "Date"
      referenceDate <- min(df[,i]) #choose a the earliest date as ref
      df$time_numeric <- as.numeric(df[,i] - referenceDate) #substraction by reference and as_numeric for every date
      } 
  } 
  return(df)
}

weight_measure_file <- dateToNum(weight_measure_file)

#creating 4 groups for easier graphic interpretations
for(i in c(1:length(weight_measure_file$cage))){
  if(any(weight_measure_file$cage[i] %in% c(1,3,5))){
    weight_measure_file$group[i] <- "50 ppm FeSO4 + DSS"
  }
  if(any(weight_measure_file$cage[i] %in% c(2,4,6))){
    weight_measure_file$group[i] <- "50 ppm FeSO4 + water"
  }
  if(any(weight_measure_file$cage[i] %in% c(7,9,11))){
    weight_measure_file$group[i] <- "500 ppm FeSO4 + DSS"
  }
  if(any(weight_measure_file$cage[i] %in% c(8,10,12))){
    weight_measure_file$group[i] <- "500 ppm FeSO4 + water"
  }
}


darker_blue = "#007577"
darker_red = "#B5544F"

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














###DSS FOLLOW UP SHEET DATA
#loading the dss followup sheet from github
github_repo <- "gut-microbiota-iron"
file <- "adult_mice_DSS1_followup.csv"
dss_followup <- read.table(text = fetchGHdata(github_repo,file), header = TRUE, sep = "\t")

#if the function is not working
setwd("D:/CHUM_git/gut-microbiota-iron/adult-DSS-exp/")
dss_followup <- read.csv("adult_mice_DSS1_followup.csv", header = TRUE, sep = "\t")

#replacing "bleeding" for "hemoccult" for one of the columns
dss_followup$Day.0.1[1] <- "hemoccult"

#replacing weight values by weight variation values (see disease index calculation)
#((Day X)/(Day 1))×100

dss_followup[,5] <- gsub("\\,", ".", dss_followup[,5]) #the fifth column is the day0 weight (used for percentage calculation)
for(i in seq(from = 8, to = ncol(dss_followup), by = 3)){
  for(n in 2:(nrow(dss_followup))){
    dss_followup[n,i] <- gsub("\\,", ".", dss_followup[n,i]) #replacing "," by "."
    dss_followup[n,i] <- ((as.numeric(dss_followup[n,i])/as.numeric(dss_followup[n,5]))*100)
  }
}

###creating 3 different dataframes for each metric (weight, hemoccult, stool consistency)

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

#removing the first row and the mouse that died before the DSS + mouse that died day before dissection
colnames(dss_weight)[1:4] <- dss_weight[1,1:4]
dss_weight <- dss_weight[-c(1,4,28),]

colnames(dss_consistency)[1:4] <- dss_consistency[1,1:4]
dss_consistency <- dss_consistency[-c(1,4,28),]

colnames(dss_hemo)[1:4] <- dss_hemo[1,1:4]
dss_hemo <- dss_hemo[-c(1,4,28),]

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
#creating a disease severity index, instructions are found in word document
combined_df$weight <- as.numeric(combined_df$weight)
combined_df$weight[1:34] <- 100 #corresponds to values for day0 (all 100%)
combined_df$weight <- combined_df$weight-100 #weight difference
combined_df$weight <- abs(combined_df$weight) #absolute difference
combined_df$index <- 0

#	Weight Variation: 0 - None, 1 - 1%-5%, 2 - 5%-10%, 3 - 10%-20%, 4 - >20%
#Stool Consistency: 0 - Normal, 2 - Loose, 4 – Diarrhea
#Fecal Bleeding: Fecal Bleeding: 0 - None, 2 - Hemoccult Positive, 4 - Gross rectal bleeding


combined_dfa <- combined_df
for(i in 1:nrow(combined_dfa)){
  if(combined_dfa$consistency[i]=="N"){
    combined_dfa$index[i] <- combined_dfa$index[i]+0
  }else{
    combined_dfa$index[i] <- combined_dfa$index[i]+2
  }
  if(combined_dfa$hemo[i]=="No"){
    combined_dfa$index[i] <- combined_dfa$index[i]+0
  }else{
    combined_dfa$index[i] <- combined_dfa$index[i]+4
  }
  if (combined_dfa$weight[i] >= 0 & combined_dfa$weight[i] < 1) {
    combined_dfa$index[i] <- combined_dfa$index[i] + 0
  } else if (combined_dfa$weight[i] >= 1 & combined_dfa$weight[i] < 5) {
    combined_dfa$index[i] <- combined_dfa$index[i] + 1
  } else if (combined_dfa$weight[i] >= 5 & combined_dfa$weight[i] < 10) {
    combined_dfa$index[i] <- combined_dfa$index[i] + 2
  } else if (combined_dfa$weight[i] >= 10 & combined_dfa$weight[i] < 20) {
    combined_dfa$index[i] <- combined_dfa$index[i] + 3
  } else if (combined_dfa$weight[i] >= 20) {
    combined_dfa$index[i] <- combined_dfa$index[i] + 4
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

#replacing all the col names because this is ridiculous
colnames(combined_dfa)[1:5] <- c("cage","diet","treatment","id","date")

#creating scatter plot with the four different treatments (diet combined with dss or control)
#this graph has a disease index score in the y column
data <- as.data.frame(combined_dfa)
data %>%
  ggplot(aes(x = date, y = index, color = diet))+
  stat_summary(aes(group = group, shape = treatment),fun ="mean", geom = "point", size = 3)+
  stat_summary(aes(group = group, linetype = ifelse(grepl("DSS", group), "DSS", "Water")),fun = "mean",geom = "line",size = 1)+
  labs(title = "Disease severity index (DSI) during DSS",
       x = "Day",
       y = "DSI",
       color = "Diet")+
  scale_linetype_manual(name = "Treatment", 
                        values = c("DSS" = "dashed", "Water" = "solid"))+
  scale_x_discrete(labels = c("0", "1", "2", "3", "4", "5"))+
  guides(shape = 'none')+
  theme_minimal()+
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
