#loading libraries
library("tidyverse")
library("ggplot2")
library("dplyr")

weight_measure_file <- gh_data

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
  labs(title = "Body weight through time",
       x = "Day",
       y = "Weight (g)")+
  theme_minimal()+
  ylim(1,30)






