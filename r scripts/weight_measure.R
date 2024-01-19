#Weight measurements for experiment 1 (DSS old mice)

setwd("I:/Chercheurs/Santos_Manuela/Thibault M/DSS/adult")
#loading libraries
library("tidyverse")
library("ggplot2")
library("dplyr")



#loading the file
weight_measure_file = read.csv2("adult_dss_weight_measurement.csv")

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

#creating scatter plot with the four different treatments (diet combined with dss or control)
ggplot(weight_measure_file,aes(x = date, y = weight))


#using this enables to verify is a variable is of type "time"
str(weight_measure_file)





