#Weight measurements for experiment 1 (DSS old mice)

setwd("I:/Chercheurs/Santos_Manuela/Thibault M/DSS/adult")
#loading libraries
library("tidyverse")
library("ggplot2")


#loading the file
weight_measure_file = read.csv2("adult_dss_weight_measurement.csv")

#removing the 3 dead mice
weight_measure_file = weight_measure_file[-c(16,27,36,37),]

#creating scatter plot with the four different treatments (diet combined with dss or control)
#redefining the data frame by creating a new column for the date of the measurement

weight_measure_file$msr_date = NA
final_file = weight_measure_file
 for (i in weight_measure_file[,c(11:length(weight_measure_file)-1)]){
   weight_measure_file$msr_date = NA
   weight_measure_file$msr_date = names(weight_measure_file[i])
   final_file = rbind(final_file,weight_measure_file)
 }
ggplot(weight_measure_file,aes(x = weight_measure_file[,c(10:23)], y = ))

str(weight_measure_file)


iojwefeewoifejoifeijoejoiwjiowejoiweijowejioweijowejiowejioweijoweijowjiofeijoejoijiofejio