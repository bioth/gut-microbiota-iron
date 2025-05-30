# loading libraries
{
  library("tidyverse") # loading bunch of packages
  library("ggplot2") # come on, everyone knows what it is used for
  library("dplyr") # arranging and manipulating data easily
  library("geepack") # library for loading GEE tests
  library("lme4") # library for loading ANOVA
  library("car") # for anova too
  library("ggsignif") # adding significance bars to ggplots
  library("readxl") # Read and write excel files
}


# Setting working directory
setwd("D:/CHUM_git/gut-microbiota-iron/")
setwd("I:/Chercheurs/Santos_Manuela/Thibault M/gut-microbiota-iron/")
setwd("C:/Users/Thibault/Documents/CHUM_git/gut-microbiota-iron/")

# Loading functions for data manipulation
source("other scripts/dataManipFunctions.R")

# Young mice + DSS experiment
{
setwd("experiments/finished exp/young-DSS-exp3")
young_weight <- read.csv("young48_weight.csv", header = TRUE, sep = ";")
young_weight <- young_weight[,-c(6:10)] # Remove repeated weight measures at start of experiment
young_weight <- weightDataManipulation(young_weight, groupInfoCols = 4, fromDay0 = TRUE)

# Scatter plot with weight evolution over time for diet exposure
young_weight_plot <- weightPlot(young_weight, percentage = TRUE, diet_only = TRUE)
young_weight_plot

young_weight <- read.csv("young48_weight_cageChanges.csv", header = TRUE, sep = ";")
young_weight <- young_weight[,-c(5:16)] # Remove weight measures before start of DSS
young_weight <- weightDataManipulation(young_weight, groupInfoCols = 4, fromDay0 = TRUE)

# Scatter plot with the four different treatments (diet combined with dss or control)
young_weight_plot <- weightPlot(young_weight, percentage = TRUE, diet_only = FALSE)
young_weight_plot

# Mice DAI evolution
young_dss_followup <- read.csv("young48_dss_followup.csv", header = TRUE, sep = ";")
young_dss_followup <- dssFollowupManipulation(df = young_dss_followup,groupInfoCols = 4,dateStart = "2024-05-29",nbrDays = 5, negativeOnly = TRUE)
young_dssflwup_plot <- dssDiseaseIndexPlot(young_dss_followup)
young_dssflwup_plot

# Statistics for last day of DAI
verifyStatsAssumptions(df = young_dss_followup[young_dss_followup$time_numeric == "5",], group = "gg_group", measure = "index")
wilcox.test(index ~ gg_group, data = young_dss_followup[young_dss_followup$time_numeric == "5",])

# Mice dissection data
dissec_young <- read.csv("young48_dss_dissection.csv", sep = ";", header = TRUE)
dissec_young <- dissectionDataManipulation(dissec_young, groupInfoCols = 4, numerical = FALSE)
dissec_young$gg_group <- factor(dissec_young$gg_group, levels = c("50 water","500 water","50 dss","500 dss"))
dissec_young$colon_length_nrm <- dissec_young$colon_length/dissec_young$body_weight

#boxplot for body weight
young_dissec_bw <- dissecBoxplot(dissec_young,"body") 
young_dissec_bw

# Statistics
verifyStatsAssumptions(df = dissec_young, group = "gg_group", measure = "body_weight")
pairwise.t.test( # Pairwise welch t tests
  x = dissec_young$body_weight,
  g = dissec_young$gg_group,         # factor defining the groups
  p.adjust.method = "BH",      # adjust for multiple testing (e.g., "bonferroni", "BH", etc.)
  pool.sd = FALSE              # FALSE = Welch's t-test
)

#boxplot for std liver weight
young_dissec_lvr <- dissecBoxplot(dissec_young,"liver") 
young_dissec_lvr

# Statistics
verifyStatsAssumptions(df = dissec_young, group = "gg_group", measure = "std_liver_weigth") 
anova <- aov(std_liver_weigth ~ diet * treatment, data = dissec_young) #Fit a ANOVA model
summary(anova)
results <- TukeyHSD(anova) # Perform Tukey's HSD test and store the results in the list
results

# boxplot for std spleen weight
young_dissec_spln <- dissecBoxplot(dissec_young,"spleen") 
young_dissec_spln

# Statistics
verifyStatsAssumptions(df = dissec_young, group = "gg_group", measure = "std_spleen_weigth")
pairwise.t.test( # Pairwise welch t tests
  x = dissec_young$std_spleen_weigth,
  g = dissec_young$gg_group,         # factor defining the groups
  p.adjust.method = "BH",      # adjust for multiple testing (e.g., "bonferroni", "BH", etc.)
  pool.sd = FALSE              # FALSE = Welch's t-test
)

#boxplot for colon length (non std)
young_dissec_cln <- dissecBoxplot(dissec_young,"colon") 
young_dissec_cln

# Statistics
verifyStatsAssumptions(df = dissec_young, group = "gg_group", measure = "colon_length")
pairwise.t.test( # Pairwise welch t tests
  x = dissec_young$colon_length,
  g = dissec_young$gg_group,         # factor defining the groups
  p.adjust.method = "BH",      # adjust for multiple testing (e.g., "bonferroni", "BH", etc.)
  pool.sd = FALSE              # FALSE = Welch's t-test
)

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

p = ironBoxplot(df, "iron_concentration", group = "diet", title = "Iron concentration in stools at day 35", y_axis_title = "yg of iron per g of stools", custom_colors = c("blue","red"))
p <- p+
  scale_x_discrete(labels = c("50 ppm","500 ppm"))+
  labs(y = "µg iron per g of stool", title = "Iron in stools at end of diet exposure")+
  guides(color = "none")
p

# Stats 
verifyStatsAssumptions(df, "diet" , "iron_concentration")
wilcox.test(iron_concentration ~ diet, data = df)

#Tf ferrozine assay for stools
sheets <- read_excel_allsheets("young48_dss_ferrozine_tF.xlsx")
df <- as.data.frame(sheets["Ferrozine Stool Tfinal"])
df <- df[4:53,c(3,14:16)]
colnames(df) <- df[1,]
colnames(df)[1:2] <- c("id","iron_concentration")
df <- df[-c(1,26,48),]
df$gg_group <- paste(df$treatment, "+", df$diet, sep = "")
df$gg_group <- factor(df$gg_group, levels = c("water+50","dss+50","water+500","dss+500"))
df$iron_concentration <- as.numeric(df$iron_concentration)

p = ironBoxplot(df, "iron_concentration", group = "gg_group", title = "Iron concentration in stools at final day", y_axis_title = "yg of iron per g of stools", custom_colors = c("blue","red","darkblue", "darkred"))
p+
  scale_x_discrete(labels = c("50 ppm\nCtrl","50 ppm\nDSS","500 ppm\nCtrl","500 ppm\nDSS"))+
  labs(y = "µg iron per g of liver", title = "")+
  guides(color = "none")
p

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




}















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




setwd("experiments/ongoing exp/young-abx-exp6/")

# Load metadata to retrieve group information
meta <- read_xlsx("dissection.xlsx")
meta <- meta[,c(1:4)]
meta$ID <- substring(meta$ID, first = 1, last = 5)
meta <- meta[-17,] # Remove dead mouse

#Tfinal ferrozine assay for liver
sheets <- read_excel_allsheets("Ferrozine liver and spleen.xlsx") # Find custom function read_excel_allsheets above
df <- as.data.frame(sheets["Ferrozine LIVER"])
df <- df[,c(3,14)]
colnames(df)[1:2] <- c("ID","iron_concentration")
df <- df[-c(1:3),] 
df$iron_concentration <- as.numeric(df$iron_concentration)
df <- merge(df, meta, by ="ID") # Bind metadata information to df
df$diet <- factor(df$diet, levels = c("50","500"))
df$treatment <- factor(df$treatment, levels = c("water", "abx"))
df$gg_group <- factor(paste0(df$diet, ":", df$treatment), levels = c("50:water", "500:water", "50:abx", "500:abx"))

p = ironBoxplot(df, "iron_concentration", display_significance_bars = F, title = "Total iron in liver", y_axis_title = "µg Fe/g of liver", custom_colors = c("blue","red","darkblue","darkred"), path = NULL)
p <- p+
  scale_x_discrete(labels = c("50 ppm","500 ppm"))+
  labs(y = "µg Fe/g of stool", title = "")+
  ylim(0,NA)+
  guides(color = "none")+
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 10))
p

# Levene's Test for homogeneity of variance
leveneTest(iron_concentration ~ gg_group, data = df)

# Shapiro test per group for normality assumptions
by(df$iron_concentration, df$gg_group, shapiro.test)

# Pairwise welch t tests
pairwise.t.test(
  x = df$iron_concentration,
  g = df$gg_group,         # factor defining the groups
  p.adjust.method = "BH",      # adjust for multiple testing (e.g., "bonferroni", "BH", etc.)
  pool.sd = FALSE              # FALSE = Welch's t-test
)

df <- as.data.frame(sheets["Ferrozine SPLEEN"])
df <- df[,c(3,14)]
colnames(df)[1:2] <- c("ID","iron_concentration")
df <- df[-c(1:3),] 
df$iron_concentration <- as.numeric(df$iron_concentration)
df <- merge(df, meta, by ="ID") # Bind metadata information to df
df$diet <- factor(df$diet, levels = c("50","500"))
df$treatment <- factor(df$treatment, levels = c("water", "abx"))
df$gg_group <- factor(paste0(df$diet, ":", df$treatment), levels = c("50:water", "500:water", "50:abx", "500:abx"))

p = ironBoxplot(df, "iron_concentration", display_significance_bars = F, title = "Total iron in spleen", y_axis_title = "µg Fe/g of liver", custom_colors = c("blue","red","darkblue","darkred"), path = NULL)
p

# Levene's Test for homogeneity of variance
leveneTest(iron_concentration ~ gg_group, data = df)

# Shapiro test per group for normality assumptions
by(df$iron_concentration, df$gg_group, shapiro.test)

#Fit a ANOVA model for the current time point
anova <- aov(iron_concentration ~ diet * treatment, data = df)
summary(anova)
# Perform Tukey's HSD test and store the results in the list
results <- TukeyHSD(anova)
print(results)
}

