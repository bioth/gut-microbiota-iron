#loading libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
BiocManager::install("DESeq2")

library(phyloseq)
library(DESeq2)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(data.table)

theme_set(theme_bw())

#for microbiota 17
#set working directory
setwd("~/Documents/CHUM_git/Microbiota_17/")
asv_table <- as.data.frame(fread("asv_table/seqtab.nochim_run_m1.csv", sep = ";"))

# Set the first column as row names and remove it from the data frame
rownames(asv_table) <- asv_table[,1]  # Use the first column as row names
asv_table <- asv_table[,-1]  # Drop the first column

#loading metadata of interest
metadata <- read.csv("metadata/metadata.csv", sep = ";", row.names = 1)

{
#associating each sample_reads_id with its metadata (diet, abx etc)
nas <- NULL
for(i in 1:nrow(df)){
  
  if(!all(!grepl(df$SampleID[i],rownames(asv_table)))){
    rownames(df)[i] <- rownames(asv_table)[grep(df$SampleID[i],rownames(asv_table))]
  }
  else{
    #remove metadata row if sampleID not found in sample_reads_id
    nas <- append(nas, i)
  }
}

#remove metadata row if sampleID not found in sample_reads_id
df <- df[-c(nas),]
#replace "f" time point with appropriate number of days 
df$timepoint <- ifelse(df$timepoint == "f", 80, df$timepoint)
#add treatment and diet information
for(i in 1:nrow(df)){
  df$treatment[i] <- ifelse(grepl("no donor",df$group[i]), "ctrl", "abx")
  df$diet[i] <- ifelse(grepl("500",df$group[i]), "500", "50")
}
df$timepoint <- as.numeric(df$timepoint)
} #paying with metadata, could be useful later

#load taxonomical assignments
taxa <- as.matrix(fread("taxonomy/taxa_annotation_m1.csv", sep = ";"))

# Set the first column as row names and remove it from the data frame
rownames(taxa) <- taxa[,1]  # Use the first column as row names
taxa <- taxa[,-1]  # Drop the first column

#creating phyloseq object
ps <- phyloseq(otu_table(asv_table, taxa_are_rows = FALSE),
               sample_data(metadata),
               tax_table(taxa))

#use short names for the asvs (eg ASV21) rather than full dna sequence name
#though we can keep these dna sequences for other purposes (?)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

#create two separate ps objects for each Claire and Samuel's datasets
ps_samuel <- subset_samples(ps, student == "Samuel")
ps_claire <- subset_samples(ps, student == "Claire")

#converting ps object to deseq object
deseq_samuel <- phyloseq_to_deseq2(ps_samuel, ~genotype+treatment) 
deseq_claire <- phyloseq_to_deseq2(ps_claire, ~diet) 




















#visualize mean alpha diversity
# Extract data and compute means
richness_data <- estimate_richness(ps, measures = c("Shannon", "Simpson"))

# Extract sample data
sample_data <- as.data.frame(sample_data(ps))

# Combine richness data with sample data
combined_data <- cbind(sample_data, richness_data)

# Compute mean values for each combination of timepoint, diet, and treatment
mean_data <- combined_data %>%
  group_by(timepoint, diet, treatment) %>%
  summarize(across(c(Shannon, Simpson), mean, na.rm = TRUE))

# Ensure appropriate column names
mean_data <- mean_data %>%
  rename(Timepoint = timepoint, Diet = diet, Treatment = treatment)

# Shannon Diversity Plot
ggplot(mean_data, aes(x = Timepoint, y = Shannon, color = Diet, shape = Treatment)) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(Diet, Treatment)), linetype = "dashed") +
  labs(title = "Mean Shannon Diversity by Timepoint", y = "Mean Shannon Index") +
  theme_minimal()

# Simpson Diversity Plot
ggplot(mean_data, aes(x = Timepoint, y = Simpson, color = Diet, shape = Treatment)) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(Diet, Treatment)), linetype = "dashed") +
  labs(title = "Mean Simpson Diversity by Timepoint", y = "Mean Simpson Index") +
  theme_minimal()

#plot_richness(ps, x = "timepoint", measures=c("Shannon","Simpson"), color = "diet", shape = "treatment")









#for microbiota 17
#set working directory
setwd("D:/CHUM_git/Microbiota_17/asv_table/")
asv_table <- read.csv("seqtab.nochim_run2.csv", sep = ";")

#loading metadata of interest
setwd("../metadata/")
metadata <- read.csv("metadata.csv", sep = ";")
df <- metadata

#associating each sample_reads_id with its metadata (diet, abx etc)
nas <- NULL
for(i in 1:nrow(df)){
  
  if(!all(!grepl(df$sample_id[i],rownames(asv_table)))){
    rownames(df)[i] <- rownames(asv_table)[grep(df$sample_id[i],rownames(asv_table))]
  }
  else{
    #remove metadata row if sampleID not found in sample_reads_id
    nas <- append(nas, i)
  }
}

#remove metadata row if sampleID not found in sample_reads_id
df <- df[-c(nas),]

#add col counting number of days
df$days <- df$week*7

#load taxonomical assignments
taxa <- read.csv("../taxonomy/taxa_annotation.csv", sep = ";")
taxa <- as.matrix(taxa)

#make df and asv_table rownames the same
rownames(asv_table) <- stringr::str_remove_all(rownames(asv_table), "_R1_filt.fastq.gz")

#arrange df rows in alphabetical order
df <- df[order(df$sample_id),]

#put sample_id as rownames
rownames(df) <- df$sample_id

#check if rownames are the same
identical(rownames(df), rownames(asv_table))

#creating phyloseq object
ps <- phyloseq(otu_table(asv_table, taxa_are_rows = FALSE),
               sample_data(df),
               tax_table(taxa))

#use short names for the asvs (eg ASV21) rather than full dna sequence name
#though we can keep these dna sequences for other purposes (?)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps



#visualize mean alpha diversity
# Extract data and compute means
richness_data <- estimate_richness(ps, measures = c("Shannon", "Simpson"))

# Extract sample data
sample_data <- as.data.frame(sample_data(ps))

# Combine richness data with sample data
combined_data <- cbind(sample_data, richness_data)


mean_cl_normal <- function(x, mult = 1.96) { #mult is 1.96 for a 95% confidence interval
  # Calculate the mean of the input vector x
  mean_val <- mean(x, na.rm = TRUE)
  
  # Calculate the standard error of the mean
  se_val <- sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
  
  # Return a data frame with the mean (y), and the lower (ymin) and upper (ymax) bounds
  data.frame(y = mean_val, ymin = mean_val - mult * se_val, ymax = mean_val + mult * se_val)
}


# Compute mean values for each combination of timepoint, diet, and treatment
samuel_data <- combined_data %>%
  filter(student == "Samuel") %>%
  mutate(gg_group = paste(genotype, treatment, sep = "_")) 

#Shannon index
anova_result <- aov(Shannon ~ genotype * treatment, data = samuel_data)

tukey_results <- TukeyHSD(anova_result)
factor_names <- names(tukey_results)

# Prepare to combine results
tukey_df_list <- lapply(factor_names, function(factor) {
  tukey_result <- tukey_results[[factor]]
  if (nrow(tukey_result) > 0) {  # Check if there are any results
    tukey_df <- as.data.frame(tukey_result)
    tukey_df$comparison <- rownames(tukey_df)
    tukey_df$factor <- factor
    return(tukey_df)
  } else {
    return(NULL)  # Return NULL if no results
  }
})

# Remove NULL elements and combine all results into one data frame
tukey_df_all <- do.call(rbind, Filter(Negate(is.null), tukey_df_list))

# Add significance stars based on p-values
tukey_df_all$significance <- cut(tukey_df_all$`p adj`,
                                 breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
                                 labels = c("***", "**", "*", ""),
                                 include.lowest = TRUE)

  
# Shannon Diversity Plot
ggplot(samuel_data, aes(x = gg_group, y = Shannon, color = treatment, shape = genotype)) +
  geom_point(size = 3) +
  labs(title = "Mean Shannon Diversity", y = "Shannon Index") +
  scale_shape_discrete(limits = c("Wt","IL-22ra1-/-"))+
  scale_color_discrete(limits = c("DSS + EcNC101 + Vehicle", "DSS + EcNC101 + Putrescine"))+
  scale_x_discrete(limits = c("Wt_DSS + EcNC101 + Vehicle","Wt_DSS + EcNC101 + Putrescine","IL-22ra1-/-_DSS + EcNC101 + Vehicle","IL-22ra1-/-_DSS + EcNC101 + Putrescine"))+
  stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..-0.25, yend=..y..), color = "black", linewidth =1)+ #adding horizontal bars representing means
  stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..+0.25, yend=..y..), color = "black", linewidth =1)+ 
  stat_summary(fun.data="mean_cl_normal", geom="errorbar", color="red", width=0.2, size = 0.7) + #adding SEM error bars
  ylim(min(samuel_data$Shannon)-0.25,max(samuel_data$Shannon)+0.25)+
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
    axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
    legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars

setwd("../../some photos/")
ggsave("shannon_diversity.png", width = 8, height = 8, dpi = 300, bg = "white")

# Simpson Diversity Plot
anova_result <- aov(Simpson ~ genotype * treatment, data = samuel_data)
summary(anova_result)
tukey_results <- TukeyHSD(anova_result)
print(tukey_results)
factor_names <- names(tukey_results)

# Prepare to combine results
tukey_df_list <- lapply(factor_names, function(factor) {
  tukey_result <- tukey_results[[factor]]
  if (nrow(tukey_result) > 0) {  # Check if there are any results
    tukey_df <- as.data.frame(tukey_result)
    tukey_df$comparison <- rownames(tukey_df)
    tukey_df$factor <- factor
    return(tukey_df)
  } else {
    return(NULL)  # Return NULL if no results
  }
})

# Remove NULL elements and combine all results into one data frame
tukey_df_all <- do.call(rbind, Filter(Negate(is.null), tukey_df_list))

# Add significance stars based on p-values
tukey_df_all$significance <- cut(tukey_df_all$`p adj`,
                                 breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
                                 labels = c("***", "**", "*", ""),
                                 include.lowest = TRUE)


# Shannon Diversity Plot
ggplot(samuel_data, aes(x = gg_group, y = Simpson, color = treatment, shape = genotype)) +
  geom_point(size = 3) +
  labs(title = "Mean Simpson Diversity", y = "Simpson Index") +
  scale_shape_discrete(limits = c("Wt","IL-22ra1-/-"))+
  scale_color_discrete(limits = c("DSS + EcNC101 + Vehicle", "DSS + EcNC101 + Putrescine"))+
  scale_x_discrete(limits = c("Wt_DSS + EcNC101 + Vehicle","Wt_DSS + EcNC101 + Putrescine","IL-22ra1-/-_DSS + EcNC101 + Vehicle","IL-22ra1-/-_DSS + EcNC101 + Putrescine"))+
  stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..-0.25, yend=..y..), color = "black", linewidth =1)+ #adding horizontal bars representing means
  stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..+0.25, yend=..y..), color = "black", linewidth =1)+ 
  stat_summary(fun.data="mean_cl_normal", geom="errorbar", color="red", width=0.2, size = 0.7) + #adding SEM error bars
  ylim(min(samuel_data$Simpson)-0.001,max(samuel_data$Simpson)+0.001)+
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
    axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
    legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars

setwd("../../some photos/")
ggsave("simpson_diversity.png", width = 8, height = 8, dpi = 300, bg = "white")



#going for a bacterial relative abundance barplot
# Extract the taxonomy table from the phyloseq object
tax_table <- tax_table(ps)

# Create a data frame from the phyloseq object
physeq_df <- psmelt(ps)




# Group by Phylum and calculate the relative abundance
physeq_df <- physeq_df %>%
  filter(student == "Samuel") %>%
  mutate(gg_group = paste(genotype, treatment, sep = "_"))%>%
  group_by(gg_group, Phylum) %>%
  summarize(Abundance = sum(Abundance)) %>%
  group_by(gg_group) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance) * 100)

# colnames(physeq_df)[1] = "sample_id"
# physeq_df <- physeq_df %>%
#   left_join(samuel_data, by = "sample_id")
# physeq_df <- physeq_df[,-c(5,7,9,11,12,13,14)]
# # Remove rows where gg_group is NA
# physeq_df <- physeq_df %>%
#   filter(!is.na(gg_group))

library(ggforce)

# Plot the relative abundance by phylum
ggplot(physeq_df, aes(x = gg_group, y = Relative_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "right")+
  labs(x = "Genotype + Treatment", y = "Relative Abundance (%)", title = "Phylum Distribution") +
  scale_fill_brewer(palette = "Paired")+ # Optional: Use a color palette for better visualization
  facet_wrap(~ gg_group, scales = "free_x")


ggplot(physeq_df, aes(x = "", y = Relative_Abundance, fill = Phylum)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  facet_wrap(~ gg_group) +
  theme_void() +
  labs(title = "Bacterial relative abundance by phyla")+
  theme(
    plot.title = element_text(size = 16, face = "bold", vjust = 5),  # Adjust title font size and style
    legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend font size
  )

setwd("../../some photos/")
ggsave("bacterial_relative_abundance.png", width = 8, height = 8, dpi = 300, bg = "white")

anova_result <- aov(Relative_Abundance ~ gg_group * Phylum, data = physeq_df)
summary(anova_result)
tukey_results <- TukeyHSD(anova_result)
print(tukey_results)
factor_names <- names(tukey_results)

# Prepare to combine results
tukey_df_list <- lapply(factor_names, function(factor) {
  tukey_result <- tukey_results[[factor]]
  if (nrow(tukey_result) > 0) {  # Check if there are any results
    tukey_df <- as.data.frame(tukey_result)
    tukey_df$comparison <- rownames(tukey_df)
    tukey_df$factor <- factor
    return(tukey_df)
  } else {
    return(NULL)  # Return NULL if no results
  }
})

# Remove NULL elements and combine all results into one data frame
tukey_df_all <- do.call(rbind, Filter(Negate(is.null), tukey_df_list))

# Add significance stars based on p-values
tukey_df_all$significance <- cut(tukey_df_all$`p adj`,
                                 breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
                                 labels = c("***", "**", "*", ""),
                                 include.lowest = TRUE)







# Compute mean values for each combination of timepoint, diet, and treatment
mean_data_claire <- combined_data %>%
  filter(student == "Claire") %>%
  group_by(diet,days) %>%
  summarize(across(c(Shannon, Simpson), mean, na.rm = TRUE))

claire_t1 <- combined_data %>%
  filter(student =="Claire") %>%
  filter(week == 3)
print(t.test(Shannon ~ diet, data = claire_t1))

claire_t2 <- combined_data %>%
  filter(student =="Claire") %>%
  filter(week == 8)
print(t.test(Shannon ~ diet, data = claire_t2))

claire_t3 <- combined_data %>%
  filter(student =="Claire") %>%
  filter(week == 10)
print(t.test(Shannon ~ diet, data = claire_t3))

claire_t4 <- combined_data %>%
  filter(student =="Claire") %>%
  filter(week == 14)
print(t.test(Shannon ~ diet, data = claire_t4))





# Shannon Diversity Plot
ggplot(mean_data_claire, aes(x = days, y = Shannon, color = diet, shape = diet)) +
  geom_point(size = 3) +
  labs(title = "Mean Shannon Diversity", y = "Mean Shannon Index") +
  # ylim(0,max(mean_data_claire$Shannon))+
  theme_minimal()

anova_shannon <- aov(Shannon ~ diet, data = combined_data)
summary(anova_result)

# Simpson Diversity Plot
ggplot(mean_data_claire, aes(x = days, y = Simpson, color = diet, shape = diet)) +
  geom_point(size = 3) +
  labs(title = "Mean Simpson Diversity", y = "Mean Simpson Index") +
  theme_minimal()

claire_t1 <- combined_data %>%
  filter(student =="Claire") %>%
  filter(week == 3)
print(t.test(Simpson ~ diet, data = claire_t1))

claire_t2 <- combined_data %>%
  filter(student =="Claire") %>%
  filter(week == 8)
print(t.test(Simpson ~ diet, data = claire_t2))

claire_t3 <- combined_data %>%
  filter(student =="Claire") %>%
  filter(week == 10)
print(t.test(Simpson ~ diet, data = claire_t3))

claire_t4 <- combined_data %>%
  filter(student =="Claire") %>%
  filter(week == 14)
print(t.test(Simpson ~ diet, data = claire_t4))









#bar plot
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

ps.top20_claire <- ps.top20 %>%
  filter(student == "Claire")

plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")

#plot_richness(ps, x = "timepoint", measures=c("Shannon","Simpson"), color = "diet", shape = "treatment")
