#loading libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
BiocManager::install("DESeq2")
BiocManager::install("DECIPHER")

devtools::install_github("ThibaultCuisiniere/StackbarExtended")



#Loading required packages
{
library(phyloseq)
library(DESeq2) #for differential abundance analysis
library(Biostrings)
library(ggplot2)
library(dplyr)
library(data.table)
library(ggsignif) #for adding significance bars to a ggplot
library(car)
library(ggpubr)
library(patchwork)
library(cowplot) #for adding texts to ggplots
library(vegan) #for beta diversity PCOAs and unifrac
library(DECIPHER) #for performing multiple sequence alignments
library(phangorn) #for building trees
library(ape)
library(tidyverse)
library(writexl)
}

#To unload all packages
{
# Get all loaded packages
loaded_pkgs <- setdiff(names(sessionInfo()$otherPkgs), "base")

# Unload each package
lapply(loaded_pkgs, function(pkg) {
  detach(paste("package", pkg, sep = ":"), character.only = TRUE, unload = TRUE)
})
}

#function to add SEM (1.96 for 95% confidence interval)
mean_cl_normal <- function(x, mult = 1.96) { #mult is 1.96 for a 95% confidence interval
  # Calculate the mean of the input vector x
  mean_val <- mean(x, na.rm = TRUE)
  
  # Calculate the standard error of the mean
  se_val <- sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
  
  # Return a data frame with the mean (y), and the lower (ymin) and upper (ymax) bounds
  data.frame(y = mean_val, ymin = mean_val - mult * se_val, ymax = mean_val + mult * se_val)
}

#Function checking if a dir exists and creating it otherwise
existingDirCheck <- function(path){
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Directory created: ", path)
  } else {
    message("Directory already exists: ", path)
  }
  
}

#Function to remove "/" and "-" characters from a string + lowercase
clean_string <- function(input_string) {
  result <- tolower(gsub("[/-]", "", input_string))
  return(result)
}

#for microbiota 17
#set working directory
setwd("~/Documents/CHUM_git/Microbiota_17/")
asv_table <- as.data.frame(fread("asv_table/seqtab.nochim_run_m1.csv", sep = ";"))

# Set the first column as row names and remove it from the data frame
rownames(asv_table) <- asv_table[,1]  # Use the first column as row names
asv_table <- asv_table[,-1]  # Drop the first column

#loading metadata of interest
metadata <- read.csv("metadata/metadata.csv", sep = ";")

#adding first col as rownames too
rownames(metadata) <- metadata$sample_id

#transforming week col from num to character
metadata$week <- as.factor(metadata$week)
metadata$diet <- as.factor(metadata$diet)

# Replaces Samuel names by shorter versions
metadata$treatment <- gsub(".*Putrescine.*", "Putrescine", metadata$treatment)
metadata$treatment <- gsub(".*Vehicle.*", "Vehicle", metadata$treatment)
metadata$genotype <- gsub(".*IL-22.*", "IL-22KO", metadata$genotype)

# Creates gg_group specific to Claire
metadata$gg_group[metadata$student == "Claire"] <- 
  paste(metadata$week[metadata$student == "Claire"], 
        metadata$diet[metadata$student == "Claire"], 
        sep = ":")

# Creates gg_group specific to Samuel
metadata$gg_group[metadata$student == "Samuel"] <- 
  paste(metadata$genotype[metadata$student == "Samuel"], 
        metadata$treatment[metadata$student == "Samuel"], 
        sep = ":")

{
# Check if the rownames of metadata match the colnames of the ASV table
if(isFALSE(all(rownames(metadata) %in% colnames(asv_table)))){
  metadata <- metadata[rownames(metadata) %in% rownames(asv_table), ]
} # Should return TRUE

# Find row names in asv_table that are not in metadata
missing_rows <- setdiff(rownames(metadata), rownames(asv_table))

# Print the missing row names
print(missing_rows)
  } #checking if they are missing samples

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
} #playing with metadata, could be useful later

#load taxonomical assignments
taxa <- as.matrix(fread("taxonomy/taxa_annotation_m1.csv", sep = ";"))

# Set the first column as row names and remove it from the data frame
rownames(taxa) <- taxa[,1]  # Use the first column as row names
taxa <- taxa[,-1]  # Drop the first column

# Load phylogenetic tree if possible
tree <- read.tree("~/Documents/CHUM_git/figures/samuel/beta_diversity/phylo_tree/phylogenetic_tree.newick")

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
names(dna) <- taxa_names(ps)

#Creating phylogenetic tree
{
#Rarefaction curve, not necessary when lots of high quality data such as here
{
  tab = otu_table(ps)
  class(tab) <- "matrix" # as.matrix() will do nothing
  tab <- t(tab)
  rarecurve(tab, step=10000, lwd=2, ylab="ASV",  label=F)
}

#adding a phylogenetic tree
# Align sequences
alignment <- AlignSeqs(DNAStringSet(dna), processors = 20)

# Build a distance matrix
dist_matrix <- DistanceMatrix(alignment, processors = 20)

#build tree using neighbor-joining method
tree  <- nj(dist_matrix)

# Export the tree as a Newick file
write.tree(tree, file = "~/Documents/CHUM_git/figures/samuel/beta_diversity/phylo_tree/phylogenetic_tree.newick")

#refinement with maximum likelihood
{
  # Convert the alignment to a phangorn-compatible format
  alignment_matrix <- as.matrix(alignment)
  
  # Estimate the substitution model (e.g., GTR)
  fitGTR <- pml(initial_tree, data = alignment_matrix)
  
  # Optimize the ML tree
  ml_tree <- optim.pml(fitGTR, model = "GTR", rearrangement = "stochastic")
}
}

# Add tree to phyloseq object
ps <- merge_phyloseq(ps, phy_tree(tree))

#create two separate ps objects for each Claire and Samuel's datasets
ps_samuel <- subset_samples(ps, student == "Samuel")
ps_claire <- subset_samples(ps, student == "Claire")

#function filtering out ASVs present in less than a threshold of the samples
ps_samuel <- prune_taxa(taxa_sums(ps_samuel) > 10, ps_samuel)
ps_claire <- prune_taxa(taxa_sums(ps_claire) > 10, ps_claire)

#for Samuel's data, put gg_group as factor and define order
sample_data(ps_samuel)$gg_group <- factor(sample_data(ps_samuel)$gg_group, levels = c("Wt:Vehicle", "Wt:Putrescine", "IL-22KO:Vehicle", "IL-22KO:Putrescine"))  # Vehicle as reference


#alpha diversity measures
{
  #Estinate richness measures for dataset
  richness_data <- estimate_richness(ps_samuel, measures = c("Shannon", "Simpson","InvSimpson","Chao1","Observed"))
  
  #Add sample metadata to richness dataframe
  richness_data <- cbind(as.data.frame(sample_data(ps_samuel)), richness_data)



# Add gg_group variable
# richness_data <- richness_data %>%
#   mutate(gg_group = paste(genotype, treatment, sep = ":"))

#for Claire
# richness_data <- richness_data %>%
#   mutate(gg_group = paste(diet, week, sep = "_"))

#testing stuff
# test <- richness_data %>%
#   select() %>%
#   summarize(mean_value = mean(Shannon))

# test <- richness_data %>%
#   filter(gg_group == "500_8") %>%
#   summarise(mean_cl_normal(Shannon))

#Alpha diversity graphs for Samuel
alpha_diversity_graph <- function(measure, data){
  
  #calculate means to set limits
  # mean_data <- data %>%
  #   group_by(gg_group) %>%
  #   summarise(mean_value = mean(measure, na.rm = TRUE))
  
  # #Choosing appropriate statistical tests
  anova_result <- aov(data[[measure]] ~ genotype * treatment, data = data)

  # Extract the residuals
  residuals <- residuals(anova_result)

  # # Shapiro-Wilk test for normality
  print("Checking for normality")
  shapiro_test <- shapiro.test(residuals)
  print(shapiro_test)

  #Verifying homoscedasticity assumption
  print("Checking for homoscedasticity")
  levene_test <- leveneTest(data[[measure]] ~ genotype * treatment, data = data)
  print(levene_test)

  cat("Which test do you want to perform ? (ANOVA: a | KRUSKAL-WALLIS: k (for non-normally distributed data))")

  # Get the user's response
  test_chosen <- readline()

  if(test_chosen == "a"){

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
    stats <- tukey_df_all[-c(1:2),]
    #return(tukey_df_all)

  }
  if(test_chosen == "k"){
    group1 <- data %>%
      filter(gg_group %in% c("Wt_DSS + EcNC101 + Vehicle", "Wt_DSS + EcNC101 + Putrescine"))
    group2 <- data %>%
      filter(gg_group %in% c("IL-22ra1-/-_DSS + EcNC101 + Vehicle", "IL-22ra1-/-_DSS + EcNC101 + Putrescine"))
    group3 <- data %>%
      filter(gg_group %in% c("Wt_DSS + EcNC101 + Vehicle", "IL-22ra1-/-_DSS + EcNC101 + Vehicle"))
    group4 <- data %>%
      filter(gg_group %in% c("Wt_DSS + EcNC101 + Putrescine", "IL-22ra1-/-_DSS + EcNC101 + Putrescine"))

    kruskal_result <- kruskal.test(group1[[measure]] ~ gg_group, data = group1)
    print(kruskal_result)
    posthoc_test <- print(kwAllPairsDunnTest(group1[[measure]] ~ gg_group, data = group1))

    kruskal_result <- kruskal.test(group2[[measure]] ~ gg_group, data = group2)
    print(kruskal_result)
    posthoc_test <- print(kwAllPairsDunnTest(group2[[measure]] ~ gg_group, data = group2))

    kruskal_result <- kruskal.test(group3[[measure]] ~ gg_group, data = group3)
    print(kruskal_result)
    posthoc_test <- print(kwAllPairsDunnTest(group3[[measure]] ~ gg_group, data = group3))

    kruskal_result <- kruskal.test(group4[[measure]] ~ gg_group, data = group4)
    print(kruskal_result)
    posthoc_test <- print(kwAllPairsDunnTest(group4[[measure]] ~ gg_group, data = group4))
  }
  #statistics


    # # Graph claire 
  {
    # ggplot(data, aes(x = (week-3)*7, y = .data[[measure]])) +
    #   stat_summary(aes(group = diet), fun= mean, geom = "point" , size = 2, color = "black")+ #mapping=aes(xend=..x..-0.25, yend=..y..),
    #   geom_point()+
    #   #stat_summary(aes(color = diet),fun = mean, geom = "line", size = 1) + # Mean lines
    #   stat_summary(aes(fill = diet), fun.data="mean_cl_normal", geom = "ribbon", alpha = 0.1)+
    #   labs(title = paste("Mean",measure,"diversity"), y = paste(measure, "index"), x = "Time (days)") +
    #   #ylim(min(mean_data$mean_value),max(mean_data$mean_value))+
    #   theme(
    #     plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    #     axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
    #     axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    #     axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust x-axis tick label font size
    #     axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
    #     legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    #     legend.text = element_text(size = 12),  # Adjust legend font size
    #     panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
    #     panel.grid.minor = element_blank(),  # Remove minor grid lines
    #     axis.line = element_line(color = "black", size = 1))# Include axis lines  # Include axis bars              # Optional: use a minimal theme
  }
  
  p <- ggplot(data, aes(x = genotype, y = data[[measure]], color = treatment)) +
    geom_point(size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) +
    
    # Error bars
    stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
                 aes(color = treatment),
                 width = 0.2, size = 0.7,
                 position = position_dodge(-0.75)) +
    
    #Mean lines
    stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
                 aes(ymin = ..y.., ymax = ..y.., group = treatment),
                 color = "black", linewidth = 0.5, width = 0.5,
                 position = position_dodge(-0.75))+
  
    labs(title = paste("Mean",measure,"diversity"), y = paste(measure, "index")) +
    scale_x_discrete(limits = c("Wt", "IL-22ra1-/-"))+
    scale_color_discrete(limits = c("DSS + EcNC101 + Vehicle", "DSS + EcNC101 + Putrescine"))+
    
    theme(
      plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
      axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust x-axis tick label font size
      axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
      legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
      legend.text = element_text(size = 12),  # Adjust legend font size
      panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
  
    table <- ggtexttable(stats[,-c(1:3)], rows = NULL, 
                       theme = ttheme(base_style = "mCyanWhite", base_size = 5))
    
    # Adding a note below the table
    # note <- ggtexttable(data.frame(Note = "Statistical test: ANOVA"),
    #                     rows = NULL,
    #                     theme = ttheme(base_size = 5))
  
  # Combine ggplot and table
    # combined_plot <- (p / table) | note + plot_layout(heights = c(2, 1), widths = c(3, 1))
    # combined_plot
    
    # Add a title to the table indicating ANOVA was performed
    # table_title <- ggdraw() + 
    #   draw_label("Statistical Test: ANOVA", 
    #              fontface = 'bold', 
    #              size = 10, 
    #              x = 0.5, y = 1, hjust = 0.5)
    
    # Combine plot and table
    combined_plot <- p  / table  # Stacking title and table
    combined_plot + plot_layout(heights = c(2, 1)) + 
      theme(plot.title = element_text(hjust = 1))  # Center the main plot title
  
  ggsave(paste(measure,"diversity.png",sep=""), width = 8, height = 7, dpi = 300, bg = "white")
} 

setwd("~/Documents/CHUM_git/figures/samuel/alpha_diversity/")

#Shannon index
alpha_diversity_graph("Shannon", data = richness_data) 

#Simpson index
alpha_diversity_graph("Simpson", data = richness_data) 

#Inverted Simpson index
alpha_diversity_graph("InvSimpson", data = richness_data) 

#Chao1 index
alpha_diversity_graph("Chao1", data = richness_data)

#saving data so that others can use it
write.csv(richness_data[,-c(1,3)], row.names = TRUE, col.names = 
            TRUE, "~/Documents/CHUM_git/figures/samuel/alpha_diversity/data/alpha_diversity_measures.csv")
}

#Beta diversity analysis
{
#testing the weighted unifrac method
wunifrac_dist <- phyloseq::distance(ps_samuel, method = "wunifrac")

#testing bray-curtis method
bray_dist <- phyloseq::distance(ps_samuel, method = "bray")

pcoaGraphAndStat <- function(dist, ps){
  
  # Perform PCoA
  pcoa_results <- ordinate(ps, method = "PCoA", distance = dist)
  
  #ordination plot
  p <- plot_ordination(ps, pcoa_results, color= "treatment", shape = "genotype") + 
    theme_classic() +
    theme(strip.background = element_blank())+
    stat_ellipse(aes(group = gg_group),      # Add ellipses grouping points by genotype
                 type = "t",  # t-distribution for better fit
                 level = 0.95,  # Confidence level for the ellipse
                 geom = "polygon", alpha = 0)+
    
    labs(title = "PCoA of Bray-Curtis distance matrix") +
    
    theme(
      plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
      axis.title.x = element_text(size = 12),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
      axis.text.x = element_text(size = 12),  # Adjust x-axis tick label font size
      axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
      legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
      legend.text = element_text(size = 12),  # Adjust legend font size
      panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
  
  
  metadata <- data.frame(sample_data(ps))
  groups = levels(sample_data(ps_samuel)$gg_group)
  
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      
      subset_metadata <- metadata[metadata$gg_group %in% c(groups[i],groups[j]),]
      subset_dist <- as.matrix(wunifrac_dist)[subset_metadata$sample_id,subset_metadata$sample_id]
      test.adonis <- adonis(subset_dist ~ gg_group, data = subset_metadata)
      test.adonis <- as.data.frame(test.adonis$aov.tab)
      print(p)
      print(paste("Comparing",groups[i],"and",groups[j], sep = " "))
      print(test.adonis)

    }}
  # cbn <- combn(x=unique(metadata$body.site), m = 2)
  # p <- c()
  # 
  # for(i in 1:ncol(cbn)){
  #   ps.subs <- subset_samples(ps.rarefied, body.site %in% cbn[,i])
  #   metadata_sub <- data.frame(sample_data(ps.subs))
  #   permanova_pairwise <- adonis(phyloseq::distance(ps.subs, method = "bray") ~ body.site, 
  #                                data = metadata_sub)
  #   p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
  # }
  # 
  # p.adj <- p.adjust(p, method = "BH")
  # p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
  # p.table
  
}

pcoaGraphAndStat(bray_dist, ps_samuel)
ggsave("~/Documents/CHUM_git/figures/samuel/beta_diversity/bray_curtis_samuel.png", dpi = 300,  width = 9, height = 7, bg = "white")

pcoaGraphAndStat(wunifrac_dist, ps_samuel)
ggsave("~/Documents/CHUM_git/figures/samuel/beta_diversity/weighted_unifrac_samuel.png", dpi = 300,  width = 9, height = 7, bg = "white")
}

#phyla distribution using relative abundance
{
#save gg_group levels
gg_grouping = levels(sample_data(ps_samuel)$gg_group)

#create ps object with relative abundance for only phyla
phyla_samuel <- ps_samuel %>%
  tax_glom("Phylum") %>%
  transform_sample_counts(function(x) {x * 100 / sum(x)})

phyla_df = cbind(otu_table(phyla_samuel), sample_data(phyla_samuel))
new_colnames <- as.character(tax_table(phyla_samuel)[colnames(phyla_df[1:6]),"Phylum"])
colnames(phyla_df)[1:6] = new_colnames

# Convert the dataframe from wide to long format
phyla_df <- phyla_df %>%
  pivot_longer(cols = 1:6, # Specify the range of columns representing the phyla
               names_to = "phyla",                 # This will create a new column 'phyla' with the phyla names
               values_to = "relab") # This will create a new column 'abundance' with their values


# Aggregate the data by phyla and gg_group, then calculate percentages
phyla_df_agg <- phyla_df %>%
  group_by(gg_group, phyla) %>%
  summarise(relab = sum(relab)) %>%
  mutate(percentage = relab / sum(relab) * 100)


#add log transformation of the data
phyla_df$log_relab <- log(phyla_df$relab + 1)

#check for homogeneity of variances with non log transformed and log transformed data
print(leveneTest(relab~gg_group*phyla, data = phyla_df)) #significantly different
print(leveneTest(log_relab~gg_group*phyla, data = phyla_df)) #significantly different

#check for normality across groups
for(group in gg_grouping){
  for(phylum in new_colnames){
    stat_df <- phyla_df[phyla_df$gg_group == group & phyla_df$phyla == phylum, ]
    print(shapiro.test(stat_df$log_relab))
  }
} #log transformed or no, many of them have the data not normally distributed

#General linear model = makes R crash
{
library(lme4)

phyla_df$sample_id <- as.factor(phyla_df$sample_id)

glmm_model <- glmer(relab ~ gg_group + phyla (1 | sample_id), family = gaussian(link = "log"), data = phyla_df)
summary(glmm_model)
}


phylaPieChart <- function(df){
  
  # Aggregate the data by phyla and gg_group, then calculate percentages
  df_agg <- df %>%
    group_by(gg_group, phyla) %>%
    summarise(relab = sum(relab)) %>%
    mutate(percentage = relab / sum(relab) * 100)
  
  ggplot(df_agg, aes(x = "", y = relab, fill = phyla)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    theme_void() +
    
    # Display percentages on the pie chart
    geom_text(
      aes(label = paste0(round(percentage, 2), "%")), 
      position = position_stack(vjust = 0.5)  # Adjust position in the middle of the slices
    ) +
    facet_wrap(~ gg_group)   # Create separate pie charts for each gg_group
    #scale_fill_brewer(palette = "Set3")  # Optional: Change color palette
  
  
  #statistics
  #shapiro.test(phyla_df)
  #testing for normality
  
  
}

phylaPieChart(phyla_df)


  
ggplot(phyla_df, aes(log_relab)) +
  geom_histogram(binwidth = 1) +
  facet_grid(phyla ~ gg_group)   # Create separate pie charts for each gg_group
  #scale_fill_brewer(palette = "Set3")  # Optional: Change color palette
  

phylaDistribution(phyla_df)


stat_df <- df[df$gg_group == group & df$phyla == phylum, ]
stat_df




#Trying Thibault's package
library(StackbarExtended)

#put cage as factor
sample_data(ps_samuel)$cage = as.factor(sample_data(ps_samuel)$cage)

subset_ps <- subset_samples(ps_samuel, gg_group == c("Wt:Vehicle","Wt:Putrescine"))

my_plot <- plot_microbiota(
  ps_object = ps_samuel,
  exp_group = 'gg_group',
  sample_name = 'sample_id',
  hues = c("Purples", "Blues", "Greens", "Oranges"),
  differential_analysis = T,
  sig_lab = T,
  fdr_threshold = 0.05,
  main_level = "Phylum",
  sub_level = "Family"
)

sample_data(subset_ps)


View(tax_table(ps_samuel))

print(my_plot$plot)
print(my_plot$significant_table_main)
print(my_plot$significant_table_sub)





#Differential expression analysis
#sample_data(ps_samuel)$genotype <- factor(sample_data(ps_samuel)$genotype, levels = c("Wt", "IL-22ra1-/-"))  # Wt as reference
#sample_data(ps_samuel)$treatment <- factor(sample_data(ps_samuel)$treatment, levels = c("DSS + EcNC101 + Vehicle", "DSS + EcNC101 + Putrescine"))  # Vehicle as reference
sample_data(ps_samuel)$gg_group <- factor(sample_data(ps_samuel)$gg_group, levels = c("Wt:Vehicle", "Wt:Putrescine", "IL-22ra1-/-:Vehicle", "IL-22ra1-/-:Putrescine"))  # Vehicle as reference

deseq_samuel <- phyloseq_to_deseq2(ps_samuel, ~ gg_group) 

deseq_samuel = DESeq(deseq_samuel, test="Wald", fitType = "parametric")

resultsNames(deseq_samuel)

#Phyla
# Pairwise comparisons
comparison_1 <- results(deseq_samuel, contrast = c("gg_group", "Wt:Vehicle", "Wt:Putrescine"))
comparison_2 <- results(deseq_samuel, contrast = c("gg_group", "Wt:Vehicle", "IL-22ra1-/-:Vehicle"))
comparison_3 <- results(deseq_samuel, contrast = c("gg_group", "Wt:Putrescine", "IL-22ra1-/-:Putrescine"))
comparison_4 <- results(deseq_samuel, contrast = c("gg_group", "IL-22ra1-/-:Vehicle", "IL-22ra1-/-:Putrescine"))

#Select for phylum only
comparison_1 <- cbind(as(comparison_1, "data.frame"), as(tax_table(ps_samuel)[rownames(comparison_1), "Phylum"], "matrix"))
comparison_2 <- cbind(as(comparison_2, "data.frame"), as(tax_table(ps_samuel)[rownames(comparison_2), "Phylum"], "matrix"))
comparison_3 <- cbind(as(comparison_3, "data.frame"), as(tax_table(ps_samuel)[rownames(comparison_3), "Phylum"], "matrix"))
comparison_4 <- cbind(as(comparison_4, "data.frame"), as(tax_table(ps_samuel)[rownames(comparison_4), "Phylum"], "matrix"))



#Aggregate p-values
aggregatePValues <- function(comparison){
  
  #replace NA pvalues by 1
  comparison$padj[is.na(comparison$padj)] <- 1
  
  #aggregate p-values
  sigtab = aggregate(comparison$log2FoldChange, by = list(comparison$Phylum), FUN = mean)
  sigtab$PValue <- aggregate(comparison$padj, by = list(comparison$Phylum), FUN = mean)$x
  return(sigtab)
}

sigtab_1 = aggregatePValues(comparison_1)
sigtab_2 = aggregatePValues(comparison_2)
sigtab_3 = aggregatePValues(comparison_3)
sigtab_4 = aggregatePValues(comparison_4)




#Differential analysis for species
#old portion of code
{
  #converting ps object to deseq object
  deseq_samuel <- phyloseq_to_deseq2(ps_samuel, ~genotype+treatment) 
  deseq_claire <- phyloseq_to_deseq2(ps_claire, ~diet*week) 
  
  deseq_claire = DESeq(deseq_claire, test="Wald", fitType = "parametric")
  
  # Compare diet 500 vs. diet 50 at week 3
  res_week3 <- results(deseq_claire, contrast = c("diet", "500", "50"), name = "week3")
  
  # Compare diet 500 vs. diet 50 at week 8
  res_week8 <- results(deseq_claire, contrast = list(c("diet500.week8")))
  
  # Compare diet 500 vs. diet 50 at week 10
  res_week10 <- results(deseq_claire, contrast = list(c("diet500.week10")))
  
  # Compare diet 500 vs. diet 50 at week 14
  res_week14 <- results(deseq_claire, contrast = list(c("diet500.week14")))
  
  #function that transforms resutls into significance table (needs ps object, results and alpha as inputs)
  sigTable <- function(ps, res, alpha = 0.01, week){
    
    #investigating results
    sigtab = res[which(res$padj < alpha), ]
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
    
    #Extract abundance data for significant ASVs
    asv_abundance <- otu_table(ps)[ ,rownames(sigtab)]
    abundance_metadata <- cbind(as.data.frame(asv_abundance), as.data.frame(sample_data(ps)))
    abundance_metadata <- abundance_metadata[abundance_metadata$week == week,]
    
    return(abundance_metadata)
  }
  
  test1 = sigTable(ps_claire, res_week3, week = 3) #this is where I observe counts for ASV25
  test2 = sigTable(ps_claire, res_week8, week = 8)
  test3 = sigTable(ps_claire, res_week10, week = 10)
  test4 = sigTable(ps_claire, res_week14, week = 14)
}

relAbundanceDietForTimepoint <- function(ps, alpha = 0.01, week){
  
  ds <- phyloseq_to_deseq2(ps, ~diet) 
  
  ds = DESeq(ds, test="Wald", fitType = "parametric")
  
  res <- results(ds)
  
  #investigating results
  sigtab = res[which(res$padj < 0.01), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  print(sigtab)
  
  #Extract abundance data for significant ASVs
  asv_abundance <- otu_table(ps)[ ,rownames(sigtab)]
  abundance_metadata <- cbind(as.data.frame(asv_abundance), as.data.frame(sample_data(ps)))
  
  #creating new directory
  dir_path <- paste("~/Documents/CHUM_git/figures/relab/week", week, "_claire", sep = "")
  
  # Check if the directory exists, and if not, create it
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Directory created: ", dir_path)
  } else {
    message("Directory already exists: ", dir_path)
  }
  
  #define colors for graph
  diet_colors <- c("50" = "blue", "500" = "red")
  
  for(asv in rownames(sigtab)){
    
    if(isFALSE(is.na(sigtab[asv, "Species"]))){
      
      # Extract the first letter of the genus and append "."
      genus_initial <- paste0(substr(sigtab[asv, "Genus"], 1, 1), ".")
      
      # Define the p-value
      p_value <- sigtab[asv, "padj"]
      
      ggplot(data = abundance_metadata, aes(x = diet, y = abundance_metadata[[asv]], color = diet))+
        geom_point(size = 3) +
        #scale_color_manual(values = diet_colors) +
        labs(title = paste(genus_initial, sigtab[asv, "Species"], "at week", week, sep = " "), y = paste(genus_initial, sigtab[asv, "Species"], "/ 16S", sep = " ")) +
        scale_color_discrete(limits = c("500", "50"))+
        scale_x_discrete(limits = c("50", "500"))+
        stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..-0.25, yend=..y..), color = "black", linewidth =1)+ #adding horizontal bars representing means
        stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..+0.25, yend=..y..), color = "black", linewidth =1)+
        stat_summary(fun.data="mean_cl_normal", geom="errorbar", aes(color = diet), width=0.2, size = 0.7) + #adding SEM error bars
        #ylim(min(data[[measure]])-0.25,max(data[[measure]])+0.25)+
        
        # Add significance bar
        geom_signif(comparisons = list(c("50", "500")),
                    annotations = paste0("p = ", format(p_value, digits = 2, scientific = TRUE)),  # Display p-value
                    y_position = max(abundance_metadata[[asv]], na.rm = TRUE) + 0.1, 
                    tip_length = 0.02, 
                    vjust = 0.5,
                    size = 1.2,  # Make the bar wider
                    color = "black") +
        
        
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
      ggsave(paste(dir_path, "/", sigtab[asv, "Species"], "_wk", week, "_relab.png",sep=""), width = 8, height = 8, dpi = 300, bg = "white")
      
    }
    else{print("No species found")}
    
    
  }
  
  return(abundance_metadata)
}

{
  ps_claire3 <- subset_samples(ps_claire, week == 3)
  ps_claire8 <- subset_samples(ps_claire, week == 8)
  ps_claire10 <- subset_samples(ps_claire, week == 10)
  ps_claire14 <- subset_samples(ps_claire, week == 14)
}
relab_claire3 = relAbundanceDietForTimepoint(ps_claire3, week = 3)
relab_claire8 = relAbundanceDietForTimepoint(ps_claire8, week = 8)
relab_claire10 = relAbundanceDietForTimepoint(ps_claire10, week = 10)
relab_claire14 = relAbundanceDietForTimepoint(ps_claire14, week = 14)

#For Samuel
#Calculating relative abundance 
relab_samuel <- transform_sample_counts(ps_samuel, function(x) x / sum(x))

deseq_samuel <- phyloseq_to_deseq2(ps_samuel, ~genotype+treatment + genotype:treatment) 

deseq_samuel = DESeq(deseq_samuel, test="Wald", fitType = "parametric")

resultsNames(deseq_samuel)

res_il22_putrescine_vs_vehicle <- results(deseq_samuel, 
                                          contrast = c("treatment","treatment_DSS...EcNC101...Vehicle", "DSS...EcNC101...Putrescine"))
# Assuming 'dds' is your DESeqDataSet object and 'condition' is a column in colData
dds_subset <- subset(dds, colData(dds)$condition == "treatmentA")


res <- results(deseq_samuel)

res_genotype <- results(deseq_samuel, name="genotype_Wt_vs_IL.22ra1...")
sigtab = cbind(as(res_genotype, "data.frame"), as(tax_table(ps_samuel)[rownames(res_genotype), ], "matrix"))
print(sigtab)
res_treatment <- results(deseq_samuel, name="treatment_DSS...EcNC101...Vehicle_vs_DSS...EcNC101...Putrescine")
res_interaction <- results(deseq_samuel, name="genotypeWt.treatmentDSS...EcNC101...Vehicle")
res_wt_putrescine_vs_vehicle <- results(deseq_samuel, contrast=list(c("genotypeWt.treatmentDSS...EcNC101...Vehicle", "treatment_DSS...EcNC101...Vehicle_vs_DSS...EcNC101...Putrescine")))
test <- results(deseq_samuel, contrast=list(c("Intercept", "treatment_DSS...EcNC101...Vehicle_vs_DSS...EcNC101...Putrescine")))
res_wt_putrescine_vs_vehicle <- results(deseq_samuel, contrast=c("treatment", "DSS...EcNC101...Putrescine", "DSS...EcNC101...Vehicle"), subset=genotype=="Wt")
res_wt_putrescine_vs_vehicle <- results(deseq_samuel, list(c("treatment_DSS...EcNC101...Vehicle_vs_DSS...EcNC101...Putrescine", "genotypeWt.treatmentDSS...EcNC101...Vehicle")))
res_il22_putrescine_vs_vehicle <- results(deseq_samuel, 
                                          contrast=c("treatment_DSS...EcNC101...Vehicle_vs_DSS...EcNC101...Putrescine", "genotypeIL.22ra1.treatmentDSS...EcNC101...Vehicle"))
res_il22_putrescine_vs_vehicle <- results(deseq_samuel, contrast=c("treatment", "DSS...EcNC101...Putrescine", "DSS...EcNC101...Vehicle"))


# Get results for the full model including interaction
res_full <- results(deseq_samuel, contrast=c("treatment", "DSS...EcNC101...Putrescine", "DSS...EcNC101...Vehicle"))

# Filter for IL.22ra1 genotype
res_il22 <- res_full[which(rownames(res_full) %in% rownames(subset(ps_samuel, genotype == "IL.22ra1"))), ]

res_interaction <- results(deseq_samuel, contrast=list(c("genotype_il22-.treatment_putrescine", "genotype_wt.treatment_vehicle")))

sigTab <- function(res, ps){
  sigtab = cbind(as(res, "data.frame"), as(tax_table(ps)[rownames(res), ], "matrix"))
  print(sigtab)
}


sigTab(test, ps_samuel)

#investigating results
sigtab = res[which(res$padj < 0.01), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_samuel)[rownames(sigtab), ], "matrix"))
print(sigtab)
write.csv(sigtab, "~/Documents/CHUM_git/figures/samuel/differential_analysis/data/significance_table.csv", row.names = TRUE, col.names = 
            TRUE)

#Extract abundance data for significant ASVs
asv_abundance <- otu_table(ps_samuel)[ ,rownames(sigtab)]
abundance_metadata <- cbind(as.data.frame(asv_abundance), as.data.frame(sample_data(ps_samuel)))
write.csv(abundance_metadata, "~/Documents/CHUM_git/figures/samuel/differential_analysis/data/asv_counts_and_metadata.csv", row.names = TRUE, col.names = 
            TRUE)

#creating new directory
dir_path <- "~/Documents/CHUM_git/figures/samuel/differential_analysis/"

# Check if the directory exists, and if not, create it
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
  message("Directory created: ", dir_path)
} else {
  message("Directory already exists: ", dir_path)
}

relAbundance2Var <- function(sigtab, abundance_metadata){
  
  for(asv in rownames(sigtab)){
    
    if(isFALSE(is.na(sigtab[asv, "Species"]))){
      
      # Extract the first letter of the genus and append "."
      genus_initial <- paste0(substr(sigtab[asv, "Genus"], 1, 1), ".")
      
      # Define the p-value
      p_value <- sigtab[asv, "padj"]
      
      p <- ggplot(data = abundance_metadata, aes(x = genotype, y = abundance_metadata[[asv]], color = treatment)) +
        geom_point(size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) + 
        
        # Error bars
        stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
                     aes(color = treatment),
                     width = 0.2, size = 0.7,
                     position = position_dodge(-0.75)) +
        
        #Mean lines
        stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
                     aes(ymin = ..y.., ymax = ..y.., group = treatment),
                     color = "black", linewidth = 0.5, width = 0.5,
                     position = position_dodge(-0.75))+
        
        
        labs(title = paste(genus_initial, sigtab[asv, "Species"], sep = " "),
             y = paste(genus_initial, sigtab[asv, "Species"], "/ 16S", sep = " ")) +
        scale_x_discrete(limits = c("Wt", "IL-22ra1-/-")) +
        scale_color_discrete(limits = c("DSS + EcNC101 + Vehicle", "DSS + EcNC101 + Putrescine")) +

      
      # Add significance bar
      # geom_signif(comparisons = list(c("50", "500")),
      #             annotations = paste0("p = ", format(p_value, digits = 2, scientific = TRUE)),  # Display p-value
      #             y_position = max(abundance_metadata[[asv]], na.rm = TRUE) + 0.1, 
      #             tip_length = 0.02, 
      #             vjust = 0.5,
      #             size = 1.2,  # Make the bar wider
      #             color = "black") +
      
      
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
      #ggsave(paste(dir_path, "/", sigtab[asv, "Species"], "_relab.png",sep=""), width = 8, height = 8, dpi = 300, bg = "white")
      
    }
    else{print("No species found")}
    
    
  }
  print(p)
  #return(abundance_metadata)
}

relab_samuel = relAbundance2Var(sigtab, abundance_metadata)
#old code playing with deseq and graphs
{
  #perform deseq2 analysis = function that takes deseq object and associated ps object as input
  #returns a significance table
  deseq2Analysis <- function(ps, ds, alpha = 0.01){
    ds = DESeq(ds, test="Wald", fitType = "parametric")
    
    #investigating results
    res = results(ds, cooksCutoff = FALSE)
    alpha = alpha #significance threshold 
    sigtab = res[which(res$padj < alpha), ]
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
    return(sigtab)
  }
  
  sigtab_claire = deseq2Analysis(ps_claire, deseq_claire)
  sigtab_samuel = deseq2Analysis(ps_samuel, deseq_samuel)
  
  
  #looking at otus that were significantly different
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  # Phylum order
  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  # Genus order
  x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
  ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
} 

}



#Claire stuff
#Estinate richness measures for dataset
richness_data <- estimate_richness(ps_claire, measures = c("Shannon", "Simpson","InvSimpson","Chao1","Observed"))

#Add sample metadata to richness dataframe
richness_data <- cbind(as.data.frame(sample_data(ps_claire)), richness_data)

#saving data so that others can use it
write.csv(richness_data, row.names = FALSE, col.names = 
            TRUE, "~/Documents/CHUM_git/figures/claire/alpha_diversity/data/alpha_diversity_measures.csv")

#Testing new pipeline for alpha diversity
source(file = "~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/alpha_diversity_graphs_and_stats.R")
alphaDiversityGgGroup(ps_samuel, path = "~/Documents/CHUM_git/pipeline_tests/", group = "gg_group")

#Beta diversity analysis for different timepoints. You must provide a filtered ps object, the timeVariable and the varToCompare (present in sample_data)
source(file = "~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/beta_diversity_graphs_and_stats.R")

#For Claire
betaDiversityTimepoint(ps_claire, "week", "diet", distMethod = "bray", customColors = c('blue','red'), "~/Documents/CHUM_git/figures/claire/test/")
betaDiversityTimepoint(ps_claire, "week", "diet", distMethod = "wunifrac", customColors = c('blue','red'), "~/Documents/CHUM_git/figures/claire/test/")

#For Samuel
pairs <- list(list("Wt:Vehicle","Wt:Putrescine"), list("IL-22ra1-/-:Vehicle","IL-22ra1-/-:Putrescine"), list("Wt:Vehicle","IL-22ra1-/-:Vehicle"), list("Wt:Putrescine","IL-22ra1-/-:Putrescine"))
customColors = list(list('black','#A22004'), list("#AB8F23","#04208D"), list("black","#AB8F23"),list("#A22004","#04208D"))
betaDiversityPairwise(ps_samuel, "gg_group", pairs, "bray", customColors, "~/Documents/CHUM_git/figures/samuel/test/")
betaDiversityPairwise(ps_samuel, "gg_group", pairs, "wunifrac", customColors, "~/Documents/CHUM_git/figures/samuel/test/")






#Testing stuff with relative abundance and differntial analysis or something like that
deseq_samuel <- phyloseq_to_deseq2(ps_samuel, ~ gg_group) 
deseq_samuel <- DESeq(deseq_samuel, test="Wald", fitType = "parametric")
resultsNames(deseq_samuel)

source(file = "~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/relab_analysis_graphs_and_stats.R")
pairs <- list(list("Wt:Vehicle","Wt:Putrescine"), list("IL-22ra1-/-:Vehicle","IL-22ra1-/-:Putrescine"), list("Wt:Vehicle","IL-22ra1-/-:Vehicle"), list("Wt:Putrescine","IL-22ra1-/-:Putrescine"))
customColors = list(list('black','#A22004'), list("#AB8F23","#04208D"), list("black","#AB8F23"),list("#A22004","#04208D"))
relabSpeciesPairwise(ps_samuel, deseq_samuel, measure = "log2fold", "gg_group", pairs, threshold = 0.01, customColors, "~/Documents/CHUM_git/figures/samuel/test_relab/")
  
#deseq_claire <- phyloseq_to_deseq2(ps_claire, ~ week*diet) 
deseq_claire <- phyloseq_to_deseq2(ps_claire, ~ week + diet:week)
deseq_claire <- DESeq(deseq_claire, test="Wald", fitType = "parametric")
resultsNames(deseq_claire)

#Testing if this can work at any taxonomic level
customColors = c("blue","red")
#At species level
relabTimeline(ps_claire, deseq_claire, measure = "log2fold", "week", "diet", "Species", threshold = 0.01, customColors,  "~/Documents/CHUM_git/figures/claire/test_timeline/")

#At other taxonomic levels
taxonomicLevels <- c("Genus","Family","Order","Class","Phylum")
taxonomicLevels <- c("Phylum")
for(txnLevel in taxonomicLevels){
  
  #Creates ps subset for taxonomical level of interest
  ps_subset <- tax_glom(ps_claire, taxrank = txnLevel)
  deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ week + diet:week) 
  deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
  relabTimeline(ps_subset, deseq_subset, measure = "log2fold", "week", "diet", txnLevel, threshold = 1, customColors, "~/Documents/CHUM_git/figures/claire/test_timeline/")
}





#Putting my whole life in perspective....
deseq_claire <- phyloseq_to_deseq2(ps_claire, ~ diet*week) 
deseq_claire <- DESeq(deseq_claire, test="Wald", fitType = "parametric")
resultsNames(deseq_claire)
plotDispEsts(deseq_claire) 
res_week3 <- results(deseq_claire, name = "diet_500_vs_50")
print(res_week3)
res_week8 <- results(deseq_claire, contrast = list("diet_500_vs_50", "week_8_vs_3"))
print(res_week8)
res_week10 <- results(deseq_claire, contrast = list("diet_500_vs_50", "week_10_vs_3"))
print(res_week10)
res_week14 <- results(deseq_claire, contrast = list("diet_500_vs_50", "week_10_vs_3"))
print(res_week14)









#Revising technique for Samuel's data
levels(colData(deseq_samuel)$genotype)
levels(colData(deseq_samuel)$treatment)

#Setting "Wt" as the baseline for genotype
colData(deseq_samuel)$genotype <- relevel(colData(deseq_samuel)$genotype, ref="Wt")

#Setting "Vehicle" as the baseline for treatment
colData(deseq_samuel)$treatment <- relevel(colData(deseq_samuel)$treatment, ref="Vehicle")

deseq_samuel <- phyloseq_to_deseq2(ps_samuel, ~ genotype + treatment + genotype:treatment) 
deseq_samuel <- DESeq(deseq_samuel, test="Wald", fitType = "parametric")
resultsNames(deseq_samuel)

source(file = "~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/relab_analysis_graphs_and_stats.R")
customColors = list('black','#A22004',"#AB8F23","#04208D")
#At species level
relabGroups(ps_samuel, deseq_samuel, measure = "log2fold", "gg_group", taxa = "Species", displayPvalue = FALSE, threshold = 0.01, customColors, "~/Documents/CHUM_git/figures/samuel/rel_ab_groups/")

#At other taxonomic levels
taxonomicLevels <- c("Genus","Family","Order","Class","Phylum")
# taxonomicLevels <- c("Phylum")
for(txnLevel in taxonomicLevels){
  
  #Creates ps subset for taxonomical level of interest
  ps_subset <- tax_glom(ps_samuel, taxrank = txnLevel)
  deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ genotype + treatment + genotype:treatment) 
  deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
  relabGroups(ps_subset, deseq_subset, measure = "log2fold", "gg_group", taxa = txnLevel, displayPvalue = FALSE, threshold = 0.01, customColors, "~/Documents/CHUM_git/figures/samuel/rel_ab_groups/")
}

