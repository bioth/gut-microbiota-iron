library(phyloseq)
library(Biostrings)
library(ggplot2)
library(dplyr)
#library(metamicrobiomeR)
#library(ggrepel)
library(car)
#library(PMCMR)
#library(nlme)
#library(ggsignif)
library(DESeq2)

#for microbiota 17 Samuel's subset
setwd("D:/CHUM_git/Microbiota_17_Samuel/")
asv_table <- read.csv("seqtab.nochim_run1.csv", sep = ";")

#loading metadata of interest
df <- as.data.frame(readxl::read_excel("../Microbiota_17/metadata/Samuel_Meta_data_grouping.xlsx"))

#put sample_id as rownames
rownames(df) <- df[,2]

#load taxonomical assignments
taxa <- as.matrix(read.csv("taxa_annotation.csv", sep = ";"))

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

#Calculate the proportion of samples containing each ASV
otu_data <- as(otu_table(ps), "matrix")
prevalence <- apply(otu_data, 2, function(x) {sum(x > 0)})
prevalence <- (prevalence / nsamples(ps))*100

# Define a prevalence threshold (e.g., 10% of samples)
threshold <- 10

# Keep only ASVs present in more than the threshold proportion of samples
keep_taxa <- prevalence > threshold
ps_flt <- prune_taxa(keep_taxa, ps)

# Trying DESeq2 package
#converting phyloseq object into deseq2 object and perform statistical testing
deseq2_ps <- phyloseq_to_deseq2(ps_flt, ~Genotype+Groups)
a <- DESeq(deseq2_ps, test = "Wald", fitType = "parametric")

#store results and cut for alpha values < 0.01
res <- results(a, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_flt)[rownames(sigtab), ], "matrix"))
head(sigtab)

#first number of OTUs for which they were signficant differences (so small rn)
dim(sigtab)

#You can plot this
#Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)

#Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

#Phyla and genuses
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))








# #Creates dataframe with relative abundance data
ps_rel_ab <- transform_sample_counts(ps, function(x) x/sum(x))

#Transform into dataframe
ps_rel_ab <- psmelt(ps_rel_ab)

#function to add SEM
mean_cl_normal <- function(x, mult = 1.96) { #mult is 1.96 for a 95% confidence interval
  # Calculate the mean of the input vector x
  mean_val <- mean(x, na.rm = TRUE)
  
  # Calculate the standard error of the mean
  se_val <- sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
  
  # Return a data frame with the mean (y), and the lower (ymin) and upper (ymax) bounds
  data.frame(y = mean_val, ymin = mean_val - mult * se_val, ymax = mean_val + mult * se_val)
}

##Visualize mean alpha diversity
# Extract data and compute means
richness_data <- estimate_richness(ps_flt, measures = c("Shannon", "Simpson","Chao1"))

# Extract sample data
sample_data <- as.data.frame(sample_data(ps_flt))

# Combine richness data with sample data
richness_data <- cbind(sample_data, richness_data)

# Compute mean values for each combination of diet, and treatment
richness_data <- richness_data %>%
  mutate(gg_group = paste(Genotype, Groups, sep = "_")) 

setwd("figures")

#Alpha diversity graphs
alpha_diversity_graph <- function(measure, data){
  
  #Choosing appropriate statistical tests
  anova_result <- aov(data[[measure]] ~ Genotype * Groups, data = data)
  
  # Extract the residuals
  residuals <- residuals(anova_result)
  
  # Shapiro-Wilk test for normality
  print("Checking for normality")
  shapiro_test <- shapiro.test(residuals)
  print(shapiro_test)
  
  # Verifying homoscedasticity assumption
  print("Checking for homoscedasticity")
  levene_test <- leveneTest(data[[measure]] ~ Genotype * Groups, data = data)
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
    View(tukey_df_all)
    
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
  
  # Graph
  ggplot(data, aes(x = gg_group, y = data[[measure]], color = Groups, shape = Genotype)) +
    geom_point(size = 3) +
    labs(title = paste("Mean",measure,"diversity"), y = paste(measure, "index")) +
    scale_shape_discrete(limits = c("Wt","IL-22ra1-/-"))+
    scale_color_discrete(limits = c("DSS + EcNC101 + Vehicle", "DSS + EcNC101 + Putrescine"))+
    scale_x_discrete(limits = c("Wt_DSS + EcNC101 + Vehicle","Wt_DSS + EcNC101 + Putrescine","IL-22ra1-/-_DSS + EcNC101 + Vehicle","IL-22ra1-/-_DSS + EcNC101 + Putrescine"))+
    stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..-0.25, yend=..y..), color = "black", linewidth =1)+ #adding horizontal bars representing means
    stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..+0.25, yend=..y..), color = "black", linewidth =1)+ 
    stat_summary(fun.data="mean_cl_normal", geom="errorbar", color="red", width=0.2, size = 0.7) + #adding SEM error bars
    ylim(min(data[[measure]])-0.25,max(data[[measure]])+0.25)+
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
  

  ggsave(paste(measure,"diversity.png",sep=""), width = 8, height = 8, dpi = 300, bg = "white")
}

#Shannon index
alpha_diversity_graph("Shannon", data = richness_data) 

#Simpson index
alpha_diversity_graph("Simpson", data = richness_data) 

#Chao1 index
alpha_diversity_graph("Chao1", data = richness_data) 

#Beta diversity



#going for bacterial relative abundance piecharts for the phylum distributions
# Extract the taxonomy table from the phyloseq object
tax_table <- tax_table(ps_flt)

# Create a data frame from the phyloseq object
ps_df <- psmelt(ps_flt)


ps_df <- taxa.meansdn(taxtab = tax_table, sumvar = "bf", groupvar = )

# Group by Phylum and calculate the relative abundance
ps_df <- ps_df %>%
  mutate(gg_group = paste(Genotype, Groups, sep = "_"))%>%
  group_by(gg_group, Phylum) %>%
  summarize(Abundance = sum(Abundance)) %>%
  group_by(gg_group) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance) * 100)

# Step 2: Identify phyla with relative abundance < 1% across all gg_groups
low_abundance_phyla <- ps_df %>%
  group_by(Phylum) %>%
  summarize(Max_Relative_Abundance = max(Relative_Abundance)) %>%
  filter(Max_Relative_Abundance < 1) %>%
  pull(Phylum)

# Step 3: Group low abundance phyla into "Others"
ps_df <- ps_df %>%
  mutate(Phylum = ifelse(Phylum %in% low_abundance_phyla, "Others", Phylum)) %>%
  group_by(gg_group, Phylum) %>%
  summarize(Relative_Abundance = sum(Relative_Abundance))

library(ggforce)
library(ggre)

#Define groups orders for graphs
graph_groups_order = c("Wt_DSS + EcNC101 + Vehicle","Wt_DSS + EcNC101 + Putrescine",

                                            "IL-22ra1-/-_DSS + EcNC101 + Vehicle","IL-22ra1-/-_DSS + EcNC101 + Putrescine")

# Reorder gg_group based on desired_order
ps_df$gg_group <- factor(ps_df$gg_group, levels = graph_groups_order)

#piechart with bacterial relative abundance by phyla
ggplot(ps_df, aes(x = "", y = Relative_Abundance, fill = Phylum)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  facet_wrap(~ gg_group) +
  theme_void() +
  labs(title = "Bacterial relative abundance by phyla")+
  #Add lines that indicate the percentages
  geom_text(aes(label = ifelse(Relative_Abundance > 1, 
                               paste0(round(Relative_Abundance, 1), "%"), "")), 
            position = position_stack(vjust = 0.5), 
            size = 7, 
            color = "black",
            hjust = 0.5) +  # Adjust hjust value as needed
  theme(
    plot.title = element_text(size = 16, face = "bold", vjust = 5),  # Adjust title font size and style
    legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend font size
  )

ggsave("phyla_distribution.png", width = 10, height = 10, dpi = 300, bg = "white")

# #stats WARNING TEST AND WAY TO PERFORM ANALYSIS HERE MAY NOT BE APPROPRIATED
# anova_result <- aov(Relative_Abundance ~ gg_group * Phylum, data = physeq_df)
# summary(anova_result)
# tukey_results <- TukeyHSD(anova_result)
# print(tukey_results)
# factor_names <- names(tukey_results)
# 
# # Prepare to combine results
# tukey_df_list <- lapply(factor_names, function(factor) {
#   tukey_result <- tukey_results[[factor]]
#   if (nrow(tukey_result) > 0) {  # Check if there are any results
#     tukey_df <- as.data.frame(tukey_result)
#     tukey_df$comparison <- rownames(tukey_df)
#     tukey_df$factor <- factor
#     return(tukey_df)
#   } else {
#     return(NULL)  # Return NULL if no results
#   }
# })
# 
# # Remove NULL elements and combine all results into one data frame
# tukey_df_all <- do.call(rbind, Filter(Negate(is.null), tukey_df_list))
# 
# # Add significance stars based on p-values
# tukey_df_all$significance <- cut(tukey_df_all$`p adj`,
#                                  breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
#                                  labels = c("***", "**", "*", ""),
#                                  include.lowest = TRUE)
































#Bacterial species relative abundance
#isolating species data
ps_species <- tax_glom(ps_flt, taxrank = "Species")

#Calculate number of unique species
ps_species_melt <- psmelt(ps_species)
number_of_species <- length(unique(ps_species_melt$Species))
print(number_of_species)

#change working directory for relative abundances graphs
setwd("relative_abundances/")

#Function that generates a species relative abundance graph for all species present 
speciesRelativeAbundanceGraph <- function(asv = 1, ps){
  
  # Check if asv index is within range
  taxa_table <- tax_table(ps)
  if (asv > nrow(taxa_table)) {
    stop("The ASV index is out of range")
  }
  
  # Isolate the ASV ID for the given species
  species_asv_id <- rownames(taxa_table)[asv]
  
  # Subset the phyloseq object to keep only the selected ASV
  ps_sub <- prune_taxa(species_asv_id, ps)
  
  # Calculate relative abundances across samples
  otu_table_sub <- otu_table(ps_sub)
  ps_rel_ab <- apply(otu_table_sub, 2, function(x) x / sum(x))
  
  # Convert to data frame and add sample metadata
  rel_ab_df <- as.data.frame(ps_rel_ab)
  rel_ab_df <- cbind(sample_data(ps_sub), rel_ab_df)
  
  #Putting names for the species and the plot name
  taxa_table <- as.data.frame(taxa_table)
  species_name <- paste(taxa_table$Genus[asv], taxa_table$Species[asv])
  plot_name <- paste(species_name,"relative_ab_plot.png", sep = "")
  
  #Adding gg_group
  rel_ab_df <- rel_ab_df %>%
    mutate(gg_group = paste(Genotype, Groups, sep = "_")) 
  
  # Convert species_asv_id to symbol
  species_asv_id <- sym(species_asv_id)
  
  # Plotting the results
  ggplot(rel_ab_df, aes(x = gg_group, y = (!!species_asv_id)*100, color = Groups, shape = Genotype)) +
    geom_point(size = 3) +
    labs(title = paste(species_name, "relative abundance"), y = "Relative abundance (%)") +
    scale_shape_discrete(limits = c("Wt","IL-22ra1-/-")) +
    scale_color_discrete(limits = c("DSS + EcNC101 + Vehicle", "DSS + EcNC101 + Putrescine")) +
    scale_x_discrete(limits = c("Wt_DSS + EcNC101 + Vehicle", "Wt_DSS + EcNC101 + Putrescine",
                                "IL-22ra1-/-_DSS + EcNC101 + Vehicle", "IL-22ra1-/-_DSS + EcNC101 + Putrescine")) +
    stat_summary(fun = "mean", geom = "segment", mapping = aes(xend = ..x.. - 0.25, yend = ..y..), color = "black", linewidth = 1) +
    stat_summary(fun = "mean", geom = "segment", mapping = aes(xend = ..x.. + 0.25, yend = ..y..), color = "black", linewidth = 1) +
    stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", color = "red", width = 0.2, size = 0.7) +
    ylim(0, max(rel_ab_df[[species_asv_id]], na.rm = TRUE)*100 + 1) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12),
      panel.grid.major = element_line(color = "gray90", size = 0.5),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 1)
    )
  
  ggsave(plot_name, width = 8, height = 8, dpi = 300, bg = "white")
  
  #Stats
  anova_result <- aov(rel_ab_df[[species_asv_id]] ~ Genotype * Groups, data = rel_ab_df)
  
  # Extract the residuals
  residuals <- residuals(anova_result)
  
  # Shapiro-Wilk test for normality
  shapiro_test <- shapiro.test(residuals)
  print(shapiro_test)
  
  # Verifying homoscedasticity assumption
  levene_test <- leveneTest(rel_ab_df[[species_asv_id]] ~ Genotype * Groups, data = rel_ab_df)
  print(levene_test)
  
  cat("Which test do you want to perform ? (ANOVA: a / GLS: g/ KRUSKAL_WALLIS: k) / WELCH T-TEST: w ")
  
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
    
  }
  if(test_chosen == "g"){
    
    # GLS in R (if variances are significantly different)
    gls_model <- gls(rel_ab_df[[species_asv_id]] ~ Genotype * Groups, data = rel_ab_df, weights = varIdent(form = ~ 1 | Genotype * Groups))
    summary(gls_model)
    
    # dunnet_results <- dunnettT3Test(welch_anova_result)
    # print(dunnet_results)
  }
  if(test_chosen == "w"){
    group1 <- rel_ab_df %>%
      filter(gg_group %in% c("Wt_DSS + EcNC101 + Vehicle", "Wt_DSS + EcNC101 + Putrescine"))
    group2 <- rel_ab_df %>%
      filter(gg_group %in% c("IL-22ra1-/-_DSS + EcNC101 + Vehicle", "IL-22ra1-/-_DSS + EcNC101 + Putrescine"))
    group3 <- rel_ab_df %>%
      filter(gg_group %in% c("Wt_DSS + EcNC101 + Vehicle", "IL-22ra1-/-_DSS + EcNC101 + Vehicle"))
    group4 <- rel_ab_df %>%
      filter(gg_group %in% c("Wt_DSS + EcNC101 + Putrescine", "IL-22ra1-/-_DSS + EcNC101 + Putrescine"))
    t_test_result <- t.test(group1[[species_asv_id]] ~ gg_group, data = group1, var.equal = FALSE)
    print(t_test_result)
    t_test_result <- t.test(group2[[species_asv_id]] ~ gg_group, data = group1, var.equal = FALSE)
    print(t_test_result)
    t_test_result <- t.test(group3[[species_asv_id]] ~ gg_group, data = group1, var.equal = FALSE)
    print(t_test_result)
    t_test_result <- t.test(group4[[species_asv_id]] ~ gg_group, data = group1, var.equal = FALSE)
    print(t_test_result)
  }
  if(test_chosen == "k"){
    group1 <- rel_ab_df %>%
      filter(gg_group %in% c("Wt_DSS + EcNC101 + Vehicle", "Wt_DSS + EcNC101 + Putrescine"))
    group2 <- rel_ab_df %>%
      filter(gg_group %in% c("IL-22ra1-/-_DSS + EcNC101 + Vehicle", "IL-22ra1-/-_DSS + EcNC101 + Putrescine"))
    group3 <- rel_ab_df %>%
      filter(gg_group %in% c("Wt_DSS + EcNC101 + Vehicle", "IL-22ra1-/-_DSS + EcNC101 + Vehicle"))
    group4 <- rel_ab_df %>%
      filter(gg_group %in% c("Wt_DSS + EcNC101 + Putrescine", "IL-22ra1-/-_DSS + EcNC101 + Putrescine"))
    kruskal_result <- kruskal.test(group1[[species_asv_id]] ~ gg_group, data = group1)
    print(kruskal_result)
    kruskal_result <- kruskal.test(group2[[species_asv_id]] ~ gg_group, data = group2)
    print(kruskal_result)
    kruskal_result <- kruskal.test(group3[[species_asv_id]] ~ gg_group, data = group3)
    print(kruskal_result)
    kruskal_result <- kruskal.test(group4[[species_asv_id]] ~ gg_group, data = group4)
    print(kruskal_result)
  }
  print(species_name)
 
}


#For loop to generate all of the graphs
for(i in 1:number_of_species){
  speciesRelativeAbundanceGraph(ps_species, asv = i)
}



speciesRelativeAbundanceGraph(ps_species, asv = 1)





library("DESeq2")

#Convert the subsetted phyloseq object to a DESeq2 dataset
deseq_data <- phyloseq_to_deseq2(ps_sub, ~ Groups + Genotype)  # `group` should be your grouping variable

#Run DESeq2 analysis
deseq_data <- DESeq(deseq_data)





















#throw this one out for Claire's data
taxa.compare(rel_ab_df, propmed.rel = "gamlss", transform = "none", zeroreplace.method = "none", comvar = ASV13, adjustvar,
  personid = rel_ab_df$`Sample ID`,
  longitudinal = "no",
  percent.filter = 0.05,
  relabund.filter = 5e-05,
  p.adjust.method = "fdr"
)




