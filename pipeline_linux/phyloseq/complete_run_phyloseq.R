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
  library(openxlsx)
  library(reshape2)
  library(Hmisc)
  library(plotly) #To plot 3D pcoas
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

#Load custom functions for microbiota analysis
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/alpha_diversity_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/beta_diversity_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/correlation_graphs_and_stats.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/relab_analysis_graphs_and_stats.R")

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

# Transform otu_table of a phyloseq object based on the chosen transformation
transformCounts <- function(ps, transformation = "rel_ab", log_base = 10) {
  if (transformation == "rel_ab") {
    
    # Check orientation of OTU table and transpose if necessary
    if (isFALSE(taxa_are_rows(otu_table(ps)))) {
      otu_table(ps) <- t(otu_table(ps))
      message("OTU table was transposed to have ASVs as rows.")
    } else {
      message("ASVs already as rows, no need to transpose OTU table.")
    }
    
    # Apply prop.table for relative abundance calculation and multiply by 100 for percentage
    otu_matrix <- apply(otu_table(ps), 2, prop.table) * 100
    
    # Ensure that the result is an otu_table object
    otu_table(ps) <- otu_table(otu_matrix, taxa_are_rows = TRUE)
    
  } else if (transformation == "log") {
    
    # Check if there are any zero counts to avoid log(0) issues
    if (any(otu_table(ps) == 0)) {
      message("Warning: Zero values detected in the OTU table. Adding a small constant to avoid log(0).")
      otu_table(ps) <- otu_table(ps) + 1e-6
    }
    
    # Apply log transformation (base 10 by default)
    otu_matrix <- log(otu_table(ps), base = log_base)
    
    # Ensure that the result is an otu_table object
    otu_table(ps) <- otu_table(otu_matrix, taxa_are_rows = TRUE)
  } else {
    stop("Transformation method not recognized. Use 'rel_ab' for relative abundance or 'log' for log transformation.")
  }
  
  return(ps) # Return the transformed phyloseq object
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
metadata$genotype <- gsub(".*IL-22.*", "IL-22ra1-/-", metadata$genotype)

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
#And adds refseq entry to the phyloseq object
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
sum(taxa_sums(ps)) #total number of reads
length(taxa_sums(ps)) #total number of ASVs

#create two separate ps objects for each Claire and Samuel's datasets
ps_samuel <- subset_samples(ps, student == "Samuel")
sum(taxa_sums(ps_samuel))
length(taxa_sums(ps_samuel))

ps_claire <- subset_samples(ps, student == "Claire")
sum(taxa_sums(ps_claire))
length(taxa_sums(ps_claire))

#Filtering should be reserved for differential abundance analysis only
{
#function filtering out ASVs for which they were in total less than a threshold count
ps_samuel <- prune_taxa(taxa_sums(ps_samuel) > 0, ps_samuel) #Actual ASVs from Samuel's data
sum(taxa_sums(ps_samuel))
length(taxa_sums(ps_samuel))

ps_samuel <- prune_taxa(taxa_sums(ps_samuel) > 10, ps_samuel)
sum(taxa_sums(ps_samuel))
length(taxa_sums(ps_samuel))

ps_claire <- prune_taxa(taxa_sums(ps_claire) > 0, ps_claire) #Actual ASVs from Samuel's data
sum(taxa_sums(ps_claire))
length(taxa_sums(ps_claire))

ps_claire <- prune_taxa(taxa_sums(ps_claire) > 10, ps_claire)
sum(taxa_sums(ps_claire))
length(taxa_sums(ps_claire))

#Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
ps_samuel <- prune_taxa(colSums(otu_table(ps_samuel) > 0) >= (0.05 * nsamples(ps_samuel)), ps_samuel)
sum(taxa_sums(ps_samuel))
length(taxa_sums(ps_samuel))

ps_claire <- prune_taxa(colSums(otu_table(ps_claire) > 0) >= (0.05 * nsamples(ps_claire)), ps_claire)
sum(taxa_sums(ps_claire))
length(taxa_sums(ps_claire))
}


#for Samuel's data, put gg_group as factor and define order
sample_data(ps_samuel)$gg_group <- factor(sample_data(ps_samuel)$gg_group, levels = c("Wt:Vehicle", "Wt:Putrescine", "IL-22ra1-/-:Vehicle", "IL-22ra1-/-:Putrescine"))  # Vehicle as reference

###Alpha diversity
#Claire stuff
#Estinate richness measures for dataset
richness_data <- estimate_richness(ps_claire, measures = c("Shannon", "Simpson","InvSimpson","Chao1","Observed"))

#Add sample metadata to richness dataframe
richness_data <- cbind(as.data.frame(sample_data(ps_claire)), richness_data)

#saving data so that others can use it
write.xlsx(richness_data, "~/Documents/CHUM_git/figures/claire/new_alpha_diversity/alpha_diversity_measures.xlsx")

#Samuel alpha diversity
customColors = list('black','#A22004',"#AB8F23","#04208D")
pairs <- list(list("Wt:Vehicle","Wt:Putrescine"), list("IL-22ra1-/-:Vehicle","IL-22ra1-/-:Putrescine"), list("Wt:Vehicle","IL-22ra1-/-:Vehicle"), list("Wt:Putrescine","IL-22ra1-/-:Putrescine"))
alphaDiversityGgGroup2(ps_samuel, path = "~/Documents/CHUM_git/figures/samuel/new_alpha_diversity/", gg_group = "gg_group", customColors = customColors)

###Beta diversity
#Beta diversity analysis for different timepoints. You must provide a filtered ps object, the timeVariable and the varToCompare (present in sample_data)
source(file = "~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/beta_diversity_graphs_and_stats.R")

#For Claire
betaDiversityTimepoint(ps_claire, "week", "diet", distMethod = "bray", customColors = c('blue','red'), font = "Arial", "~/Documents/CHUM_git/figures/claire/beta_diversity/")
betaDiversityTimepoint(ps_claire, "week", "diet", distMethod = "wunifrac", customColors = c('blue','red'), font = "Arial", "~/Documents/CHUM_git/figures/claire/beta_diversity/")


#Samuel beta diversity
pairs <- list(list("Wt:Vehicle","Wt:Putrescine"), list("IL-22ra1-/-:Vehicle","IL-22ra1-/-:Putrescine"), list("Wt:Vehicle","IL-22ra1-/-:Vehicle"), list("Wt:Putrescine","IL-22ra1-/-:Putrescine"))
customColors = list(list('black','#A22004'), list("#AB8F23","#04208D"), list("black","#AB8F23"),list("#A22004","#04208D"))
#Dim represents respectively height and width 
betaDiversityPairwise(ps_samuel, "gg_group", pairs, "bray", customColors, font = "Times New Roman", displayPValue = FALSE, dim = c(5,5), transform = "rel_ab", path = "~/Documents/CHUM_git/figures/filtering/samuel/beta_diversity/")
betaDiversityPairwise(ps_samuel, "gg_group", pairs, "wunifrac", customColors, font = "Times New Roman", displayPValue = FALSE, c(5,5), transform = "log", path = "~/Documents/CHUM_git/figures/filtering/samuel/beta_diversity/")

#For all groups at the same time
customColors = c('black','#A22004',"#AB8F23","#04208D")
betaDiversityAll(ps_samuel, "gg_group", "bray", customColors, font = "Times New Roman", displayPValue = FALSE, dim = c(5,5), transform = "rel_ab", variancePlot = FALSE, path = "~/Documents/CHUM_git/figures/filtering/samuel/beta_diversity/")
betaDiversityAll(ps_samuel, "gg_group", "wunifrac", customColors, font = "Times New Roman", displayPValue = FALSE, dim = c(5,5), transform = "log", variancePlot = FALSE, path = "~/Documents/CHUM_git/figures/filtering/samuel/beta_diversity/")

#To look at 3 PCs (3D representation)
betaDiversityAll(ps_samuel, "gg_group", "bray", customColors, font = "Times New Roman", displayPValue = FALSE, dim = c(5,5), transform = "rel_ab", variancePlot = FALSE, display3PCs = TRUE, path = "~/Documents/CHUM_git/figures/filtering/samuel/beta_diversity/")


###Differential abundance accross different taxonomic levels
#Filtering should be reserved for differential abundance analysis only
{
  #function filtering out ASVs for which they were in total less than a threshold count
  ps_samuel <- prune_taxa(taxa_sums(ps_samuel) > 0, ps_samuel) #Actual ASVs from Samuel's data
  sum(taxa_sums(ps_samuel))
  length(taxa_sums(ps_samuel))
  
  ps_samuel <- prune_taxa(taxa_sums(ps_samuel) > 10, ps_samuel)
  sum(taxa_sums(ps_samuel))
  length(taxa_sums(ps_samuel))
  
  ps_claire <- prune_taxa(taxa_sums(ps_claire) > 0, ps_claire) #Actual ASVs from Samuel's data
  sum(taxa_sums(ps_claire))
  length(taxa_sums(ps_claire))
  
  ps_claire <- prune_taxa(taxa_sums(ps_claire) > 10, ps_claire)
  sum(taxa_sums(ps_claire))
  length(taxa_sums(ps_claire))
  
  #Filtering out ASVs that are present in less than a chosen fraction of samples (here 5%)
  ps_samuel <- prune_taxa(colSums(otu_table(ps_samuel) > 0) >= (0.05 * nsamples(ps_samuel)), ps_samuel)
  sum(taxa_sums(ps_samuel))
  length(taxa_sums(ps_samuel))
  
  ps_claire <- prune_taxa(colSums(otu_table(ps_claire) > 0) >= (0.05 * nsamples(ps_claire)), ps_claire)
  sum(taxa_sums(ps_claire))
  length(taxa_sums(ps_claire))
}
source(file = "~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/relab_analysis_graphs_and_stats.R")
pairs <- list(list("Wt:Vehicle","Wt:Putrescine"), list("IL-22ra1-/-:Vehicle","IL-22ra1-/-:Putrescine"), list("Wt:Vehicle","IL-22ra1-/-:Vehicle"), list("Wt:Putrescine","IL-22ra1-/-:Putrescine"))
customColors = list(list('black','#A22004'), list("#AB8F23","#04208D"), list("black","#AB8F23"),list("#A22004","#04208D"))
relabSpeciesPairwise(ps_samuel, deseq_samuel, measure = "log2fold", "gg_group", pairs, threshold = 0.01, customColors, "~/Documents/CHUM_git/figures/samuel/test_relab/")

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



#Differential abundance Samuel
deseq_samuel <- phyloseq_to_deseq2(ps_samuel, ~ genotype + treatment+ genotype:treatment) 

#Setting "Wt" as the baseline for genotype
colData(deseq_samuel)$genotype <- relevel(colData(deseq_samuel)$genotype, ref="Wt")

#Setting "Vehicle" as the baseline for treatment
colData(deseq_samuel)$treatment <- relevel(colData(deseq_samuel)$treatment, ref="Vehicle")

deseq_samuel <- DESeq(deseq_samuel, test="Wald", fitType = "parametric")

resultsNames(deseq_samuel)

customColors = list('black','#A22004',"#AB8F23","#04208D")
pairs <- list(list("Wt:Vehicle","Wt:Putrescine"), list("IL-22ra1-/-:Vehicle","IL-22ra1-/-:Putrescine"), list("Wt:Vehicle","IL-22ra1-/-:Vehicle"), list("Wt:Putrescine","IL-22ra1-/-:Putrescine"))

relabGroups(ps_samuel, deseq_samuel, measure = "log2fold", "gg_group", taxa = "Species", displayPvalue = FALSE, returnSigAsvs = FALSE, normalizeCounts = FALSE, threshold = 0.01, customColors, pairs, "~/Documents/CHUM_git/figures/filtering/samuel/relative_abundance_all_groups/")

#At other taxonomic levels
taxonomicLevels <- c("Genus","Family","Order","Class","Phylum")
# taxonomicLevels <- c("Phylum")
for(txnLevel in taxonomicLevels){
  
  #Creates ps subset for taxonomical level of interest
  ps_subset <- tax_glom(ps_samuel, taxrank = txnLevel)
  deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ genotype + treatment + genotype:treatment) 
  deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
  relabGroups(ps_subset, deseq_subset, measure = "log2fold", "gg_group", taxa = txnLevel, displayPvalue = FALSE, returnSigAsvs = FALSE, normalizeCounts = FALSE, threshold = 0.01, customColors, pairs, "~/Documents/CHUM_git/figures/filtering/samuel/relative_abundance_all_groups/")
}


#Testing new approach for Claire's timepoints#
#Load custom function to make the graphs for each species
source(file = "~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/relab_analysis_graphs_and_stats.R")

#Path where to save graphs
pathToSave <- "~/Documents/CHUM_git/figures/claire/relative_abundance_by_timepoint/"

#customColors for graph display
customColors = c("blue","red")

#Iterate through timepoints
for(timePoint in levels(sample_data(ps_claire)$week)){
  
  #New path created for each week
  newPath <- paste(pathToSave, "week_", timePoint, "/", sep = "")
  existingDirCheck(newPath)
  
  #Creating phyloseq objects for each timepoint
  ps_subset <- subset_samples(ps_claire, week == timePoint)
  
  #Simple deseq object only accounting for the differences in diet
  deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ diet) 
  
  #Performing the deseq analysis
  deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
  
  #For a given taxononical levels, creates graph for each timepoint, displaying which species were found to be differentially abundant
  relabSingleTimepoint(ps_subset, deseq_subset, measure = "log2fold", "diet", timePoint = timePoint, taxa = "Species", threshold = 0.01, customColors, newPath)  
  
}

#Additionnal approach to account for all timepoints, but same species, because we keep the species for which at least one point had significant differential abundances
source(file = "~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/relab_analysis_graphs_and_stats.R")
customColors = c("blue","red")

#At species level
relabTimelineRevised(ps_claire, measure = "log2fold", "week", "diet", "Species", threshold = 0.01, customColors,  "~/Documents/CHUM_git/figures/claire/test_timeline/")

#At other taxonomic levels
taxonomicLevels <- c("Genus","Family","Order","Class","Phylum")
for(txnLevel in taxonomicLevels){
  
  #Creates ps subset for taxonomical level of interest
  ps_subset <- tax_glom(ps_claire, taxrank = txnLevel)
  relabTimelineRevised(ps_subset, measure = "log2fold", "week", "diet", txnLevel, threshold = 0.01, customColors,  "~/Documents/CHUM_git/figures/claire/test_timeline/")
}

deseq_claire <- phyloseq_to_deseq2(ps_claire, ~ diet*week) 
deseq_claire <- DESeq(deseq_claire, test="Wald", fitType = "parametric")


###Correlations between relative abundances and other metrics
#Setting a correlation matrix workflow starting with Samuel's data
#Preparing deseq object needed for the function
deseq_samuel <- phyloseq_to_deseq2(ps_samuel, ~ genotype + treatment+ genotype:treatment) 

#Setting "Wt" as the baseline for genotype
colData(deseq_samuel)$genotype <- relevel(colData(deseq_samuel)$genotype, ref="Wt")

#Setting "Vehicle" as the baseline for treatment
colData(deseq_samuel)$treatment <- relevel(colData(deseq_samuel)$treatment, ref="Vehicle")

deseq_samuel <- DESeq(deseq_samuel, test="Wald", fitType = "parametric")

#Preparing dataframe for correlation
variables <- read.xlsx("~/Documents/CHUM_git/figures/samuel/correlation/Samuel's Combine Data_Correlation.xlsx")
rownames(variables) <- variables$X3
variables <- variables[,c(8:11)]
colnames(variables) <- c("lcn-2","il-6","tnf-alpha","colon_length")


customColors = list('black','#A22004',"#AB8F23","#04208D")
pairs <- list(list("Wt:Vehicle","Wt:Putrescine"), list("IL-22ra1-/-:Vehicle","IL-22ra1-/-:Putrescine"), list("Wt:Vehicle","IL-22ra1-/-:Vehicle"), list("Wt:Putrescine","IL-22ra1-/-:Putrescine"))
#For species level
#One heatmap per gg_group
correlationGroups(ps_samuel, deseq_samuel, measure = "log2fold", "gg_group", taxa = "Species", displayPvalue = FALSE, threshold = 0.01, customColors, pairs, "~/Documents/CHUM_git/figures/filtering/samuel/correlation/", df = variables, global = FALSE, showIndivCor = FALSE, normalizedCountsOnly = FALSE)
#One heatmap for all groups
correlationGroups(ps_samuel, deseq_samuel, measure = "log2fold", "gg_group", taxa = "Species", displayPvalue = FALSE, threshold = 0.01, customColors, pairs, "~/Documents/CHUM_git/figures/filtering/samuel/correlation/", df = variables, global = TRUE, showIndivCor = TRUE, normalizedCountsOnly = FALSE)

#For other taxonomical levels of interest
taxonomicLevels <- c("Genus","Family","Class","Order","Phylum")
for(txnLevel in taxonomicLevels){
  
  #Creates ps subset for taxonomical level of interest
  ps_subset <- tax_glom(ps_samuel, taxrank = txnLevel)
  
  #Preparing deseq object needed for the function
  deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ genotype + treatment+ genotype:treatment) 
  
  #Setting "Wt" as the baseline for genotype
  colData(deseq_subset)$genotype <- relevel(colData(deseq_subset)$genotype, ref="Wt")
  
  #Setting "Vehicle" as the baseline for treatment
  colData(deseq_subset)$treatment <- relevel(colData(deseq_subset)$treatment, ref="Vehicle")
  
  deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
  
  #One heatmap per gg_group
  correlationGroups(ps_subset, deseq_subset, measure = "log2fold", "gg_group", taxa = txnLevel, displayPvalue = FALSE, threshold = 0.01, customColors, pairs, "~/Documents/CHUM_git/figures/filtering/samuel/correlation/", df = variables, global = FALSE, showIndivCor = FALSE, normalizedCountsOnly = FALSE)
  #One heatmap for all groups
  correlationGroups(ps_subset, deseq_subset, measure = "log2fold", "gg_group", taxa = txnLevel, displayPvalue = FALSE, threshold = 0.01, customColors, pairs, "~/Documents/CHUM_git/figures/filtering/samuel/correlation/", df = variables, global = TRUE, showIndivCor = FALSE, normalizedCountsOnly = FALSE)
  
}


#Setting a correlation matrix workflow starting with Claire's data
#Preparing deseq object needed for the function
#Technically, only last timepoint should be useful
ps_claire_last_timepoint <- prune_samples(sample_data(ps_claire)$week == 14, ps_claire)
#replace tax names in the otu_table so that they are corresponding to the excel file
sample_names(ps_claire_last_timepoint) <- gsub("_.*", "", sample_names(ps_claire_last_timepoint))
sample_names(ps_claire_last_timepoint)

deseq_claire <- phyloseq_to_deseq2(ps_claire_last_timepoint, ~ diet) 
deseq_claire <- DESeq(deseq_claire, test="Wald", fitType = "parametric")



#Preparing dataframe for correlation
variables <- read.xlsx("~/Documents/CHUM_git/figures/claire/correlation/Claire_Data for heatmap_All_El_1 2024.xlsx")
rownames(variables) <- variables$ID
variables <- variables[,c(6:10)]


#For species level
#One heatmap for all groups
correlation2Var(ps_claire_last_timepoint, deseq_claire, measure = "log2fold", "gg_group", taxa = "Species", displayPvalue = FALSE, threshold = 0.01, "~/Documents/CHUM_git/figures/claire/correlation/week_14/", df = variables, global = TRUE, showIndivCor = FALSE, normalizedCountsOnly = FALSE)

#Week 8
ps_claire_last_timepoint <- prune_samples(sample_data(ps_claire)$week == 8, ps_claire)
#replace tax names in the otu_table so that they are corresponding to the excel file
sample_names(ps_claire_last_timepoint) <- gsub("_.*", "", sample_names(ps_claire_last_timepoint))
sample_names(ps_claire_last_timepoint)

deseq_claire <- phyloseq_to_deseq2(ps_claire_last_timepoint, ~ diet) 
deseq_claire <- DESeq(deseq_claire, test="Wald", fitType = "parametric")



#Preparing dataframe for correlation
variables <- read.xlsx("~/Documents/CHUM_git/figures/claire/correlation/Claire_Data for heatmap_All_El_1 2024.xlsx")
rownames(variables) <- variables$ID
variables <- variables[,c(6:10)]


#For species level
#One heatmap for all groups
correlation2Var(ps_claire_last_timepoint, deseq_claire, measure = "log2fold", "gg_group", taxa = "Species", displayPvalue = FALSE, threshold = 0.01, "~/Documents/CHUM_git/figures/claire/correlation/week_8/", df = variables, global = TRUE, showIndivCor = FALSE, normalizedCountsOnly = FALSE)

#Week 10
ps_claire_last_timepoint <- prune_samples(sample_data(ps_claire)$week == 10, ps_claire)
#replace tax names in the otu_table so that they are corresponding to the excel file
sample_names(ps_claire_last_timepoint) <- gsub("_.*", "", sample_names(ps_claire_last_timepoint))
sample_names(ps_claire_last_timepoint)

deseq_claire <- phyloseq_to_deseq2(ps_claire_last_timepoint, ~ diet) 
deseq_claire <- DESeq(deseq_claire, test="Wald", fitType = "parametric")



#Preparing dataframe for correlation
variables <- read.xlsx("~/Documents/CHUM_git/figures/claire/correlation/Claire_Data for heatmap_All_El_1 2024.xlsx")
rownames(variables) <- variables$ID
variables <- variables[,c(6:10)]


#For species level
#One heatmap for all groups
correlation2Var(ps_claire_last_timepoint, deseq_claire, measure = "log2fold", "gg_group", taxa = "Species", displayPvalue = FALSE, threshold = 0.01, "~/Documents/CHUM_git/figures/claire/correlation/week_10/", df = variables, global = TRUE, showIndivCor = FALSE, normalizedCountsOnly = FALSE)



#For other taxonomical levels of interest
taxonomicLevels <- c("Genus","Family","Class","Order","Phylum")
for(txnLevel in taxonomicLevels){
  
  #Creates ps subset for taxonomical level of interest
  ps_subset <- tax_glom(ps_samuel, taxrank = txnLevel)
  
  #Preparing deseq object needed for the function
  deseq_subset <- phyloseq_to_deseq2(ps_subset, ~ genotype + treatment+ genotype:treatment) 
  
  #Setting "Wt" as the baseline for genotype
  colData(deseq_subset)$genotype <- relevel(colData(deseq_subset)$genotype, ref="Wt")
  
  #Setting "Vehicle" as the baseline for treatment
  colData(deseq_subset)$treatment <- relevel(colData(deseq_subset)$treatment, ref="Vehicle")
  
  deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
  
  #One heatmap per gg_group
  correlationGroups(ps_subset, deseq_subset, measure = "log2fold", "gg_group", taxa = txnLevel, displayPvalue = FALSE, threshold = 0.01, customColors, pairs, "~/Documents/CHUM_git/figures/filtering/samuel/correlation/", df = variables, global = FALSE, showIndivCor = FALSE, normalizedCountsOnly = FALSE)
  #One heatmap for all groups
  correlationGroups(ps_subset, deseq_subset, measure = "log2fold", "gg_group", taxa = txnLevel, displayPvalue = FALSE, threshold = 0.01, customColors, pairs, "~/Documents/CHUM_git/figures/filtering/samuel/correlation/", df = variables, global = TRUE, showIndivCor = FALSE, normalizedCountsOnly = FALSE)
  
}

