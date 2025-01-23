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
  library(StackbarExtended) # Thibault C. package
  
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
# and for picrust2
{
#function filtering out ASVs for which they were in total less than a threshold count
ps_samuel <- prune_taxa(taxa_sums(ps_samuel) > 0, ps_samuel) #Actual ASVs from Samuel's data
sum(taxa_sums(ps_samuel))
length(taxa_sums(ps_samuel))

ps_samuel <- prune_taxa(taxa_sums(ps_samuel) > 10, ps_samuel)
sum(taxa_sums(ps_samuel))
length(taxa_sums(ps_samuel))

ps_claire <- prune_taxa(taxa_sums(ps_claire) > 0, ps_claire) #Actual ASVs from Claire's data
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

# Function to produce picrust2 required inputs
producePicrust2Inputs <- function(ps, output_dir){
  
  # Creates picrust2 output folder if does not exist yet
  existingDirCheck(paste0(output_dir, "/picrust2/"))
  existingDirCheck(paste0(output_dir, "/picrust2/input"))
  
  # Define the output file path
  output_file <- file.path(paste0(output_dir, "/picrust2/input"), "seqs.fna")
  
  # Open a connection to write the file
  file_conn <- file(output_file, open = "w")
  
  for (i in seq_along(refseq(ps))){
    
    seq = as.character(refseq(ps)[i]) # Sequence
    id = as.character(names(refseq(ps)[i])) # ID (ASV number associated with sequence)
    
    # Write the FASTA entry to the file
    writeLines(paste0(">", id), file_conn) # FASTA header
    writeLines(seq, file_conn)            # Sequence
    
  }
  
  # Close the file connection
  close(file_conn)
  message("Produced fasta file successfully.")
  
  # Produce the biom file
  # Get otu_table
  otu_table <- as.data.frame(t(otu_table(ps)))
  
  # Convert table to desired format
  otu_table <- tibble::rownames_to_column(otu_table, var = "#OTU ID")
  
  # Add the header
  header <- "# Constructed from biom file"
  
  # Write the OTU table to a text file
  output_file <- file.path(paste0(output_dir, "/picrust2/input"), "otu_table.txt")
  write(header, file = output_file)
  write.table(otu_table, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)
  message("Produced otu_table file successfully.")
}

producePicrust2Inputs(ps_samuel, "~/Documents/CHUM_git/Microbiota_17")

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










# Stackbar extended at family level for Claire and Samuel's data 
# might be better to do separate graphs if you want to see the stats separately
# because logically they can only be displayed if you do one comparaison but not 
# accross multiple groups
# adding a table under the thing might be easier for this for sure
# Samuel's data
#put cage as factor
sample_data(ps_samuel)$cage = as.factor(sample_data(ps_samuel)$cage)

#for Samuel's data, put gg_group as factor and define order
sample_data(ps_samuel)$gg_group <- factor(sample_data(ps_samuel)$gg_group, levels = c("Wt:Vehicle", "Wt:Putrescine", "IL-22ra1-/-:Vehicle", "IL-22ra1-/-:Putrescine"))  # Vehicle as reference
sample_data(ps_samuel)$genotype <- factor(sample_data(ps_samuel)$genotype, levels = c("Wt", "IL-22ra1-/-"))
sample_data(ps_samuel)$treatment <- factor(sample_data(ps_samuel)$treatment, levels = c("Vehicle", "Putrescine"))

il22_exp_family <- plot_microbiota(
  ps_object = ps_samuel,
  exp_group = 'gg_group',
  sample_name = 'sample_id',
  hues = c("Purples", "Blues", "Greens", "Oranges", "Reds"),
  differential_analysis = T,
  sig_lab = T,
  n_row = 2,
  n_col = 2,
  fdr_threshold = 0.05,
  main_level = "Family",
  n_phy = 5, # number of taxa to show 
  mult_comp = T, # pairwise comparaisons for diff ab analysis
  selected_comparisons = list(c("Wt:Vehicle", "Wt:Putrescine"), c("IL-22ra1-/-:Vehicle", "IL-22ra1-/-:Putrescine"), c("Wt:Vehicle","IL-22ra1-/-:Vehicle"), c("Wt:Putrescine","IL-22ra1-/-:Putrescine"))
  )

print(il22_exp_family$plot)
print(il22_exp_family$significant_table_main)
#Nothing for Wt:Putrescine_vs_IL-22ra1-/-:Putrescine at the Family level.
# When there was nothing significant, there is not entry for the significance
# table for the comparaison of interest 

# Extract plot and customize it
plot <- il22_exp_family$plot
p <- plot + theme(
  text = element_text(family = "Times New Roman"),      # Global text settings
  strip.text = element_text(size = 14, face = "bold"),  # Facet titles
  plot.title = element_text(size = 20, face = "bold"),  # Main title
  axis.title = element_text(size = 15, face = "bold"),  # Axis titles
  axis.text = element_text(size = 12, face = "bold"),   # Axis text
  legend.title = element_text(face = "bold", size = 14)  # Legend title  # Legend text
) +
  guides(fill = "none")+
  labs(x = "Sample ID")
p

# Save plot and associated stats
existingDirCheck("../figures/samuel/stackbar")
ggsave(plot = p, filename = "../figures/samuel/stackbar/family_stackbar.png", width = 8, height = 8, dpi = 300)

# Function to write and save stackbarExtended sig_table
writeStackbarExtendedSigTable <- function(SckbarExtObject,filepath){
  
  # Initialize empty dataframe to append tables to
  table_to_write <- data.frame()
  
  # Iterate over the list of tables
  for (i in seq_along(SckbarExtObject$significant_table_main)) {
    
    # Extract the table
    table <- SckbarExtObject$significant_table_main[[i]]
    
    # Add a column with the name of the current table
    table$comparaison <- names(SckbarExtObject$significant_table_main)[i]
    
    # Append to the master table
    table_to_write <- rbind(table_to_write, table)
  }
  
  write_xlsx(x = table_to_write, path = filepath)

}
writeStackbarExtendedSigTable(il22_exp_family, filepath = "../figures/samuel/stackbar/family_stackbar.xlsx")


# Claire's data
#for Claire's data, put gg_group as factor and define order
sample_data(ps_claire)$gg_group <- factor(sample_data(ps_claire)$gg_group, levels = c("3:50", "8:50", "10:50", "14:50", "3:500", "8:500", "10:500", "14:500"))  
# Enables to rename factor levels (this is more convenient to change the names for each subgraph
# in the facet_wrap.)
library(forcats)
sample_data(ps_claire)$gg_group <- fct_recode(sample_data(ps_claire)$gg_group,
                          "50 ppm / 3w"="3:50",
                          "50 ppm / 8w"="8:50",
                          "50 ppm / 10w"="10:50",
                          "50 ppm / 14w"="14:50",
                          "500 ppm / 3w"="3:500",
                          "500 ppm / 8w"="8:500",
                          "500 ppm / 10w"="10:500",
                          "500 ppm / 14w"="14:500")

iron_exp_family <- plot_microbiota(
  ps_object = ps_claire,
  exp_group = 'gg_group',
  sample_name = 'sample_id',
  hues = c("Purples", "Blues", "Greens", "Oranges", "Reds"),
  differential_analysis = T,
  sig_lab = T,
  n_row = 2,
  n_col = 4,
  fdr_threshold = 0.05,
  main_level = "Family",
  n_phy = 5, # number of taxa to show 
  mult_comp = T, # pairwise comparaisons for diff ab analysis
  selected_comparisons = list(c( "50 ppm / 3w",  "500 ppm / 3w"), c( "50 ppm / 8w",  "500 ppm / 8w"), c( "50 ppm / 10w",  "500 ppm / 10w"), c("50 ppm / 14w", "500 ppm / 14w"))
)

print(iron_exp_family$plot)
print(iron_exp_family$significant_table_main)
plot <- iron_exp_family$plot

# Custom the plot
p <- plot + theme(
  text = element_text(family = "Arial"),      # Global text settings
  strip.text = element_text(size = 14, face = "bold"),  # Facet titles
  plot.title = element_text(size = 20, face = "bold"),  # Main title
  axis.title = element_text(size = 15, face = "bold"),  # Axis titles
  axis.text = element_text(size = 12, face = "bold"),   # Axis text
  legend.title = element_text(face = "bold", size = 14)  # Legend title  # Legend text
) +
  scale_x_discrete(labels = function(x) substr(x, 1, 5))+
  guides(fill = "none")+
  labs(x = "Sample ID")
p

# Saving the plot and the associated stats
existingDirCheck("../figures/claire/stackbar")
ggsave(plot = p, filename = "../figures/claire/stackbar/family_stackbar.png", width = 14, height = 8, dpi = 300)
writeStackbarExtendedSigTable(iron_exp_family, filepath = "../figures/claire/stackbar/family_stackbar.xlsx")











# Saving the phyla data with the stats
#phyla distribution using relative abundance
# Samuel's data
#for Samuel's data, put gg_group as factor and define order
sample_data(ps_samuel)$gg_group <- factor(sample_data(ps_samuel)$gg_group, levels = c("Wt:Vehicle", "Wt:Putrescine", "IL-22ra1-/-:Vehicle", "IL-22ra1-/-:Putrescine"))  # Vehicle as reference
sample_data(ps_samuel)$genotype <- factor(sample_data(ps_samuel)$genotype, levels = c("Wt", "IL-22ra1-/-"))
sample_data(ps_samuel)$treatment <- factor(sample_data(ps_samuel)$treatment, levels = c("Vehicle", "Putrescine"))

il22_exp_family <- plot_microbiota(
  ps_object = ps_samuel,
  exp_group = 'gg_group',
  sample_name = 'sample_id',
  hues = c("Purples", "Blues", "Greens", "Oranges", "Reds"),
  differential_analysis = T,
  sig_lab = T,
  n_row = 2,
  n_col = 2,
  fdr_threshold = 0.05,
  main_level = "Phylum",
  n_phy = 5, # number of taxa to show 
  mult_comp = T, # pairwise comparaisons for diff ab analysis
  selected_comparisons = list(c("Wt:Vehicle", "Wt:Putrescine"), c("IL-22ra1-/-:Vehicle", "IL-22ra1-/-:Putrescine"), c("Wt:Vehicle","IL-22ra1-/-:Vehicle"), c("Wt:Putrescine","IL-22ra1-/-:Putrescine"))
)

print(il22_exp_family$plot)
print(il22_exp_family$significant_table_main)
#Nothing for Wt:Putrescine_vs_IL-22ra1-/-:Putrescine at the Family level.
# When there was nothing significant, there is not entry for the significance
# table for the comparaison of interest 

# Extract plot and customize it
plot <- il22_exp_family$plot
p <- plot + theme(
  text = element_text(family = "Times New Roman"),      # Global text settings
  strip.text = element_text(size = 14, face = "bold"),  # Facet titles
  plot.title = element_text(size = 20, face = "bold"),  # Main title
  axis.title = element_text(size = 15, face = "bold"),  # Axis titles
  axis.text = element_text(size = 12, face = "bold"),   # Axis text
  legend.title = element_text(face = "bold", size = 14)  # Legend title  # Legend text
) +
  #guides(fill = "none")+
  labs(x = "Sample ID")
p

# Save plot and associated stats
existingDirCheck("../figures/samuel/stackbar")
ggsave(plot = p, filename = "../figures/samuel/stackbar/family_stackbar.png", width = 8, height = 8, dpi = 300)

# Function to write and save stackbarExtended sig_table
writeStackbarExtendedSigTable <- function(SckbarExtObject,filepath){
  
  # Initialize empty dataframe to append tables to
  table_to_write <- data.frame()
  
  # Iterate over the list of tables
  for (i in seq_along(SckbarExtObject$significant_table_main)) {
    
    # Extract the table
    table <- SckbarExtObject$significant_table_main[[i]]
    
    # Add a column with the name of the current table
    table$comparaison <- names(SckbarExtObject$significant_table_main)[i]
    
    # Append to the master table
    table_to_write <- rbind(table_to_write, table)
  }
  
  write_xlsx(x = table_to_write, path = filepath)
  
}
writeStackbarExtendedSigTable(il22_exp_family, filepath = "../figures/samuel/stackbar/family_stackbar.xlsx")














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

