library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)

source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/utilities.R")

# Analysis of picrust2 results for last timepoint
setwd("~/Documents/CHUM_git/Microbiota_18/picrust2/picrust2_out_pipeline/")

ko_metagenome <- list(read.table("KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE),
                      read.table("KO_metagenome_out/seqtab_norm.tsv.gz", sep = "\t", header = TRUE),
                      read.table("KO_metagenome_out/weighted_nsti.tsv.gz", sep = "\t", header = TRUE))

View(ko_metagenome[[1]]) # ko numbers
# Remove the X at the start of the names
colnames(ko_metagenome[[1]])[2:ncol(ko_metagenome[[1]])] = substring(colnames(ko_metagenome[[1]][2:ncol(ko_metagenome[[1]])]), first = 2, 6)
colnames(ko_metagenome[[1]])[1] <- "ko_id"

View(ko_metagenome[[2]]) # osef
View(ko_metagenome[[3]]) # osef


ec_metagenome <- list(read.table("EC_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE, row.names = 1),
                      read.table("EC_metagenome_out/seqtab_norm.tsv.gz", sep = "\t", header = TRUE),
                      read.table("EC_metagenome_out/weighted_nsti.tsv.gz", sep = "\t", header = TRUE))

View(ec_metagenome[[1]]) # ec numbers
View(ec_metagenome[[2]]) # osef
View(ec_metagenome[[3]]) # osef


marker <- read.table("marker_predicted_and_nsti.tsv.gz", sep = "\t", header = TRUE)
View(marker)

ko_predicted <- read.table("KO_predicted.tsv.gz", sep = "\t", header = TRUE, row.names = 1)


# colnames(ko_metagenome[[1]])[1] = "KO_id"

# Metadata handling
# Metadata handling
{
  #loading metadata of interest
  metadata <- read.csv("metadata/metadata.csv", sep = ";")
  
  # Remove the non-metadata stuff (liver measures and stuff)
  metadata <- metadata[,-c(5:8)]
  
  # Remove the letter at the end of id
  metadata$id <- substring(metadata$id, 1, 5)
  
  #adding id col as rownames too
  rownames(metadata) <- metadata$id
  
  # Remove dead mouse
  metadata <- metadata[-46,]
  
  # Extract 16S reads sample ids
  samples <- read.xlsx("metadata/Microbiota_18_samples_2025-01-13.xlsx")
  samples <- as.data.frame(samples$Nom)
  colnames(samples) <- "sample_id"
  samples$id <- substring(samples$sample_id, 1, 5)
  samples$timepoint <- substring(samples$sample_id, 8, nchar(samples$sample_id)) 
  
  # Bind both metadata df to link timepoints with their metadata (diet and treatment)
  metadata <- merge(samples, metadata, by = "id")
  
  # Consider timepoint 53 similar as timepoint 54
  metadata[metadata$timepoint=="53","timepoint"] <- "54"
  
  # Add week column 
  metadata$week <- ifelse(metadata$timepoint == "final", "18", as.character(round(as.numeric(metadata$timepoint)/7, 1)+3))
  
  # Adding gg_group variable (combination of time, diet and treatment)
  metadata$gg_group <- 
    paste(metadata$timepoint, 
          metadata$diet,
          metadata$treatment, 
          sep = ":")
  
  # Another gg_group variable (diet and treatment only)
  metadata$gg_group2 <- 
    paste(metadata$diet,
          metadata$treatment, 
          sep = ":")
  
  # Put full_id as rownames
  rownames(metadata) <- metadata$sample_id
}

# Only Samuel's metadata
meta <- metadata[metadata$student == "Samuel",]
meta <- meta[,-c(2,4)]
rownames(meta) <- 1:nrow(meta)
View(meta)

# Testing ggpicrust2 workflow
# Run ggpicrust2 with input data
{
  results_data_input <- ggpicrust2(#file <-"KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
    data = ko_metagenome[[1]],
    metadata = meta,
    group = "gg_group", # For example dataset, group = "Environment"
    pathway = "KO",
    daa_method = "ALDEx2",
    ko_to_kegg = TRUE,
    order = "pathway_class",
    p_values_bar = TRUE,
    x_lab = "pathway_name")
  
  results_data_input[[1]]$plot
  results_data_input[[2]]$plot
  results_data_input[[1]]$results
  
  data("metadata")
  data("ko_abundance")
} # This is pretty much awful and does not make any sense at all, I guess we are doing it on our own

# methods <- c("ALDEx2", "DESeq2", "edgeR")
daa_pathways_heatmap_pairwise <- function(pathways, map_df, custom_colors, metadata, group, threshold = 0.05, daa_method = "DESeq2"){
  
  # Map pathway names
  pathways <- merge(pathways, map_df, by = "pathway")
  pathways <- cbind(pathways[,ncol(pathways)], pathways[,-c(1,ncol(pathways))]) # Put full pathways names as first column
  colnames(pathways)[1] <- "pathway"
  
  # Differently formatted with pathway full names as row.names
  path_ab <- pathways
  row.names(path_ab) <- pathways[,1]
  path_ab <- path_ab[,-1]
  
  # First, identify differentially abundant pathways
  daa_results <- pathway_daa(abundance = path_ab, metadata = metadata, group = group, daa_method = daa_method, select = NULL, p.adjust = "BH", reference = "NULL")
  daa_results$p_adjust[is.na(daa_results$p_adjust)] <- 1 # Replace NAs values by 1
  daa_res_sub <- daa_results[daa_results$p_adjust < threshold,] # Keep only p-values < 0.05
  path_sig <- pathways[pathways$pathway %in% daa_res_sub$feature,] # Keep only features that are differentially abundant
  
  # Enrich daa_res_sub with additionnal columns for p_value heatmaps
  daa_res_sub$comparison <- paste(daa_res_sub$group1, "VS", daa_res_sub$group2) 
  daa_res_sub$significance <- cut(daa_res_sub$p_adjust,
                                  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                  labels = c("<0.001", "<0.01", "<0.05", "n.s."),
                                  right = FALSE)
  
  # Plot heatmap
  row.names(path_sig) <- path_sig[,1]
  path_sig <- path_sig[,-1]
  ptwy_hmap <- pathway_heatmap(path_sig, metadata, group, colors = custom_colors)+
    theme(
      panel.grid.major = element_line(color = "black", linewidth  = 20),
      axis.text.y = element_text(face = "italic"),
      strip.text = element_text(color = "white"))
  
  # Build a heatmap for p-values
  p_val_hmap <- ggplot(daa_res_sub, aes(y = feature, x = comparison, fill = significance)) +
    geom_tile(color = "black", lwd = 0.5, linetype = 1) +
    coord_fixed() + # Makes thing squared
    scale_fill_manual(
      values = c("n.s." = "white", "<0.05" = "#F4A3A8", "<0.01" = "#E04B54", "<0.001" = "#A40000"))+
    theme_minimal() +
    labs(x = "", y = "", fill = "P value") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  return(list(ptwy_hmap, p_val_hmap))
  
}

# Function to prepare the sub dataframes according to the given pair of groups
# Returns a list of dataframes, first df is pathways df for stats, second is same but for graph display and third is sub metadata
df_preparation <- function(groups){
  
  # First sub pathway df for statistical analysis only
  # pathways1 <- read.table("pathways_out/path_abun_unstrat.tsv.gz", sep = "\t", header = TRUE, row.names = 1)
  # colnames(pathways1) = substring(colnames(pathways1), first = 2, 6) # Remove the X at the start of the names
  # path1_sub <- pathways1[, colnames(pathways1) %in% meta$sample_id[meta$gg_group %in% groups]]
  
  # Second sub pathway df for graph
  pathways <- read.table("pathways_out/path_abun_unstrat.tsv.gz", sep = "\t", header = TRUE)
  colnames(pathways)[2:ncol(pathways)] = substring(colnames(pathways[2:ncol(pathways)]), first = 2, 6) # Remove the X at the start of the names
  path_sub <- pathways[, c(TRUE, colnames(pathways)[2:ncol(pathways)] %in% meta$sample_id[meta$gg_group %in% groups])] # Keep first col
  
  # Sub for metadata
  meta_sub <- meta[meta$gg_group %in% groups,]
  
  return(list(path_sub, meta_sub))
  
}

# Mapping for metacyc id and pathway names
path_names <- readLines("../../../../picrust2 database/pathway_names.txt")
path_ids <- readLines("../../../../picrust2 database/metacyc_ids.txt")
map_df <- data.frame(pathway_name = path_names, pathway = path_ids)

meta$gg_group <- factor(meta$gg_group, levels = c("Wt:Vehicle", "Wt:Putrescine", "IL-22ra1-/-:Vehicle", "IL-22ra1-/-:Putrescine"))
customColors = c('black','#A22004',"#AB8F23","#04208D")

dfs <- df_preparation(groups = c("Wt:Vehicle", "Wt:Putrescine"))
customColors = c('black','#A22004')
plots <- daa_pathways_heatmap_pairwise(pathways = dfs[[1]], map_df = map_df, custom_colors = customColors, metadata = dfs[[2]], group = "gg_group", threshold = 0.05)
plots[[1]]
plots[[2]]


dfs <- df_preparation(groups = c("IL-22ra1-/-:Vehicle", "IL-22ra1-/-:Putrescine"))
customColors = c("#AB8F23","#04208D")
plots <- daa_pathways_heatmap_pairwise(pathways =  dfs[[1]], map_df = map_df, custom_colors = customColors, metadata =  dfs[[2]], group = "gg_group", threshold = 0.05)
plots[[1]]
plots[[2]]

# # Load file with KOs annotations for compound/pathways of interest
# ko_annotations <- readxl::read_excel("~/Documents/CHUM_git/picrust2 database/compound_to_kos.xlsx")
# 
# # Function that iterates through list of compounds/pathways, subtract KO_abundance dataframe for KOs associated, and returns multiple graphs
# # ko_abundance needs to have kos as rownames
# KOsToCompoundAbundanceGraphs <- function(ko_abundance, ko_annotations, metadata, group, group_order, sample_id_col, path){
#   
#   existingDirCheck(path)
#   
#   # Iterate through the compounds/pathways
#   for(i in 1:nrow(ko_annotations)){
#     compound <- ko_annotations$`compound/pathway`[i]
#     kos <- unlist(strsplit(ko_annotations$Kos[i], ";"))
#     abundance <- ko_abundance[rownames(ko_abundance) %in% kos, ]
#     total_ab <- as.data.frame(colSums(abundance))
#     
#     # Associate gg_group with each sample_id using metadata information
#     for(i in 1:nrow(total_ab)){
#     total_ab[[group]][i] <- metadata[metadata[[sample_id_col]]==rownames(total_ab)[i],group]
#     }
#     
#     # Replacing a colname and adding a sample_id column, + gg_plot as a factor
#     colnames(total_ab)[1] <- compound
#     total_ab$sample_id <- rownames(total_ab)
#     total_ab[[group]] <- factor(total_ab[[group]], levels = group_order)
#     
#     p <- ggplot(total_ab, aes(x = .data[[group]], y = .data[[compound]], colors = .data[[group]])) +
#       geom_bar(stat = "identity") +
#       theme_minimal() +
#       labs(title = paste("Predicted", clean_string(compound), "levels", sep = " "), y = paste("Normalized", clean_string(compound), sep = " "))
#     p
#     
#     # Save plot
#     ggsave(plot = p, filename = paste0(path, "/", clean_string(compound), "_predicted.png"), dpi = 300, width = 5, height = 5, bg = "white")
#   }
# }
# 
# KOsToCompoundAbundanceGraphs(ko_abundance = ko_metagenome[[1]],
#                             ko_annotations = ko_annotations,
#                             metadata = meta,
#                             group = "gg_group",
#                             group_order = c("Wt:Vehicle", "Wt:Putrescine", "IL-22ra1-/-:Vehicle", "IL-22ra1-/-:Putrescine"),
#                             sample_id_col="sample_id",
#                             path = "~/Documents/CHUM_git/figures/samuel/picrust2")




