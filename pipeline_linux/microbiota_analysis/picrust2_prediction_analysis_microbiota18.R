library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(readxl)
library(openxlsx)

source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/utilities.R")

# Analysis of picrust2 results for last timepoint
setwd("~/Documents/CHUM_git/Microbiota_18/picrust2/picrust2_out_pipeline/")

ko_metagenome <- list(read.table("KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE),
                      read.table("KO_metagenome_out/seqtab_norm.tsv.gz", sep = "\t", header = TRUE),
                      read.table("KO_metagenome_out/weighted_nsti.tsv.gz", sep = "\t", header = TRUE))

View(ko_metagenome[[1]]) # ko numbers
# Remove the X at the start of the names
colnames(ko_metagenome[[1]])[2:ncol(ko_metagenome[[1]])] = substring(colnames(ko_metagenome[[1]][2:ncol(ko_metagenome[[1]])]), first = 2, 6)
rownames(ko_metagenome[[1]]) <- ko_metagenome[[1]]$function.
ko_metagenome[[1]] <- ko_metagenome[[1]][,-1]

View(ko_metagenome[[2]]) # osef
View(ko_metagenome[[3]]) # osef


ec_metagenome <- list(read.table("EC_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE, row.names = 1),
                      read.table("EC_metagenome_out/seqtab_norm.tsv.gz", sep = "\t", header = TRUE),
                      read.table("EC_metagenome_out/weighted_nsti.tsv.gz", sep = "\t", header = TRUE))

View(ec_metagenome[[1]]) # ec numbers
# Remove the X at the start of the names
colnames(ec_metagenome[[1]])[2:ncol(ec_metagenome[[1]])] = substring(colnames(ec_metagenome[[1]][2:ncol(ec_metagenome[[1]])]), first = 2, 6)
rownames(ec_metagenome[[1]]) <- ec_metagenome[[1]]$function.
ec_metagenome[[1]] <- ec_metagenome[[1]][,-1]
View(ec_metagenome[[2]]) # osef
View(ec_metagenome[[3]]) # osef


marker <- read.table("marker_predicted_and_nsti.tsv.gz", sep = "\t", header = TRUE)
View(marker)

ko_predicted <- read.table("KO_predicted.tsv.gz", sep = "\t", header = TRUE, row.names = 1)




# colnames(ko_metagenome[[1]])[1] = "KO_id"

# Metadata handling
{
  #loading metadata of interest
  metadata <- read.csv("../../metadata/metadata.csv", sep = ";")
  
  # Remove the non-metadata stuff (liver measures and stuff)
  metadata <- metadata[,-c(5:8)]
  
  # Remove the letter at the end of id
  metadata$id <- substring(metadata$id, 1, 5)
  
  #adding id col as rownames too
  rownames(metadata) <- metadata$id
  
  # Remove dead mouse
  metadata <- metadata[-46,]
  
  # Extract 16S reads sample ids
  samples <- read.xlsx("../../metadata/Microbiota_18_samples_2025-01-13.xlsx")
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

# Only final timepoint metadata
meta <- metadata[metadata$timepoint == "final",]
meta <- meta[,-c(2,4)]
rownames(meta) <- 1:nrow(meta)
rownames(meta) <- meta$id
meta$diet <- factor(meta$diet, levels = c("50", "500"))
meta$treatment <- factor(meta$treatment, levels = c("water", "dss"))
meta$gg_group2 <- factor(meta$gg_group2, levels = c("50:water", "500:water", "50:dss", "500:dss"))
View(meta)

# Rearranges ko_ab df cols order so that it corresponds to metadata
# target_order <- rownames(meta)
# reorder_indices <- match(target_order, colnames(ko_metagenome[[1]]))
# ko_metagenome[[1]] <- ko_metagenome[[1]][, reorder_indices]



# # DESEq2 analysis on ko abundance table output
# library(DESeq2)
# dds <- DESeqDataSetFromMatrix(countData = ko_metagenome[[1]],
#                               colData = meta,
#                               design = ~ treatment*diet)
# 
# library(ALDEx2)
# # CLR transformation of data
# conds <- c(meta$gg_group2)
# data.clr <- aldex.clr(ko_metagenome[[1]], conds, mc.samples = 128, denom = "all")

# Load file with KOs annotations for compound/pathways of interest
ko_annotations <- readxl::read_excel("~/Documents/CHUM_git/picrust2 database/compound_to_kos.xlsx")
ko_annotations <- readxl::read_excel("~/Documents/CHUM_git/picrust2 database/compound_to_kos2.xlsx")


# Focus on proprionate, butyrate and acetate 
ko_annotations <- ko_annotations[1:3,]

# Function that iterates through list of compounds/pathways, subtract KO_abundance dataframe for KOs associated, and returns multiple graphs
# ko_abundance needs to have kos as rownames
KOsToCompoundAbundanceGraphs <- function(ko_abundance, ko_annotations, metadata, group, customColors, sample_id_col, path){

  existingDirCheck(path)

  # Iterate through the compounds/pathways
  for(i in 1:nrow(ko_annotations)){
    compound <- ko_annotations$`compound/pathway`[i]
    kos <- unlist(strsplit(ko_annotations$Kos[i], ";"))
    abundance <- ko_abundance[rownames(ko_abundance) %in% kos, ]
    total_ab <- as.data.frame(colSums(abundance))

    # # Associate gg_group with each sample_id using metadata information
    # for(i in 1:nrow(total_ab)){
    # total_ab[[group]][i] <- metadata[metadata[[sample_id_col]]==rownames(total_ab)[i],group]
    # }
    
    total_ab <- merge(total_ab, meta, by = "row.names")

    # Replacing a colname and adding a sample_id column
    colnames(total_ab)[2] <- compound

    p <- ggplot(total_ab, aes(x = .data[[group]], y = .data[[compound]], color = .data[[group]])) +
      geom_point(size = 1, 
                 position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) + 
      
      #Error bars
      stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
                   aes(color = .data[[group]]),
                   width = 0.2, size = 0.7,
                   position = position_dodge(-0.75)) +
      
      #Mean lines
      stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
                   aes(ymin = ..y.., ymax = ..y.., group = .data[[group]]),
                   color = "black", linewidth = 0.5, width = 0.5,
                   position = position_dodge(-0.75))+
      
      #Connecting mean points with lines
      # stat_summary(fun = mean, geom = "line", size = 1.2) +  # Connecting means with lines
      
      
      
      labs(title = compound, y = paste("Normalized Abundance"), color = "Group", x = "") +
      scale_color_manual(values = customColors)+
      scale_y_continuous(limits = c(0, NA))+
      
      theme(
        plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
        axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
        axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust x-axis tick label font size
        axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
        legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
        legend.text = element_text(size = 12),  # Adjust legend font size
        panel.grid.major = element_blank(),  # Add major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black", size = 1),
        panel.background = element_blank()) # Include axis lines  # Include axis bar
    # Save plot
    ggsave(plot = p, filename = paste0(path, "/", clean_string(compound), "_predicted.png"), dpi = 300, width = 5, height = 5, bg = "white")
  }
}

KOsToCompoundAbundanceGraphs(ko_abundance = ko_metagenome[[1]],
                            ko_annotations = ko_annotations,
                            metadata = meta,
                            group = "gg_group2",
                            customColors = c("blue","red","darkblue","darkred"),
                            sample_id_col="id",
                            path = "~/Documents/CHUM_git/figures/thibault_new/picrust2")


# Stats with mixed effect models
# Iterate through the compounds/pathways
for(i in 1:nrow(ko_annotations)){
  compound <- ko_annotations$`compound/pathway`[i]
  kos <- unlist(strsplit(ko_annotations$Kos[i], ";"))
  abundance <- ko_metagenome[[1]][rownames(ko_metagenome[[1]]) %in% kos, ]
  total_ab <- as.data.frame(colSums(abundance))
  
  # # Associate gg_group with each sample_id using metadata information
  # for(i in 1:nrow(total_ab)){
  # total_ab[[group]][i] <- metadata[metadata[[sample_id_col]]==rownames(total_ab)[i],group]
  # }
  
  total_ab <- merge(total_ab, meta, by = "row.names")
  
  # Replacing a colname and adding a sample_id column
  colnames(total_ab)[2] <- compound
  
}


library(lme4)
library(car)
# Assuming your data frame 'df' includes:
# - scfa: predicted normalized SCFA abundance
# - treatment: factor with levels (e.g., "water" and "DSS")
# - diet: factor with levels (e.g., "50ppm" and "500ppm")
# - cage: a factor indicating cage identity
model <- glmer(Propionate ~ treatment * diet, data = total_ab,
              family = Gamma(link = "log"))

total_ab$log_scfa <- log(total_ab$Propionate)
model <- lm(Propionate ~ treatment * diet, data = total_ab)
model <- glm(Propionate ~ treatment * diet, data = total_ab, family = poisson)
model <- glm(Propionate ~ treatment * diet, data = total_ab, family = Gamma)
summary(model)

# Plot residuals vs fitted plot = check if non-linear patterm
plot(model, which = 1)

qqnorm(resid(model))
qqline(resid(model))

hist(total_ab$Propionate[total_ab$gg_group2 == "50:dss"], breaks = 20, main = "Histogram of Residuals")

summary(model)
Anova(model, type = "III")

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
daa_pathways_heatmap_pairwise <- function(pathways, map_df, custom_colors, metadata, group, threshold = 0.05, daa_method = "DESeq2", pAdjust = TRUE){
  
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
  if(pAdjust){
    daa_results$p_adjust[is.na(daa_results$p_adjust)] <- 1 # Replace NAs values by 1
    daa_res_sub <- daa_results[daa_results$p_adjust < threshold,]} # Keep only p-values < 0.05
  else{
    daa_results$p_value[is.na(daa_results$p_value)] <- 1 # Replace NAs values by 1
    daa_res_sub <- daa_results[daa_results$p_value < threshold,] # Keep only p-values < 0.05
  }
  path_sig <- pathways[pathways$pathway %in% daa_res_sub$feature,] # Keep only features that are differentially abundant
  
  # Enrich daa_res_sub with additionnal columns for p_value heatmaps
  daa_res_sub$comparison <- paste(daa_res_sub$group1, "VS", daa_res_sub$group2) 
  if(pAdjust){
      daa_res_sub$significance <- cut(daa_res_sub$p_adjust,
                                  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                  labels = c("<0.001", "<0.01", "<0.05", "n.s."),
                                  right = FALSE)
  }
  else{
    daa_res_sub$significance <- cut(daa_res_sub$p_value,
                                    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                    labels = c("<0.001", "<0.01", "<0.05", "n.s."),
                                    right = FALSE)
  }

  
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
  path_sub <- pathways[, c(TRUE, colnames(pathways)[2:ncol(pathways)] %in% meta$id[meta$gg_group2 %in% groups])] # Keep first col
  
  # Sub for metadata
  meta_sub <- meta[meta$gg_group2 %in% groups,]
  
  return(list(path_sub, meta_sub))
  
}

# Mapping for metacyc id and pathway names
path_names <- readLines("../../../picrust2 database/pathway_names.txt")
path_ids <- readLines("../../../picrust2 database/metacyc_ids.txt")
map_df <- data.frame(pathway_name = path_names, pathway = path_ids)

dfs <- df_preparation(groups = c("50:dss", "500:dss"))
customColors = c('darkred','darkblue')
plots <- daa_pathways_heatmap_pairwise(pathways = dfs[[1]], map_df = map_df, custom_colors = customColors, metadata = dfs[[2]], group = "gg_group2", threshold = 0.05, pAdjust = FALSE)
plots[[1]]
plots[[2]]


dfs <- df_preparation(groups = c("IL-22ra1-/-:Vehicle", "IL-22ra1-/-:Putrescine"))
customColors = c("#AB8F23","#04208D")
plots <- daa_pathways_heatmap_pairwise(pathways =  dfs[[1]], map_df = map_df, custom_colors = customColors, metadata =  dfs[[2]], group = "gg_group", threshold = 0.05)
plots[[1]]
plots[[2]]






