library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(readxl)
library(openxlsx)

source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/microbiota_analysis/utilities.R")
source("~/Documents/CHUM_git/gut-microbiota-iron/pipeline_linux/picrust2/picrust2_utilities.R")

# Metadata handling
{
  #loading metadata of interest
  metadata <- read.csv("~/Documents/CHUM_git/Microbiota_18/metadata/metadata.csv", sep = ";")
  
  # Remove the non-metadata stuff (liver measures and stuff)
  metadata <- metadata[,-c(5:8)]
  
  # Remove the letter at the end of id
  metadata$id <- substring(metadata$id, 1, 5)
  
  #adding id col as rownames too
  rownames(metadata) <- metadata$id
  
  # Remove dead mouse
  metadata <- metadata[-46,]
  
  # Extract 16S reads sample ids
  samples <- read.xlsx("~/Documents/CHUM_git/Microbiota_18/metadata/Microbiota_18_samples_2025-01-13.xlsx")
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
KOsToCompoundAbundanceGraphs <- function(ko_abundance, ko_annotations, metadata, group, customColors, sample_id_col, path, dim = c(6,6), additionalAes){

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
    
    if(isFALSE(is.null(additionalAes))){
      # p <- do.call("+", c(list(p), additionnalAes))
      
      p <- Reduce("+", c(list(p), additionalAes))
    }
    
    # Save plot
    ggsave(plot = p, filename = paste0(path, "/", clean_string(compound), "_predicted.png"), dpi = 300, width = dim[1], height = dim[2], bg = "white")
  }
}

KOsToCompoundAbundanceGraphs(ko_abundance = ko_metagenome[[1]],
                            ko_annotations = ko_annotations,
                            metadata = meta,
                            group = "gg_group2",
                            customColors = c("blue","red","darkblue","darkred"),
                            sample_id_col="id",
                            path = "~/Documents/CHUM_git/figures/thibault_new/picrust2",
                            dim = c(5,4),
                            additionalAes =
                              list(scale_x_discrete(labels = c("50 ppm\ncontrol","500 ppm\ncontrol","50 ppm\nDSS","500 ppm\nDSS")),
                                   theme(
                                     plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
                                     axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
                                     axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
                                     axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),  # Adjust x-axis tick label font size
                                     axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
                                     legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
                                     legend.text = element_text(size = 12),  # Adjust legend font size
                                     panel.grid.major = element_blank(),  # Add major grid lines
                                     panel.grid.minor = element_blank(),  # Remove minor grid lines
                                     axis.line = element_line(color = "black", size = 1)),
                                   labs(color = "", x="")))


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


# Look at unstrat results
setwd("~/Documents/CHUM_git/Microbiota_18/picrust2/picrust2_out_pipeline2/")
pred_ko_contrib <- read.table("KO_metagenome_out/pred_metagenome_contrib.tsv.gz", sep = "\t", header = TRUE)
ko_annotations <- readxl::read_excel("~/Documents/CHUM_git/picrust2 database/compound_to_kos2.xlsx")

# For butyrate - reference ko K00929
df <- pred_ko_contrib[pred_ko_contrib$function. == "K00929",]
df <- df[order(df$norm_taxon_function_contrib, decreasing = TRUE), ]
df[which.max(df$norm_taxon_function_contrib), ]
unique(df$taxon)

sub <- df[df$taxon == "ASV66",]






















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































# Utilizing ggpicrust2 tools to perform statistical analyis and find differentially abundant KOs
ko_abundance <- ko_metagenome[[1]]

# Load metadata (groups must be as factors)
meta_sub <- meta[colnames(ko_abundance), ] # reorder the metadata df with same order as colnames of the ko_abundance table

# Subsets of metadata and ko_abundance matrix for only dss groups
meta_sub <- meta_sub[meta_sub$treatment == "dss",]
ko_abundance <- ko_abundance[,meta_sub$id]

daa_res <- pathway_daa(abundance = ko_abundance, metadata = meta_sub, group = "diet", daa_method = "DESeq2", p.adjust = "fdr")
daa_res <- pathway_daa(abundance = ko_abundance, metadata = meta_sub, group = "diet", daa_method = "ALDEx2")
View(daa_res)


butyrate <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Butyrate1"]
butyrate <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Butyrate2"]
propionate <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate1"]
propionate <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate2"]
propionate <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate3"]
acetate <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Acetate"]
print(daa_res[daa_res$feature == butyrate,])
print(daa_res[daa_res$feature == propionate,])
print(daa_res[daa_res$feature == acetate,])




# Aldex2 workflow for analysis of picrust2 output
# Load KO abundance data
library(ALDEx2)
ko_abundance <- ko_metagenome[[1]]

# Round to nearest integer
ko_abundance <- round(ko_abundance)

# Load metadata (groups must be as factors)
meta_sub <- metadata[colnames(ko_abundance), ] # reorder the metadata df with same order as colnames of the ko_abundance table

# Subsets of metadata and ko_abundance matrix for only dss groups
meta_sub <- meta_sub[meta_sub$treatment == "dss",]
ko_abundance <- ko_abundance[,meta_sub$id]

# Create the model matrix
# model_matrix <- model.matrix(~ diet, data = meta_sub)

# # Create group column for CLR to work
# meta_sub$group <- meta_sub[,"gg_group2"] 
# meta_sub$group <- as.character(meta_sub$group)

meta_sub$group <- as.character(meta_sub$diet)

# Generate the clr-transformed values
aldex_clr <- aldex.clr(ko_abundance, conds = meta_sub$group, mc.samples = 128, denom = "all", verbose = FALSE, useMC = TRUE)
# aldex_clr <- aldex.clr(ko_abundance, conds = meta_sub$group, mc.samples = 128, denom = "all", verbose = FALSE, useMC = TRUE)

# Perform Welch’s t-test and Wilcoxon rank-sum test
aldex_tt_results <- aldex.ttest(aldex_clr, hist.plot = FALSE, paired.test = FALSE, verbose = FALSE)
View(aldex_tt_results["K00929",]) # Butyrate1
View(aldex_tt_results["K01034",]) # Butyrate2
View(aldex_tt_results["K01026",]) # Propionate1



# Perform the Kruskal-Wallis test
aldex_kw_results <- aldex.kw(aldex_clr, useMC = TRUE)

# Perform glm 
aldex_glm <- ALDEx2::aldex.glm(aldex_clr, model_matrix)

# Extract contrasts for pairwise comparisons
# For example, comparing 'treatmentdss' vs 'treatment50' within 'dietwater'
contrast_matrix <- makeContrasts(
  treatmentdss_vs_50_water = treatmentdss - treatment50,
  levels = colnames(model_matrix)
)

# Apply the contrast to the GLM results
pairwise_results <- coef(aldex_glm) %*% t(contrast_matrix)

butyrate <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Butyrate"]
propionate <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate"]
acetate <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Acetate"]
View(aldex_glm[butyrate,])
print(aldex_glm[propionate,])
print(aldex_glm[acetate,])


# Calculate effect sizes
aldex_effect <- aldex.glm.effect(aldex_clr)

























# Analysis of picrust2 results for timepoint t35
setwd("~/Documents/CHUM_git/Microbiota_18/picrust2/t35/picrust2/picrust2_out_pipeline/")
ko_metagenome <- read.table("KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE) # Load data

# Remove the X at the start of the names
colnames(ko_metagenome)[2:ncol(ko_metagenome)] = substring(colnames(ko_metagenome[2:ncol(ko_metagenome)]), first = 2, 6)
rownames(ko_metagenome) <- ko_metagenome$function.
ko_metagenome <- ko_metagenome[,-1]
ko_predicted <- read.table("KO_predicted.tsv.gz", sep = "\t", header = TRUE, row.names = 1) 

# Only t35  metadata
meta <- metadata[metadata$timepoint == "35",]
rownames(meta) <- meta$id
meta$diet <- factor(meta$diet, levels = c("50", "500"))
meta$treatment <- factor(meta$treatment, levels = c("water", "dss"))
meta$gg_group2 <- factor(meta$gg_group2, levels = c("50:water", "500:water", "50:dss", "500:dss"))

# Load file with KOs annotations for compound/pathways of interest
ko_annotations <- readxl::read_excel("~/Documents/CHUM_git/picrust2 database/compound_to_kos.xlsx")
ko_annotations <- readxl::read_excel("~/Documents/CHUM_git/picrust2 database/compound_to_kos2.xlsx")

# Graphs of SCFAs predicted abundance at t35
KOsToCompoundAbundanceGraphs(ko_abundance = ko_metagenome,
                             ko_annotations = ko_annotations,
                             metadata = meta,
                             group = "diet",
                             customColors = c("blue","red"),
                             sample_id_col="id",
                             path = "~/Documents/CHUM_git/figures/thibault_new/picrust2/t35",
                             dim = c(5,4),
                             additionalAes =
                               list(scale_x_discrete(labels = c("50 ppm","500 ppm")),
                                    theme(
                                      plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
                                      axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
                                      axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
                                      axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),  # Adjust x-axis tick label font size
                                      axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
                                      legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
                                      legend.text = element_text(size = 12),  # Adjust legend font size
                                      panel.grid.major = element_blank(),  # Add major grid lines
                                      panel.grid.minor = element_blank(),  # Remove minor grid lines
                                      axis.line = element_line(color = "black", size = 1)),
                                    labs(color = "", x="")))

# Stats with mixed effect models
{
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
}

# Utilizing ggpicrust2 tools to perform statistical analysis and find differentially abundant KOs
ko_abundance <- ko_metagenome

# Load metadata (groups must be as factors)
meta <- meta[colnames(ko_abundance), ] # reorder the metadata df with same order as colnames of the ko_abundance table
daa_res <- pathway_daa(abundance = ko_abundance, metadata = meta, group = "diet", daa_method = "DESeq2", p.adjust = "fdr")
daa_res <- pathway_daa(abundance = ko_abundance, metadata = meta, group = "diet", daa_method = "ALDEx2")
View(daa_res)

butyrate1 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Butyrate1"]
butyrate2 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Butyrate2"]
propionate1 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate1"]
propionate2 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate2"]
propionate3 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate3"]
acetate <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Acetate"]
print(daa_res[daa_res$feature == butyrate1,])
print(daa_res[daa_res$feature == butyrate2,])
print(daa_res[daa_res$feature == propionate1,])
print(daa_res[daa_res$feature == propionate2,])
print(daa_res[daa_res$feature == propionate3,])
print(daa_res[daa_res$feature == acetate,])



# Analysis of picrust2 results for timepoint t49
setwd("~/Documents/CHUM_git/Microbiota_18/picrust2/t49/picrust2/picrust2_out_pipeline/")
ko_metagenome <- read.table("KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE) # Load data

# Remove the X at the start of the names
colnames(ko_metagenome)[2:ncol(ko_metagenome)] = substring(colnames(ko_metagenome[2:ncol(ko_metagenome)]), first = 2, 6)
rownames(ko_metagenome) <- ko_metagenome$function.
ko_metagenome <- ko_metagenome[,-1]
ko_predicted <- read.table("KO_predicted.tsv.gz", sep = "\t", header = TRUE, row.names = 1) 

# Only t49  metadata
meta <- metadata[metadata$timepoint == "49",]
rownames(meta) <- meta$id
meta$diet <- factor(meta$diet, levels = c("50", "500"))
meta$treatment <- factor(meta$treatment, levels = c("water", "dss"))
meta$gg_group2 <- factor(meta$gg_group2, levels = c("50:water", "500:water", "50:dss", "500:dss"))

# Load file with KOs annotations for compound/pathways of interest
ko_annotations <- readxl::read_excel("~/Documents/CHUM_git/picrust2 database/compound_to_kos.xlsx")
ko_annotations <- readxl::read_excel("~/Documents/CHUM_git/picrust2 database/compound_to_kos2.xlsx")


# Test of kegg2ko abundance workflow
# Convert KO abundance to KEGG pathway abundance
setwd("~/Documents/CHUM_git/Microbiota_18/picrust2/t49/picrust2/picrust2_out_pipeline/")
kegg_abundance <- ko2kegg_abundance(file = "KO_metagenome_out/pred_metagenome_unstrat.tsv")
colnames(kegg_abundance)[1:ncol(kegg_abundance)] = substring(colnames(kegg_abundance[1:ncol(kegg_abundance)]), first = 1, 5)
butyrate <- t(kegg_abundance["ko00650",]) # https://www.kegg.jp/entry/map00650
propionate <- t(kegg_abundance["ko00640",])
colnames(butyrate) <- "abundance"
colnames(propionate) <- "abundance"
meta <- metadata[metadata$timepoint == "35",]
row.names(meta) <- meta$id
meta$diet <- factor(meta$diet, levels = c("50", "500"))
butyrate <- merge(butyrate, meta, by = "row.names")
propionate <- merge(propionate, meta, by = "row.names")

p <- ggplot(propionate, aes(x = diet, y = abundance, color = diet)) +
  geom_point(size = 1, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) + 
  
  #Error bars
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
               aes(color = diet),
               width = 0.2, size = 0.7,
               position = position_dodge(-0.75)) +
  
  #Mean lines
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
               aes(ymin = ..y.., ymax = ..y.., group = diet),
               color = "black", linewidth = 0.5, width = 0.5,
               position = position_dodge(-0.75))+
  
  #Connecting mean points with lines
  # stat_summary(fun = mean, geom = "line", size = 1.2) +  # Connecting means with lines
  
  
  
  labs(title = "Predicted propionate abundance", y = paste("Normalized Abundance"), color = "", x = "") +
  scale_color_manual(values = c("blue", "red"))+
  scale_y_continuous(limits = c(0, NA))+
  scale_x_discrete(labels = c("50 ppm","500 ppm"))+
  theme(
    plot.title = element_text(size = 13, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_text(size = 13, face = "bold"),  # Adjust x-axis label font size and style
    axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
    legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_blank(),  # Add major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),
    axis.line = element_line(color = "black", size = 1))

ggsave(plot = p, filename = "~/Documents/CHUM_git/figures/thibault_new/icm_seminar/butyrate_t35.png", bg = "white", width = 4, height = 3, dpi = 300)

# Statistics
# Load metadata (groups must be as factors)
meta <- meta[colnames(kegg_abundance), ] # reorder the metadata df with same order as colnames of the ko_abundance table
daa_res <- pathway_daa(abundance = kegg_abundance, metadata = meta, group = "diet", daa_method = "DESeq2", p.adjust = "fdr")
daa_res <- pathway_daa(abundance = kegg_abundance, metadata = meta, group = "diet", daa_method = "ALDEx2")
View(daa_res[daa_res$feature == "ko00650",])
View(daa_res[daa_res$feature == "ko00640",])



# Analysis of picrust2 results for timepoint t35
setwd("~/Documents/CHUM_git/Microbiota_18/picrust2/t35/picrust2/picrust2_out_pipeline/")
# Test of kegg2ko abundance workflow
# Convert KO abundance to KEGG pathway abundance
kegg_abundance <- ko2kegg_abundance(file = "KO_metagenome_out/pred_metagenome_unstrat.tsv")
colnames(kegg_abundance)[1:ncol(kegg_abundance)] = substring(colnames(kegg_abundance[1:ncol(kegg_abundance)]), first = 1, 5)
butyrate <- t(kegg_abundance["ko00650",]) # https://www.kegg.jp/entry/map00650
colnames(butyrate) <- "abundance"
butyrate <- merge(butyrate, meta, by = "row.names")

ggplot(butyrate, aes(x = diet, y = abundance, color = diet)) +
  geom_point(size = 1, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) + 
  
  #Error bars
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
               aes(color = diet),
               width = 0.2, size = 0.7,
               position = position_dodge(-0.75)) +
  
  #Mean lines
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
               aes(ymin = ..y.., ymax = ..y.., group = diet),
               color = "black", linewidth = 0.5, width = 0.5,
               position = position_dodge(-0.75))+
  
  #Connecting mean points with lines
  # stat_summary(fun = mean, geom = "line", size = 1.2) +  # Connecting means with lines
  
  
  
  labs(title = "Predicted butyrate abundance", y = paste("Normalized Abundance"), color = "", x = "") +
  scale_color_manual(values = c("blue", "red"))+
  scale_y_continuous(limits = c(0, NA))+
  scale_x_discrete(labels = c("50 ppm","500 ppm"))+
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
    axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
    legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_blank(),  # Add major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),
    axis.line = element_line(color = "black", size = 1))

# Statistics
# Load metadata (groups must be as factors)
meta <- meta[colnames(kegg_abundance), ] # reorder the metadata df with same order as colnames of the ko_abundance table
daa_res <- pathway_daa(abundance = kegg_abundance, metadata = meta, group = "diet", daa_method = "DESeq2", p.adjust = "fdr")
daa_res <- pathway_daa(abundance = ko_abundance, metadata = meta, group = "diet", daa_method = "ALDEx2")
View(daa_res[daa_res$feature == "ko00650",])














# Graphs of SCFAs predicted abundance at t49
KOsToCompoundAbundanceGraphs(ko_abundance = ko_metagenome,
                             ko_annotations = ko_annotations,
                             metadata = meta,
                             group = "diet",
                             customColors = c("blue","red"),
                             sample_id_col="id",
                             path = "~/Documents/CHUM_git/figures/thibault_new/picrust2/t49",
                             dim = c(5,4),
                             additionalAes =
                               list(scale_x_discrete(labels = c("50 ppm","500 ppm")),
                                    theme(
                                      plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
                                      axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
                                      axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
                                      axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),  # Adjust x-axis tick label font size
                                      axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
                                      legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
                                      legend.text = element_text(size = 12),  # Adjust legend font size
                                      panel.grid.major = element_blank(),  # Add major grid lines
                                      panel.grid.minor = element_blank(),  # Remove minor grid lines
                                      axis.line = element_line(color = "black", size = 1)),
                                    labs(color = "", x="")))

# Stats with mixed effect models
{
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
}

# Utilizing ggpicrust2 tools to perform statistical analysis and find differentially abundant KOs
ko_abundance <- ko_metagenome

# Load metadata (groups must be as factors)
meta <- meta[colnames(ko_abundance), ] # reorder the metadata df with same order as colnames of the ko_abundance table
daa_res <- pathway_daa(abundance = ko_abundance, metadata = meta, group = "diet", daa_method = "DESeq2", p.adjust = "fdr")
daa_res <- pathway_daa(abundance = ko_abundance, metadata = meta, group = "diet", daa_method = "ALDEx2")
View(daa_res)

butyrate1 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Butyrate1"]
butyrate2 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Butyrate2"]
propionate1 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate1"]
propionate2 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate2"]
propionate3 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate3"]
acetate <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Acetate"]
print(daa_res[daa_res$feature == butyrate1,])
print(daa_res[daa_res$feature == butyrate2,])
print(daa_res[daa_res$feature == propionate1,])
print(daa_res[daa_res$feature == propionate2,])
print(daa_res[daa_res$feature == propionate3,])
print(daa_res[daa_res$feature == acetate,])




# Analysis of picrust2 results for timepoint t54
setwd("~/Documents/CHUM_git/Microbiota_18/picrust2/t54/picrust2/picrust2_out_pipeline/")
ko_metagenome <- read.table("KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE) # Load data

# Remove the X at the start of the names
colnames(ko_metagenome)[2:ncol(ko_metagenome)] = substring(colnames(ko_metagenome[2:ncol(ko_metagenome)]), first = 2, 6)
rownames(ko_metagenome) <- ko_metagenome$function.
ko_metagenome <- ko_metagenome[,-1]
ko_predicted <- read.table("KO_predicted.tsv.gz", sep = "\t", header = TRUE, row.names = 1) 

# Only t54  metadata
meta <- metadata[metadata$timepoint == "54",]
meta <- meta[!meta$sample_id %in% c("10959_T54","33115_T54", # Remove samples for which t54 failed
                                    "33105_T54","10994_T54","10994_d53","10961_T54"),]

rownames(meta) <- meta$id
meta$diet <- factor(meta$diet, levels = c("50", "500"))
meta$treatment <- factor(meta$treatment, levels = c("water", "dss"))
meta$gg_group2 <- factor(meta$gg_group2, levels = c("50:water", "500:water", "50:dss", "500:dss"))

# Load file with KOs annotations for compound/pathways of interest
ko_annotations <- readxl::read_excel("~/Documents/CHUM_git/picrust2 database/compound_to_kos.xlsx")
ko_annotations <- readxl::read_excel("~/Documents/CHUM_git/picrust2 database/compound_to_kos2.xlsx")

# Graphs of SCFAs predicted abundance at t54
KOsToCompoundAbundanceGraphs(ko_abundance = ko_metagenome,
                             ko_annotations = ko_annotations,
                             metadata = meta,
                             group = "gg_group2",
                             customColors = c("blue","red", "darkblue", "darkred"),
                             sample_id_col="id",
                             path = "~/Documents/CHUM_git/figures/thibault_new/picrust2/t54",
                             dim = c(5,4),
                             additionalAes =
                               list(scale_x_discrete(labels = c("50:water", "500:water", "50:dss", "500:dss")),
                                    theme(
                                      plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
                                      axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
                                      axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
                                      axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),  # Adjust x-axis tick label font size
                                      axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
                                      legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
                                      legend.text = element_text(size = 12),  # Adjust legend font size
                                      panel.grid.major = element_blank(),  # Add major grid lines
                                      panel.grid.minor = element_blank(),  # Remove minor grid lines
                                      axis.line = element_line(color = "black", size = 1)),
                                    labs(color = "", x="")))

# Stats with mixed effect models
{
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
}

# Utilizing ggpicrust2 tools to perform statistical analysis and find differentially abundant KOs
ko_abundance <- ko_metagenome

# Load metadata (groups must be as factors)
meta <- meta[colnames(ko_abundance), ] # reorder the metadata df with same order as colnames of the ko_abundance table
daa_res <- pathway_daa(abundance = ko_abundance, metadata = meta, group = "diet", daa_method = "DESeq2", p.adjust = "fdr")
daa_res <- pathway_daa(abundance = ko_abundance, metadata = meta, group = "diet", daa_method = "ALDEx2")
View(daa_res)

butyrate1 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Butyrate1"]
butyrate2 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Butyrate2"]
propionate1 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate1"]
propionate2 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate2"]
propionate3 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate3"]
acetate <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Acetate"]
print(daa_res[daa_res$feature == butyrate1,])
print(daa_res[daa_res$feature == butyrate2,])
print(daa_res[daa_res$feature == propionate1,])
print(daa_res[daa_res$feature == propionate2,])
print(daa_res[daa_res$feature == propionate3,])
print(daa_res[daa_res$feature == acetate,])





# Analysis of picrust2 results for timepoint t49 with all groups
setwd("~/Documents/CHUM_git/Microbiota_18/picrust2/t49/picrust2/picrust2_out_pipeline/")
ko_metagenome <- read.table("KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", sep = "\t", header = TRUE) # Load data

# Remove the X at the start of the names
colnames(ko_metagenome)[2:ncol(ko_metagenome)] = substring(colnames(ko_metagenome[2:ncol(ko_metagenome)]), first = 2, 6)
rownames(ko_metagenome) <- ko_metagenome$function.
ko_metagenome <- ko_metagenome[,-1]
ko_predicted <- read.table("KO_predicted.tsv.gz", sep = "\t", header = TRUE, row.names = 1) 

# Only t54  metadata
meta <- metadata[metadata$timepoint == "49",]
rownames(meta) <- meta$id
meta$diet <- factor(meta$diet, levels = c("50", "500"))
meta$treatment <- factor(meta$treatment, levels = c("water", "dss"))
meta$gg_group2 <- factor(meta$gg_group2, levels = c("50:water", "500:water", "50:dss", "500:dss"))

# Load file with KOs annotations for compound/pathways of interest
ko_annotations <- readxl::read_excel("~/Documents/CHUM_git/picrust2 database/compound_to_kos.xlsx")
ko_annotations <- readxl::read_excel("~/Documents/CHUM_git/picrust2 database/compound_to_kos2.xlsx")

# Graphs of SCFAs predicted abundance at t54
KOsToCompoundAbundanceGraphs(ko_abundance = ko_metagenome,
                             ko_annotations = ko_annotations,
                             metadata = meta,
                             group = "gg_group2",
                             customColors = c("blue","red", "darkblue", "darkred"),
                             sample_id_col="id",
                             path = "~/Documents/CHUM_git/figures/thibault_new/picrust2/t49_all_groups",
                             dim = c(5,4),
                             additionalAes =
                               list(scale_x_discrete(labels = c("50:water", "500:water", "50:dss", "500:dss")),
                                    theme(
                                      plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
                                      axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
                                      axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
                                      axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),  # Adjust x-axis tick label font size
                                      axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
                                      legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
                                      legend.text = element_text(size = 12),  # Adjust legend font size
                                      panel.grid.major = element_blank(),  # Add major grid lines
                                      panel.grid.minor = element_blank(),  # Remove minor grid lines
                                      axis.line = element_line(color = "black", size = 1)),
                                    labs(color = "", x="")))

# Stats with mixed effect models
{
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
}

# Utilizing ggpicrust2 tools to perform statistical analysis and find differentially abundant KOs
ko_abundance <- ko_metagenome

# Load metadata (groups must be as factors)
meta <- meta[colnames(ko_abundance), ] # reorder the metadata df with same order as colnames of the ko_abundance table
daa_res <- pathway_daa(abundance = ko_abundance, metadata = meta, group = "diet", daa_method = "DESeq2", p.adjust = "fdr")
daa_res <- pathway_daa(abundance = ko_abundance, metadata = meta, group = "diet", daa_method = "ALDEx2")
View(daa_res)

butyrate1 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Butyrate1"]
butyrate2 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Butyrate2"]
propionate1 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate1"]
propionate2 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate2"]
propionate3 <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Propionate3"]
acetate <- ko_annotations$Kos[ko_annotations$`compound/pathway`== "Acetate"]
print(daa_res[daa_res$feature == butyrate1,])
print(daa_res[daa_res$feature == butyrate2,])
print(daa_res[daa_res$feature == propionate1,])
print(daa_res[daa_res$feature == propionate2,])
print(daa_res[daa_res$feature == propionate3,])
print(daa_res[daa_res$feature == acetate,])





























# Kegg pathways at last timepoint
setwd("~/Documents/CHUM_git/Microbiota_18/picrust2/tfinal/picrust2_out_pipeline2/")
setwd("~/Documents/CHUM_git/Microbiota_18/picrust2/t54/picrust2/picrust2_out_pipeline//")
# Test of kegg2ko abundance workflow
# Convert KO abundance to KEGG pathway abundance
kegg_abundance <- ko2kegg_abundance(file = "KO_metagenome_out/pred_metagenome_unstrat.tsv")
colnames(kegg_abundance)[1:ncol(kegg_abundance)] = substring(colnames(kegg_abundance[1:ncol(kegg_abundance)]), first = 1, 5)
butyrate <- t(kegg_abundance["ko00650",]) # https://www.kegg.jp/entry/map00650
propionate <- t(kegg_abundance["ko00640",])
colnames(butyrate) <- "abundance"
colnames(propionate) <- "abundance"
meta <- metadata[metadata$timepoint == "final",]
meta <- metadata[metadata$timepoint == "54",]
meta <- meta[!meta$sample_id %in% c("10959_T54","33115_T54", # Remove samples for which t54 failed
                                    "33105_T54","10994_T54","10994_d53","10961_T54"),]
rownames(meta) <- meta$id
meta$gg_group2 <- factor(meta$gg_group2, levels = c("50:water", "500:water", "50:dss", "500:dss"))
meta$diet <- factor(meta$diet, levels = c("50", "500"))
meta$treatment <- factor(meta$treatment, levels = c("water", "dss"))
butyrate <- merge(butyrate, meta, by = "row.names")
propionate <- merge(propionate, meta, by = "row.names")

ggplot(butyrate, aes(x = gg_group2, y = abundance, color = gg_group2)) +
  geom_point(size = 1, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) + 
  
  #Error bars
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
               aes(color = gg_group2),
               width = 0.2, size = 0.7,
               position = position_dodge(-0.75)) +
  
  #Mean lines
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
               aes(ymin = ..y.., ymax = ..y.., group = gg_group2),
               color = "black", linewidth = 0.5, width = 0.5,
               position = position_dodge(-0.75))+
  
  #Connecting mean points with lines
  # stat_summary(fun = mean, geom = "line", size = 1.2) +  # Connecting means with lines
  
  
  
  labs(title = "Predicted butyrate abundance", y = paste("Normalized Abundance"), color = "", x = "") +
  scale_color_manual(values = c("blue", "red", "darkblue", "darkred"))+
  scale_y_continuous(limits = c(0, NA))+
  scale_x_discrete(labels = c("50 ppm\nctrl", "500 ppm\nctrl", "50 ppm\nDSS", "500 ppm\nDSS"))+
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
    axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
    legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_blank(),  # Add major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),
    axis.line = element_line(color = "black", size = 1))

# Statistics, only between dss groups
# Load metadata (groups must be as factors)
meta_sub <- meta[meta$treatment == "dss",]
meta_sub <- meta[meta$treatment == "water",]
meta_sub <- meta[meta$diet == "50",]
meta_sub <- meta[meta$diet == "500",]
kegg_ab_sub <- kegg_abundance[,colnames(kegg_abundance) %in% unique(meta_sub$id)]
meta_sub <- meta_sub[colnames(kegg_ab_sub), ] # reorder the metadata df with same order as colnames of the ko_abundance table
daa_res <- pathway_daa(abundance = kegg_ab_sub, metadata = meta_sub, group = "diet", daa_method = "DESeq2", p.adjust = "fdr")
daa_res <- pathway_daa(abundance = kegg_ab_sub, metadata = meta_sub, group = "treatment", daa_method = "DESeq2", p.adjust = "fdr")
daa_res <- pathway_daa(abundance = kegg_ab_sub, metadata = meta_sub, group = "diet", daa_method = "ALDEx2")
daa_res <- pathway_daa(abundance = kegg_ab_sub, metadata = meta_sub, group = "treatment", daa_method = "ALDEx2")
print(daa_res[daa_res$feature == "ko00650",])
print(daa_res[daa_res$feature == "ko00640",])
