#For design with 4 groups based on 2 conditions - this the latest version used
#gg_group must be order with correct order prior to that (as a factor)
correlationGroups <- function(ps, deseq, measure = "log2fold", gg_group, taxa = "Species", threshold = 0.01, displayPvalue = FALSE, customColors, pairs, path, df, global = TRUE, showIndivCor = FALSE, transformation = "CLR", displayOnlySig = FALSE, saveFig=TRUE, displaySpeciesASVNumber = TRUE){
  
  # Checks if gg_group is a factor
  checkIfFactor(sample_data(ps)[[gg_group]])
  
  #Creates directory for taxonomic level
  dir <- paste(path, taxa, sep = "")
  existingDirCheck(path = dir)
  
  #Define empty list that will contain pairs comparaisons names
  vs <- c()

  #Save pairs comparaisons names for displaying them in the final sigtable
  for(k in seq_along(pairs)){
    
    #Save pair
    pair <- unlist(pairs[k])
  
    #Save specific pair comparaison
    vs <- c(vs, paste(clean_string(pair[1]), "vs", clean_string(pair[2]), sep="_"))
    
    }
  
  #Partition results for specific pairwise comparaisons
  res_subset1 <- results(deseq, name = resultsNames(deseq)[3]) #wt putrescine vs vehicle
  sigtab_1 <- cbind(as(res_subset1, "data.frame"), as(tax_table(ps)[rownames(res_subset1), ], "matrix"))
  sigtab_1$comparaison <- 1
  sigtab_1$vs <- vs[1]

  res_subset2 <- results(deseq, contrast=list(c(resultsNames(deseq)[3], resultsNames(deseq)[4]))) #il22 ko putrescine vs vehicle
  sigtab_2 <- cbind(as(res_subset2, "data.frame"), as(tax_table(ps)[rownames(res_subset2), ], "matrix"))
  sigtab_2$comparaison <- 2
  sigtab_2$vs <- vs[2]

  res_subset3 <- results(deseq, name=resultsNames(deseq)[2]) #vehicle wt vs il22 ko
  sigtab_3 <- cbind(as(res_subset3, "data.frame"), as(tax_table(ps)[rownames(res_subset3), ], "matrix"))
  sigtab_3$comparaison <- 3
  sigtab_3$vs <- vs[3]
  
  res_subset4 <- results(deseq, contrast=list(c(resultsNames(deseq)[2], resultsNames(deseq)[4]))) #putrescine wt vs il22 ko
  sigtab_4 <- cbind(as(res_subset4, "data.frame"), as(tax_table(ps)[rownames(res_subset4), ], "matrix"))
  sigtab_4$comparaison <- 4
  sigtab_4$vs <- vs[4]

  #Append the sigtabs together
  sigtab <- bind_rows(sigtab_1, sigtab_2, sigtab_3, sigtab_4)
  
  #Add ASV variable col to sigtab (enables to store asv names, not only as rownames, because they will be changed when using rowbind)
  sigtab["asv"] <- gsub("\\..*", "", rownames(sigtab))

  #Keeping only ASVs for which they were taxa found at the taxonomical level of interest
  sigtab <- sigtab[!is.na(sigtab[[taxa]]),]

  #Replacing NA padj by 1 (they correspond to this anyways)
  sigtab$padj[is.na(sigtab$padj)] <- 1

  #Add column that adds symbols for the significance
  # Define significance levels
  sigtab$significance <- as.character(cut(sigtab$padj,
                                          breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                          labels = c("***", "**", "*", "NS")))
  
  #Find asvs that have at least one significant value at some timepoint
  asvList <- unique(sigtab[(sigtab$padj)<threshold,"asv"])
  
  #Check if ASV list is empty and stop here if it is
  if(is_empty(asvList)){
    
    message(paste("No differential abundant taxa found for ",taxa))
    return(NULL)
  }
  
  ps <- transformCounts(ps, transformation = transformation)
  counts <- otu_table(ps)
  
  # # Ensure that the result is an otu_table object
  # otu_table(ps) <- otu_table(otu_matrix, taxa_are_rows = TRUE)
  
  #Keep only ASVs of interest (those that were significantly different)
  counts <- counts[rownames(counts) %in% asvList, ]
  
  #Transpose relab dataframe to get ASVs as columns and sample_id as rows
  counts <- t(counts)
  
  # Bind dataframe with values for correlation and relative abundance
  data_combined <- merge(counts, df, by = "row.names", sort = FALSE)
  
  if(global){
    
    #Creates directory for groups at which we are looking
    dir <- paste(dir, "/all_groups", sep = "")
    existingDirCheck(path = dir)
    
    # Set row names and drop the "Row.names" column
    rownames(data_combined) <- data_combined$Row.names
    data_combined <- data_combined[,-1]
    
    # Extract the ASV and variable matrices
    asv_matrix <- as.matrix(data_combined[, asvList, drop = FALSE]) # ASVs
    var_matrix <- as.matrix(data_combined[, !colnames(data_combined) %in% asvList, drop = FALSE]) # Non-ASV variables
  
    # Calculate correlation matrix
    cor_results <- rcorr(as.matrix(data_combined), type = "spearman")
    cor_matrix <- cor_results$r  # Extract correlation coefficients
    p_values <- cor_results$P    # Extract p-values
    
    # Adjust p-values using FDR
    p_adjusted <- matrix(
      p.adjust(as.vector(p_values), method = "fdr"),
      nrow = nrow(p_values),
      ncol = ncol(p_values),
      dimnames = dimnames(p_values))
    
    cor_matrix <- cor_matrix[rownames(cor_matrix) %in% asvList, !colnames(cor_matrix) %in% asvList]
    p_values <- p_adjusted[rownames(p_adjusted) %in% asvList, !colnames(p_adjusted) %in% asvList]
  
    # Melt correlation and adjusted p-value matrices
    cor_melt <- melt(cor_matrix, varnames = c("ASV", "Variables"))
    p_melt <- melt(p_values, varnames = c("ASV", "Variables"))
    
    # Combine melted data and filter out NA correlations (if any)
    cor_melt$p_value <- p_melt$value
    cor_melt <- cor_melt[!is.na(cor_melt$value), ]  # Remove NA rows
    cor_melt$significance <- ifelse(cor_melt$p_value < 0.001, "***", 
                                    ifelse(cor_melt$p_value < 0.01, "**", 
                                           ifelse(cor_melt$p_value < 0.05, "*", "")))
    
    # Add taxa_name and taxa_name_asv variable (for species only)
    taxonomy <- as.data.frame(tax_table(ps))
    taxonomy <- taxonomy[asvList,]
    
    # Retrieve taxonomic information
    cor_melt <- merge(cor_melt, taxonomy, by.x = "ASV", by.y = "row.names")
    
    if(taxa=="Species"){
      cor_melt$taxa_name <- paste(cor_melt$Genus,cor_melt$Species)
      if(displaySpeciesASVNumber){
        cor_melt$taxa_name <- paste(cor_melt$taxa_name, " (", cor_melt$ASV, ")", sep = "")
      }
      
    # Display only ASVs that have at least one significant correlation, by subsetting cor_melt
    if(displayOnlySig){
      asvList <- as.character(unique(cor_melt$ASV[cor_melt$significance != ""]))
      cor_melt <- cor_melt[cor_melt$ASV %in% asvList,]
    }
    
    # Create the heatmap
    p <- ggplot(cor_melt, aes(x = Variables, y = taxa_name, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                           midpoint = 0, limit = c(-1, 1), space = "Lab",
                           name = "Correlation") +
      geom_text(aes(label = significance), color = "black") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(face = "italic"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold")) +
      labs(x = "Variables", y = taxa)
    
    # Save correlation dataframe with stats and correlation coefficients 
    write_xlsx(cor_melt, path = paste0(dir,"/all_groups_correlation.xlsx"))
    
    if(saveFig){
      ggsave(plot = p, filename = paste(dir,"/correlation_all_groups.png", sep = ""), dpi = 300, height = 6, width = 10, bg = 'white')
    }
    else{
      return(p)
    }
    
    if (showIndivCor) {
      # Create directory for individual scatterplots
      dir_indiv <- paste0(dir, "/scatterplots/")
      existingDirCheck(path = dir_indiv)

      # Iterate over each ASV and variable pair (or other each sigAsvs if we are only interested in those)
      for (asv in asvList){
        for (var in setdiff(colnames(df), asvList)) {

          # Prepare data for scatterplot
          scatter_data <- data_combined %>%
            dplyr::select(all_of(c(asv, var))) %>%
            dplyr::rename(counts_var = all_of(asv), variable = all_of(var))
          
          # Convert row names to column in scatter_data
          scatter_data$sample_id <- row.names(scatter_data)
          
          # Extract gg_group from sample_data(ps)
          gg_group_data <- data.frame(sample_id = row.names(sample_data(ps)), 
                                      group = sample_data(ps)[[gg_group]])
          
          # Merge with scatter_data
          scatter_data <- merge(scatter_data, gg_group_data, by = "sample_id", all.x = TRUE)
          
          scatter_data$group <- factor(scatter_data$group, levels = levels(sample_data(ps)[[gg_group]]))
          
          # Get taxa name (if taxonomic level of interest is not species)
          if(taxa != "Species"){
            taxon <-  as.data.frame(tax_table(ps))[asv,taxa]
          }
          
          # Skip if there are no data points
          if (nrow(scatter_data) == 0) next
          
          # Create scatterplot
          scatter_plot <- ggplot(scatter_data, aes(x = variable, y = counts_var, color = group)) + # You can do log(variable) for less noisy look at the relationship between the variables
            geom_point() +
            geom_smooth(method = "lm", color = "red", se = TRUE) +
            theme_minimal() +
            labs(
              title = paste("Scatterplot:", asv, "vs", var),
              x = var,
              y = paste(transformation, "transformed abundance")
            )+
            scale_color_manual(values = customColors)
          
          # Save scatterplot
          ggsave(
            filename = paste0(dir_indiv, "scatter_", ifelse(taxa == "Species", asv, taxon), "_vs_", var, ".png"),
            plot = scatter_plot,
            width = 8, height = 6, dpi = 300,
            bg = "white"
          )
        }
      }
      message("Scatterplots saved to: ", dir_indiv)
    }
    
    
  }else{ #Produce multiple heatmaps for each group
    
    for(group in levels(sample_data(ps)[[gg_group]])){
      
      #Creates directory for groups at which we are looking
      dir_group <- paste(dir, "/", clean_string(group), sep = "")
      existingDirCheck(path = dir_group)
      
      # Set row names as sample_ids (Row.names)
      rownames(data_combined) <- data_combined$Row.names
      
      #Save list of samples linked to each group and creates prepare matrices for correlation
      samples <- sample_data(ps)[sample_data(ps)[[gg_group]] == group, "sample_id"]
      samples <- samples$sample_id
      data_grp <- data_combined[as.character(data_combined[,"Row.names"]) %in% samples, ]
      data_grp <- data_grp[,-1] #Remove first col which is Row.names because it is not a variable for the correlation
      asv_matrix <- as.matrix(data_grp[, asvList, drop = FALSE]) # ASVs
      var_matrix <- as.matrix(data_grp[, !colnames(data_grp) %in% asvList, drop = FALSE]) # Non-ASV variables
      
      # Calculate correlation matrix
      cor_results <- rcorr(as.matrix(data_grp), type = "spearman")
      cor_matrix <- cor_results$r  # Extract correlation coefficients
      p_values <- cor_results$P    # Extract p-values
      
      # print(cor_results)

      # Adjust p-values using FDR
      p_adjusted <- matrix(
        p.adjust(as.vector(p_values), method = "fdr"),
        nrow = nrow(p_values),
        ncol = ncol(p_values),
        dimnames = dimnames(p_values))

      cor_matrix <- cor_matrix[rownames(cor_matrix) %in% asvList, !colnames(cor_matrix) %in% asvList]
      p_values <- p_adjusted[rownames(p_adjusted) %in% asvList, !colnames(p_adjusted) %in% asvList]

      # Melt correlation and adjusted p-value matrices
      cor_melt <- melt(cor_matrix, varnames = c("ASV", "Variables"))
      p_melt <- melt(p_values, varnames = c("ASV", "Variables"))

      # Combine melted data and filter out NA correlations (if any)
      cor_melt$p_value <- p_melt$value
      cor_melt <- cor_melt[!is.na(cor_melt$value), ]  # Remove NA rows
      cor_melt$significance <- ifelse(cor_melt$p_value < 0.001, "***",
                                      ifelse(cor_melt$p_value < 0.01, "**",
                                             ifelse(cor_melt$p_value < 0.05, "*", "")))

      # Add taxa_name and taxa_name_asv variable (for species only)
      taxonomy <- as.data.frame(tax_table(ps))
      taxonomy <- taxonomy[asvList,]

      if(taxa=="Species"){
        cor_melt$taxa_name <- paste(taxonomy[cor_melt$ASV, "Genus"],taxonomy[cor_melt$ASV, "Species"])
        cor_melt$taxa_name <- paste(cor_melt$taxa_name, " (", cor_melt$ASV, ")", sep = "")
      }else{
        cor_melt$taxa_name <- taxonomy[cor_melt$ASV, taxa]
      }
      
      # Create the heatmap
      p <- ggplot(cor_melt, aes(x = Variables, y = taxa_name, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                             midpoint = 0, limit = c(-1, 1), space = "Lab",
                             name = "Correlation") +
        geom_text(aes(label = significance), color = "black") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_text(face = "italic"),
              axis.title.x = element_text(size = 12, face = "bold"),
              axis.title.y = element_text(size = 12, face = "bold")) +
        labs(x = "Variables", y = taxa)

      ggsave(plot = p, filename = paste0(dir_group,"/",clean_string(group),"_group_correlation.png"), dpi = 300, height = 6, width = 10, bg = 'white')
      
      # Save correlation dataframe with stats and correlation coefficients 
      write_xlsx(cor_melt, path = paste0(dir_group,"/",clean_string(group),"_group_correlation.xlsx"))
    }

  }
  
}
}


#For design with like 2 groups but multiple timepoints
#the time variable must be ordered as a factor
#identifying significant ASVs based on relabTimelineRevised, meaning deseq analysis is done within the function and then results are combined to identify significant ASVs across multiple timepoints
#Logically, we should only be using last timepoint!
# DEPRECIATED 
correlationTimepoints <- function(ps, measure = "log2fold", timeVariable, varToCompare, taxa = "Species", threshold = 0.01, displayPvalue = FALSE, customColors, path, df, global = TRUE){
  
  #Creates directory for taxonomic level
  dir <- paste(path, taxa, sep = "")
  existingDirCheck(path = dir)
  
  #Create empty sigtabCombined dataframe to merge sigtabs together
  sigtabCombined <- data.frame()
  
  #Create empty normalized_counts list to store normalized_counts matrices
  nrm_counts_all <- list()
  
  #Iterate through timepoints
  for(timePoint in levels(sample_data(ps)[[timeVariable]])){
    
    #Creating phyloseq objects for each timepoint
    ps_subset <- prune_samples(sample_data(ps)[[timeVariable]] == timePoint, ps)
    
    #Simple deseq object only accounting for the differences in varToCompare
    deseq_subset <- phyloseq_to_deseq2(ps_subset, as.formula(paste("~", varToCompare))) 
    
    #Performing the deseq analysis
    deseq_subset <- DESeq(deseq_subset, test="Wald", fitType = "parametric")
    
    #Get normalized counts from DESeq2 object
    normalized_counts <- counts(deseq_subset, normalized = TRUE)
    
    #Add the single normalized count to df to df with all
    nrm_counts_all <- append(nrm_counts_all, list(normalized_counts))
    
    #Save results at single timepoint
    res <- results(deseq_subset, name = resultsNames(deseq_subset)[2])
    
    #Save significance table
    sigtab <- cbind(as(res, "data.frame"), as(tax_table(ps_subset)[rownames(res),], "matrix"))
    
    #Add timeVariable column to the sigtab
    sigtab[[timeVariable]] <- timePoint
    
    #Add ASV variable col to sigtab (enables to store asv names, not only as rownames, because they will be changed when using rowbind)
    sigtab["asv"] <- rownames(sigtab)
    
    #Replacing NA padj by 1 (they correspond to this anyways)
    sigtab$padj[is.na(sigtab$padj)] <- 1
    
    #Add column that adds symbols for the significance 
    # Define significance levels
    sigtab$significance <- cut(sigtab$padj,
                               breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                               labels = c("***", "**", "*", "NS"))
    
    #Append timepoint-specific dataframe to combined dataframe
    sigtabCombined <- rbind(sigtabCombined, sigtab)
    
  }
  
  #Keeping only ASVs for which they were taxa found at the taxonomical level of interest
  sigtabCombined <- sigtabCombined[!is.na(sigtabCombined[[taxa]]),]
  
  #Find asvs that have at least one significant value at some timepoint
  asvList <- unique(sigtabCombined[(sigtabCombined$padj)<threshold,"asv"])
  
  #Check if ASV list is empty and stop here if it is
  if(is_empty(asvList)){
    
    message(paste("No differential abundant taxa found for ",taxa))
    return(NULL)
  }
  
  # print(normalized_counts)
  
  #Calculate relative abundance as a percentage for each sample
  # relative_abundance <- normalized_counts * 100 / colSums(normalized_counts)
  
  # Apply prop.table for relative abundance calculation and multiply by 100 for percentage
  relative_abundance <- apply(normalized_counts, 2, prop.table) * 100
  
  # # Ensure that the result is an otu_table object
  # otu_table(ps) <- otu_table(otu_matrix, taxa_are_rows = TRUE)
  
  #Keep only ASVs of interest (those that were significantly different)
  relative_abundance <- relative_abundance[rownames(relative_abundance) %in% asvList, ]
  
  #Transpose relab dataframe to get ASVs as columns and sample_id as rows
  relative_abundance <- t(relative_abundance)
  
  # Bind dataframe with values for correlation and relative abundance
  data_combined <- merge(relative_abundance, df, by = "row.names", sort = FALSE)
  
  if(global){
    
    # Set row names and drop the "Row.names" column
    rownames(data_combined) <- data_combined$Row.names
    data_combined <- data_combined[,-1]
    
    # Extract the ASV and variable matrices
    asv_matrix <- as.matrix(data_combined[, asvList, drop = FALSE]) # ASVs
    var_matrix <- as.matrix(data_combined[, !colnames(data_combined) %in% asvList, drop = FALSE]) # Non-ASV variables
    
    # Calculate correlation matrix
    cor_results <- rcorr(as.matrix(data_combined), type = "spearman")
    cor_matrix <- cor_results$r  # Extract correlation coefficients
    p_values <- cor_results$P    # Extract p-values
    
    # Adjust p-values using FDR
    p_adjusted <- matrix(
      p.adjust(as.vector(p_values), method = "fdr"),
      nrow = nrow(p_values),
      ncol = ncol(p_values),
      dimnames = dimnames(p_values))
    
    cor_matrix <- cor_matrix[rownames(cor_matrix) %in% asvList, !colnames(cor_matrix) %in% asvList]
    p_values <- p_adjusted[rownames(p_adjusted) %in% asvList, !colnames(p_adjusted) %in% asvList]
    
    # Melt correlation and adjusted p-value matrices
    cor_melt <- melt(cor_matrix, varnames = c("ASV", "Variables"))
    p_melt <- melt(p_values, varnames = c("ASV", "Variables"))
    
    # Combine melted data and filter out NA correlations (if any)
    cor_melt$p_value <- p_melt$value
    cor_melt <- cor_melt[!is.na(cor_melt$value), ]  # Remove NA rows
    cor_melt$significance <- ifelse(cor_melt$p_value < 0.001, "***", 
                                    ifelse(cor_melt$p_value < 0.01, "**", 
                                           ifelse(cor_melt$p_value < 0.05, "*", "")))
    
    # Add taxa_name and taxa_name_asv variable (for species only)
    taxonomy <- as.data.frame(tax_table(ps))
    taxonomy <- taxonomy[asvList,]
    
    if(taxa=="Species"){
      cor_melt$taxa_name <- paste(taxonomy[cor_melt$ASV, "Genus"],taxonomy[cor_melt$ASV, "Species"])
      cor_melt$taxa_name <- paste(cor_melt$taxa_name, " (", cor_melt$ASV, ")", sep = "")
    }else{
      cor_melt$taxa_name <- taxonomy[cor_melt$ASV, taxa]
    }
    
    # Create the heatmap
    p <- ggplot(cor_melt, aes(x = Variables, y = taxa_name, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                           midpoint = 0, limit = c(-1, 1), space = "Lab",
                           name = "Correlation") +
      geom_text(aes(label = significance), color = "black") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(face = "italic"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold")) +
      labs(x = "Variables", y = taxa)
    
    ggsave(plot = p, filename = paste(dir,"/correlation_all_groups.png", sep = ""), dpi = 300, height = 6, width = 10, bg = 'white')
    
  }else{ #Produce multiple heatmaps for each group
    
    for(group in levels(sample_data(ps)[[gg_group]])){
      
      # Set row names as sample_ids (Row.names)
      rownames(data_combined) <- data_combined$Row.names
      
      #Save list of samples linked to each group and creates prepare matrices for correlation
      samples <- sample_data(ps)[sample_data(ps)[[gg_group]] == group, "sample_id"]
      samples <- samples$sample_id
      data_grp <- data_combined[as.character(data_combined[,"Row.names"]) %in% samples, ]
      data_grp <- data_grp[,-1] #Remove first col which is Row.names because it is not a variable for the correlation
      asv_matrix <- as.matrix(data_grp[, asvList, drop = FALSE]) # ASVs
      var_matrix <- as.matrix(data_grp[, !colnames(data_grp) %in% asvList, drop = FALSE]) # Non-ASV variables
      
      # Calculate correlation matrix
      cor_results <- rcorr(as.matrix(data_grp), type = "spearman")
      cor_matrix <- cor_results$r  # Extract correlation coefficients
      p_values <- cor_results$P    # Extract p-values
      
      # print(cor_results)
      
      # Adjust p-values using FDR
      p_adjusted <- matrix(
        p.adjust(as.vector(p_values), method = "fdr"),
        nrow = nrow(p_values),
        ncol = ncol(p_values),
        dimnames = dimnames(p_values))
      
      cor_matrix <- cor_matrix[rownames(cor_matrix) %in% asvList, !colnames(cor_matrix) %in% asvList]
      p_values <- p_adjusted[rownames(p_adjusted) %in% asvList, !colnames(p_adjusted) %in% asvList]
      
      # Melt correlation and adjusted p-value matrices
      cor_melt <- melt(cor_matrix, varnames = c("ASV", "Variables"))
      p_melt <- melt(p_values, varnames = c("ASV", "Variables"))
      
      # Combine melted data and filter out NA correlations (if any)
      cor_melt$p_value <- p_melt$value
      cor_melt <- cor_melt[!is.na(cor_melt$value), ]  # Remove NA rows
      cor_melt$significance <- ifelse(cor_melt$p_value < 0.001, "***",
                                      ifelse(cor_melt$p_value < 0.01, "**",
                                             ifelse(cor_melt$p_value < 0.05, "*", "")))
      
      # Add taxa_name and taxa_name_asv variable (for species only)
      taxonomy <- as.data.frame(tax_table(ps))
      taxonomy <- taxonomy[asvList,]
      
      if(taxa=="Species"){
        cor_melt$taxa_name <- paste(taxonomy[cor_melt$ASV, "Genus"],taxonomy[cor_melt$ASV, "Species"])
        cor_melt$taxa_name <- paste(cor_melt$taxa_name, " (", cor_melt$ASV, ")", sep = "")
      }else{
        cor_melt$taxa_name <- taxonomy[cor_melt$ASV, taxa]
      }
      
      # Create the heatmap
      p <- ggplot(cor_melt, aes(x = Variables, y = taxa_name, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                             midpoint = 0, limit = c(-1, 1), space = "Lab",
                             name = "Correlation") +
        geom_text(aes(label = significance), color = "black") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_text(face = "italic"),
              axis.title.x = element_text(size = 12, face = "bold"),
              axis.title.y = element_text(size = 12, face = "bold")) +
        labs(x = "Variables", y = taxa)
      
      ggsave(plot = p, filename = paste0(dir,"/",clean_string(group),"_group_correlation.png"), dpi = 300, height = 6, width = 10, bg = 'white')
    }
    
  }
  
}






correlation2Var <- function(ps, deseq, measure = "log2fold", varToCompare, taxa = "Species", threshold = 0.01, displayPvalue = FALSE, path, df, global = TRUE, showIndivCor = FALSE, transformation = "CLR", displayOnlySig = FALSE, returnMainFig = FALSE, displaySpeciesASVNumber = TRUE, colorsHmap = c("blue","red")){
  
  #Creates directory for taxonomic level
  dir <- paste(path, taxa, sep = "")
  existingDirCheck(path = dir)
  
  # Transform counts (CLR transformation by default because improves correlation results)
  ps <- transformCounts(ps, transformation = transformation)
  counts <- otu_table(ps)
  
  # #Define empty list that will contain pairs comparaisons names
  # vs <- c()
  # 
  # #Save pairs comparaisons names for displaying them in the final sigtable
  # for(k in seq_along(pairs)){
  #   
  #   #Save pair
  #   pair <- unlist(pairs[k])
  #   
  #   #Save specific pair comparaison
  #   vs <- c(vs, paste(clean_string(pair[1]), "vs", clean_string(pair[2]), sep="_"))
  #   
  # }
  
  #Partition results for specific pairwise comparaisons
  res <- results(deseq, name = resultsNames(deseq)[2]) #50 vs 500
  sigtab <- cbind(as(res, "data.frame"), as(tax_table(ps)[rownames(res), ], "matrix"))
  
  #Add ASV variable col to sigtab (enables to store asv names, not only as rownames, because they will be changed when using rowbind)
  sigtab["asv"] <- gsub("\\..*", "", rownames(sigtab))
  
  #Keeping only ASVs for which they were taxa found at the taxonomical level of interest
  sigtab <- sigtab[!is.na(sigtab[[taxa]]),]
  
  #Replacing NA padj by 1 (they correspond to this anyways)
  sigtab$padj[is.na(sigtab$padj)] <- 1
  
  #Add column that adds symbols for the significance
  # Define significance levels
  sigtab$significance <- as.character(cut(sigtab$padj,
                                          breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                          labels = c("***", "**", "*", "NS")))
  
  #Find asvs that have at least one significant value at some timepoint
  asvList <- unique(sigtab[(sigtab$padj)<threshold,"asv"])
  
  #Check if ASV list is empty and stop here if it is
  if(is_empty(asvList)){
    
    message(paste("No differential abundant taxa found for ",taxa))
    return(NULL)
  }
  
  # print(normalized_counts)
  
  #Calculate relative abundance as a percentage for each sample
  # relative_abundance <- normalized_counts * 100 / colSums(normalized_counts)
  # if(normalizedCountsOnly){
  #   
  #   relative_abundance <- normalized_counts
  #   message("Using only normalized counts and not transforming them into relative abundance.")
  # }else{
  #   
  #   # Apply prop.table for relative abundance calculation and multiply by 100 for percentage
  #   relative_abundance <- apply(normalized_counts, 2, prop.table) * 100
  # }
  
  
  # # Ensure that the result is an otu_table object
  # otu_table(ps) <- otu_table(otu_matrix, taxa_are_rows = TRUE)
  
  #Keep only ASVs of interest (those that were significantly different)
  counts <- counts[rownames(counts) %in% asvList,]
  
  # #Transpose relab dataframe to get ASVs as columns and sample_id as rows
  counts <- t(counts)
  
  # Bind dataframe with values for correlation and relative abundance
  data_combined <- merge(counts, df, by = "row.names", sort = FALSE)
  
  if(global){
    
    #Creates directory for groups at which we are looking
    dir <- paste(dir, "/all_groups", sep = "")
    existingDirCheck(path = dir)
    
    # Set row names and drop the "Row.names" column
    rownames(data_combined) <- data_combined$Row.names
    data_combined <- data_combined[,-1]
    
    # Extract the ASV and variable matrices
    asv_matrix <- as.matrix(data_combined[, asvList, drop = FALSE]) # ASVs
    var_matrix <- as.matrix(data_combined[, !colnames(data_combined) %in% asvList, drop = FALSE]) # Non-ASV variables
    
    # Calculate correlation matrix
    cor_results <- rcorr(as.matrix(data_combined), type = "spearman")
    cor_matrix <- cor_results$r  # Extract correlation coefficients
    p_values <- cor_results$P    # Extract p-values
    
    # Adjust p-values using FDR
    p_adjusted <- matrix(
      p.adjust(as.vector(p_values), method = "fdr"),
      nrow = nrow(p_values),
      ncol = ncol(p_values),
      dimnames = dimnames(p_values))
    
    cor_matrix <- cor_matrix[rownames(cor_matrix) %in% asvList, !colnames(cor_matrix) %in% asvList]
    p_values <- p_adjusted[rownames(p_adjusted) %in% asvList, !colnames(p_adjusted) %in% asvList]
    
    # If only one sifgnificant ASV, fix or else the code does not work (doesnt work for now)
    if (length(asvList) == 1) {
      cor_matrix <- as.data.frame(cor_matrix)
      cor_matrix$ASV <- rownames(cor_matrix)
      cor_melt <- pivot_longer(cor_matrix, cols = -ASV, names_to = "Variables", values_to = "value")
    } else {
      cor_melt <- melt(cor_matrix, varnames = c("ASV", "Variables"))
    }
    
    
    # # Melt correlation and adjusted p-value matrices
    # cor_melt <- melt(cor_matrix, varnames = c("ASV", "Variables"))
    p_melt <- melt(p_values, varnames = c("ASV", "Variables"))
    
    # Combine melted data and filter out NA correlations (if any)
    cor_melt$p_value <- p_melt$value
    cor_melt <- cor_melt[!is.na(cor_melt$value), ]  # Remove NA rows
    cor_melt$significance <- ifelse(cor_melt$p_value < 0.001, "***", 
                                    ifelse(cor_melt$p_value < 0.01, "**", 
                                           ifelse(cor_melt$p_value < 0.05, "*", "")))
    
    # Add taxa_name and taxa_name_asv variable (for species only)
    taxonomy <- as.data.frame(tax_table(ps))
    taxonomy <- taxonomy[asvList,]
    
    # Retrieve taxonomic information
    cor_melt <- merge(cor_melt, taxonomy, by.x = "ASV", by.y = "row.names")
    
    if(taxa=="Species"){
      cor_melt$taxa_name <- paste(cor_melt$Genus,cor_melt$Species)
      if(displaySpeciesASVNumber){
        cor_melt$taxa_name <- paste(cor_melt$taxa_name, " (", cor_melt$ASV, ")", sep = "")
      }
    
    # if(taxa=="Species"){
    #   cor_melt$taxa_name <- paste(taxonomy[cor_melt$ASV, "Genus"],taxonomy[cor_melt$ASV, "Species"])
    #   if(displaySpeciesASVNumber){
    #     cor_melt$taxa_name <- paste(cor_melt$taxa_name, " (", cor_melt$ASV, ")", sep = "")
    #   }
      
    }else{
      cor_melt$taxa_name <- taxonomy[cor_melt$ASV, taxa]
    }
    
    # Display only ASVs that have at least one significant correlation, by subsetting cor_melt
    if(displayOnlySig){
      asvList <- as.character(unique(cor_melt$ASV[cor_melt$significance != ""]))
      cor_melt <- cor_melt[cor_melt$ASV %in% asvList,]
    }
    
    # Create the heatmap
    p <- ggplot(cor_melt, aes(x = Variables, y = taxa_name, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = colorsHmap[1], high = colorsHmap[2], mid = "white",
                           midpoint = 0, limit = c(-1, 1), space = "Lab",
                           name = "Correlation") +
      geom_text(aes(label = significance), color = "black") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(face = "italic"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold")) +
      labs(x = "Variables", y = taxa)
    
    if(returnMainFig){
      return(p)
    }else{
      ggsave(plot = p, filename = paste(dir,"/correlation_all_groups.png", sep = ""), dpi = 300, height = 6, width = 10, bg = 'white')
      # Save correlation dataframe with stats and correlation coefficients 
      write_xlsx(cor_melt, path = paste0(dir,"/all_groups_correlation.xlsx"))
    }
    
    
    if (showIndivCor) {
      
      # Create directory for individual scatterplots
      dir_indiv <- paste0(dir, "/scatterplots/")
      existingDirCheck(path = dir_indiv)
      
      # Iterate over each ASV and variable pair
      for (asv in asvList) {
        for (var in setdiff(colnames(df), asvList)) {
          
          #Prepare data for scatterplot
          scatter_data <- data_combined %>%
            dplyr::select(all_of(c(asv, var))) %>%
            dplyr::rename(counts_var = all_of(asv), variable = all_of(var))
          
          # Convert row names to column in scatter_data
          scatter_data$sample_id <- row.names(scatter_data)
          
          # Extract gg_group from sample_data(ps)
          gg_group_data <- data.frame(sample_id = row.names(sample_data(ps)), 
                                      group = sample_data(ps)[[varToCompare]])
          
          # Merge with scatter_data
          scatter_data <- merge(scatter_data, gg_group_data, by = "sample_id", all.x = TRUE)
          
          scatter_data$group <- factor(scatter_data$group, levels = levels(sample_data(ps)[[varToCompare]]))
          
          # Get taxa name (if taxonomic level of interest is not species)
          taxon <-  as.data.frame(tax_table(ps))[asv,taxa]
          
          # Skip if there are no data points
          if (nrow(scatter_data) == 0) next
          
          # Create scatterplot
          scatter_plot <- ggplot(scatter_data, aes(x = variable , y = counts_var, color = group)) +
            geom_point() +
            geom_smooth(method = "lm", color = "red", se = TRUE) +
            theme_minimal() +
            labs(
              title = paste("Scatterplot:", asv, "vs", var),
              x = "Relative Abundance (%)",
              y = var,
              color = "Group"
            )
          
          # Save scatterplot
          ggsave(
            filename = paste0(dir_indiv, "scatter_", ifelse(taxa == "Species", paste0(taxon, asv), taxon), "_vs_", var, ".png"),
            plot = scatter_plot,
            width = 8, height = 6, dpi = 300,
            bg = "white"
          )
          
          # Save scatter data
          write_xlsx(x = scatter_data, path = paste0(dir_indiv, "scatter_data_",ifelse(taxa == "Species", paste0(taxon, "_", asv), taxon), "_", var, ".xlsx"))
        }
      }
      message("Scatterplots saved to: ", dir_indiv)
    }
    
    
  }else{ #Produce multiple heatmaps for each group
    
    for(group in levels(sample_data(ps)[[gg_group]])){
      
      #Creates directory for groups at which we are looking
      dir_group <- paste(dir, "/", clean_string(group), sep = "")
      existingDirCheck(path = dir_group)
      
      # Set row names as sample_ids (Row.names)
      rownames(data_combined) <- data_combined$Row.names
      
      #Save list of samples linked to each group and creates prepare matrices for correlation
      samples <- sample_data(ps)[sample_data(ps)[[gg_group]] == group, "sample_id"]
      samples <- samples$sample_id
      data_grp <- data_combined[as.character(data_combined[,"Row.names"]) %in% samples, ]
      data_grp <- data_grp[,-1] #Remove first col which is Row.names because it is not a variable for the correlation
      asv_matrix <- as.matrix(data_grp[, asvList, drop = FALSE]) # ASVs
      var_matrix <- as.matrix(data_grp[, !colnames(data_grp) %in% asvList, drop = FALSE]) # Non-ASV variables
      
      # Calculate correlation matrix
      cor_results <- rcorr(as.matrix(data_grp), type = "spearman")
      cor_matrix <- cor_results$r  # Extract correlation coefficients
      p_values <- cor_results$P    # Extract p-values
      
      # print(cor_results)
      
      # Adjust p-values using FDR
      p_adjusted <- matrix(
        p.adjust(as.vector(p_values), method = "fdr"),
        nrow = nrow(p_values),
        ncol = ncol(p_values),
        dimnames = dimnames(p_values))
      
      cor_matrix <- cor_matrix[rownames(cor_matrix) %in% asvList, !colnames(cor_matrix) %in% asvList]
      p_values <- p_adjusted[rownames(p_adjusted) %in% asvList, !colnames(p_adjusted) %in% asvList]
      
      # Melt correlation and adjusted p-value matrices
      cor_melt <- melt(cor_matrix, varnames = c("ASV", "Variables"))
      p_melt <- melt(p_values, varnames = c("ASV", "Variables"))
      
      # Combine melted data and filter out NA correlations (if any)
      cor_melt$p_value <- p_melt$value
      cor_melt <- cor_melt[!is.na(cor_melt$value), ]  # Remove NA rows
      cor_melt$significance <- ifelse(cor_melt$p_value < 0.001, "***",
                                      ifelse(cor_melt$p_value < 0.01, "**",
                                             ifelse(cor_melt$p_value < 0.05, "*", "")))
      
      # Add taxa_name and taxa_name_asv variable (for species only)
      taxonomy <- as.data.frame(tax_table(ps))
      taxonomy <- taxonomy[asvList,]
      
      if(taxa=="Species"){
        cor_melt$taxa_name <- paste(taxonomy[cor_melt$ASV, "Genus"],taxonomy[cor_melt$ASV, "Species"])
        cor_melt$taxa_name <- paste(cor_melt$taxa_name, " (", cor_melt$ASV, ")", sep = "")
      }else{
        cor_melt$taxa_name <- taxonomy[cor_melt$ASV, taxa]
      }
      
      print(cor_melt)
      
      # Create the heatmap
      p <- ggplot(cor_melt, aes(x = Variables, y = taxa_name, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                             midpoint = 0, limit = c(-1, 1), space = "Lab",
                             name = "Correlation") +
        geom_text(aes(label = significance), color = "black") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_text(face = "italic"),
              axis.title.x = element_text(size = 12, face = "bold"),
              axis.title.y = element_text(size = 12, face = "bold")) +
        labs(x = "Variables", y = taxa)
      
      ggsave(plot = p, filename = paste0(dir_group,"/",clean_string(group),"_group_correlation.png"), dpi = 300, height = 6, width = 10, bg = 'white')
      
      # Save correlation dataframe with stats and correlation coefficients 
      write_xlsx(cor_melt, path = paste0(dir_group,"/",clean_string(group),"_group_correlation.xlsx"))
    }
    
  }
  
}
