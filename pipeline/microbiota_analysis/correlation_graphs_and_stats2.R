# Make Spearman correlation between variables from a dataframe and counts of taxa abundance 
# for selected taxa (as ASV number), ps otu_table should have same rownames as the variables dataframe

correlationSelectedASVs <-  function(df_var, ps, selected_ASVs, taxrank = "Species",
                        timepoint = "T112"){
  
  if(taxa_are_rows(ps)){
    counts <- t(otu_table(ps))[,c(selected_ASVs)]
  }else{
    counts <- otu_table(ps)[,c(selected_ASVs)]
  } # Get asv counts of taxa of interest
  
  row.names(counts) <- str_replace(row.names(counts), pattern = paste0("_",timepoint), replacement = "")
  df_cor<- merge(df_var, counts, by = "row.names")
  rownames(df_cor) <- df_cor$Row.names
  df_cor <- df_cor[,-1]
  
  cor_results <- rcorr(as.matrix(df_cor), type = "spearman")
  cor_matrix <- cor_results$r  # Extract correlation coefficients
  # cor_matrix <- as.data.frame(cor_matrix[(1+ncol(df_var)):(ncol(df_var)+length(selected_ASVs)),
  #                                        1:(ncol(df_var))])
  
  cor_matrix <- cor_matrix[1:(ncol(df_var)),
                                         (1+ncol(df_var)):(ncol(df_var)+length(selected_ASVs))]
  
  p_values <- cor_results$P    # Extract p-values
  p_values <- p_values[1:(ncol(df_var)),
                                     (1+ncol(df_var)):(ncol(df_var)+length(selected_ASVs))]
  
  # FDR correction
  p_adjusted <- matrix(
    p.adjust(as.vector(p_values), method = "fdr"),
    nrow = nrow(p_values),
    ncol = ncol(p_values),
    dimnames = dimnames(p_values))
  print(p_adjusted)
  
  # Preparing dataframe to plot
  cor_melt <- melt(cor_matrix, varnames = c("Variables", "Taxa"))
  p_melt <- melt(p_adjusted, varnames = c("Variables","Taxa"))
  cor_melt$p_value <- p_melt$value # Combine melted data
  cor_melt <- cor_melt[!is.na(cor_melt$value), ]  # Remove NA rows
  cor_melt$significance <- ifelse(cor_melt$p_value < 0.001, "***", 
                                  ifelse(cor_melt$p_value < 0.01, "**", 
                                         ifelse(cor_melt$p_value < 0.05, "*", "")))
  
  # Prepare taxa names to be displayed
  taxaNames <- paste0(as.character(tax_table(ps)[selected_ASVs,taxrank]),"\n(", timepoint, ")")
  
  # Plot correlation heatmap graph
  p <- ggplot(cor_melt, aes(x = Variables, y = Taxa, fill = value)) + 
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#d1762d", high = "#639381", mid = "white",
                         midpoint = 0, limit = c(-1, 1), space = "Lab",
                         name = "Spearman's\nCorrelation") +
    geom_text(aes(label = significance), color = "black") +
    labs(y = taxrank, x = "Variables")+
    scale_y_discrete(labels  = taxaNames)+
    my_theme()+
    theme(axis.line = element_blank(),
          legend.background = element_blank())
  
  return(p)
  
  
}