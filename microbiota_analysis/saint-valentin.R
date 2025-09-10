# Function that takes ps object as input, subset at a particular taxonomic level
# and returns an excel file with the statistics for comparaisons of interest
# and the actual percentage of relative abundance for each species

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

taxGlomResAndStats <- function(ps, taxrank, exp_group, selected_comparaisons, path){
  
  # Subset for taxon of interest
  ps_taxa <- tax_glom(ps, taxrank = taxrank)
  
  # Transform values in relative abundance
  ps_relab <- transformCounts(ps_taxa, transformation = "rel_ab")
  
  # Save otu_table and save it 
  relab_table <- otu_table(ps_relab)
  existingDirCheck(path)
  write_xlsx(relab_table, path = paste0(path, "/", clean_string(taxrank), "_relative_abundance.xlsx"))
  
  # Convert to deseq2 object
  deseq <- phyloseq_to_deseq2(ps_taxa, formula(paste("~", exp_group)))
  
  # Run deseq analysis
  deseq <- DESeq(deseq, test="Wald", fitType = "parametric")
  
  # Save results and add taxonomical information to them
  results <- as.data.frame(cbind(results(deseq),
                                 tax_table(ps_taxa)[rownames(tax_table(ps_taxa)) %in% rownames(results(deseq))]))
  
  meta <- as.data.frame(sample_data(ps_taxa)) #sample metadata
  comparisons <- combn(levels(meta[[exp_group]]), 2, simplify = FALSE) # writes all possible comparaisons between yout experimental groups
  
  # Filter only selected comparisons
  comparisons <- Filter(function(cmp) {
    any(sapply(selected_comparisons, function(sel) all(sel == cmp)))
  }, comparisons)
  
  # Generate names for the results based on comparisons
  comparison_names <- sapply(comparisons, function(cmp) paste(cmp, collapse = "_vs_"))
  
  # Using correctly formatted contrasts based on results names
  results_list <- list()
  results_list <- lapply(comparisons, function(cmp) {
    contrast_vec <- c(exp_group, cmp[1], cmp[2])  # cmp[1] is the numerator, cmp[2] is the denominator
    tryCatch({
      res <- results(deseq, contrast = contrast_vec)
      return(res)
    }, error = function(e) {
      message("Failed to compute results for comparison: ", paste(cmp, collapse = " vs "), "\nError: ", e$message)
      return(NULL)  
    })
  })
  
  # Set names based on comparisons
  results_list <- setNames(results_list, comparison_names)
  
  res_final <- list()
  
  # Iterate over the results_list to process each comparison
  for (i in seq_along(results_list)) {
    res <- results_list[[i]]
    comparison_name <- names(results_list)[i]
    
    # Apply the significance threshold
    # significant_features <- subset(res, padj < fdr_threshold)
    
    
    # Link significant features with their taxonomy information
    res <- cbind(
      res,
      tax_table(ps_taxa)[rownames(tax_table(ps_taxa)) %in% rownames(res), drop = FALSE]
    )
    
    # Convert to data frame and add to the list
    res <- as.data.frame(res)
    res_final[[comparison_name]] <- res
    
  }
  
  writeStackbarExtendedSigTable(res_final, filepath = paste0(path, "/", taxrank, "_stats"))

}


