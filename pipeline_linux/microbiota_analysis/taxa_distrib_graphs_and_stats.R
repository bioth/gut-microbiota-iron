# Function that takes ps object as input, subset at a particular taxonomic level
# and returns an excel file with the statistics for comparaisons of interest
# and the actual percentage of relative abundance for each species
# You can add to that a visualization with pie charts, but it's optional

taxGlomResAndStats <- function(ps, taxrank, exp_group, selected_comparisons, include_graph, path){
  
  # Subset for taxon of interest
  ps_taxa <- tax_glom(ps, taxrank = taxrank)
  
  # Transform values in relative abundance
  ps_relab <- transformCounts(ps_taxa, transformation = "rel_ab")
  
  # Extract otu_table (rows as ASVs)
  relab_table <- as.data.frame(otu_table(ps_relab))
  
  # Create sub_tables for each exp group, extract means and combine everything into a new dataframe
  full_table <- data.frame(row.names = rownames(relab_table))
  for(group in levels(sample_data(ps)[[exp_group]])){
    samples_group <- rownames(sample_data(ps)[sample_data(ps)[,exp_group]==group,])
    sub_otu_table <- relab_table[,colnames(relab_table) %in% samples_group]
    full_table[[group]] <- rowMeans(sub_otu_table)
  }
  relab_table <- full_table
  
  # Add taxonomic information
  relab_table <- cbind(
    relab_table,
    tax_table(ps_taxa)[rownames(tax_table(ps_taxa)) %in% rownames(relab_table), drop = FALSE]
  )
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
  
  # Using correctly formatted contrasts based on results names
  results_list <- list()
  results_list <- lapply(comparisons, function(cmp) {
    contrast_vec <- c(exp_group, cmp[1], cmp[2])  # cmp[1] is the numerator, cmp[2] is the denominator
    tryCatch({
      
      res <- results(deseq, contrast = contrast_vec)
      
      # Link features with their taxonomy information
      res <- cbind(
        res,
        tax_table(ps_taxa)[rownames(tax_table(ps_taxa)) %in% rownames(res), drop = FALSE]
      )
      
      # Convert to data frame and add to the list
      res <- as.data.frame(res)
      return(res)
      
    }, error = function(e) {
      message("Failed to compute results for comparison: ", paste(cmp, collapse = " vs "), "\nError: ", e$message)
      return(NULL)  
    })
  })
  
  # Names for each table of the list
  names(results_list) <- sapply(comparisons, function(cmp) paste(cmp, collapse = "_vs_"))
  
  # Save combined significance table
  writeStackbarExtendedSigTable(results_list, filepath = paste0(path, "/", clean_string(taxrank), "_stats.xlsx"))
  
  # Generate pie charts for taxa distrubutions across groups
  if(include_graph){
    
    # Transforms data in long format
    data <- relab_table  %>% 
      pivot_longer(
        cols = 1:length(levels(sample_data(ps)[[exp_group]])), 
        names_to = exp_group,
        values_to = "relab"
      )
    data = as.data.frame(data) 
    data[[exp_group]] <- factor(data[[exp_group]], levels = c(levels(sample_data(ps)[[exp_group]])))  # Reorder the levels of exp_group
    
    # Plot
    p <- ggplot(data, aes(x = "", y = relab, fill = data[[taxrank]])) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      theme_void() +
      
      # Display percentages on the pie chart
      geom_text(
        aes(label = paste0(round(relab, 2), "%")), 
        position = position_stack(vjust = 0.5)  # Adjust position in the middle of the slices
      ) +
      facet_wrap(as.formula(paste("~", exp_group)))   # Create separate pie charts for each gg_group
    #scale_fill_brewer(palette = "Set3")  # Optional: Change color palette
    
    # Save plot
    ggsave(plot = p, filename = paste0(path, "/", clean_string(taxrank), "_plot.png"), bg = "white", width = 10, height = 5, dpi = 300)
  }
}


