log2foldChangeGraphMultipleGroups <- function(ps, deseq, gg_group, taxa = "Species", threshold = 0.05, FDR = TRUE,
                                              customColors, pairs, path, single_factor_design = FALSE,additionnalAes = NULL,
                                              dim = c(6,6)){
  
  #Creates directory for taxonomic level
  dir <- paste(path, taxa, sep = "")
  existingDirCheck(path = dir)
  
  counts <- t(otu_table(ps)) # Get counts from otu_table

  vs <- c() #Define empty list that will contain pairs comparaisons names
  
  #Save pairs comparaisons names for displaying them in the final sigtable
  for(k in seq_along(pairs)){
    
    #Save pair
    pair <- unlist(pairs[k])
    
    #Save specific pair comparaison
    vs <- c(vs, paste(clean_string(pair[1]), "vs", clean_string(pair[2]), sep="_"))
  }
  
  # If you performed deseq with a single variable containing the four groups
  if(single_factor_design){
    
    #Partition results for specific pairwise comparaisons
    res_subset1 <- results(deseq, contrast = list(c(resultsNames(deseq)[2]))) #wt putrescine vs vehicle
    sigtab_1 <- cbind(as(res_subset1, "data.frame"), as(tax_table(ps)[rownames(res_subset1), ], "matrix"))
    sigtab_1$comparaison <- 1
    sigtab_1$vs <- vs[1]
    
    res_subset2 <- results(deseq, contrast = list(c(resultsNames(deseq)[3], resultsNames(deseq)[4]))) #il22 ko putrescine vs vehicle
    sigtab_2 <- cbind(as(res_subset2, "data.frame"), as(tax_table(ps)[rownames(res_subset2), ], "matrix"))
    sigtab_2$comparaison <- 2
    sigtab_2$vs <- vs[2]
    
    res_subset3 <- results(deseq, contrast = list(c(resultsNames(deseq)[3]))) #vehicle wt vs il22 ko
    sigtab_3 <- cbind(as(res_subset3, "data.frame"), as(tax_table(ps)[rownames(res_subset3), ], "matrix"))
    sigtab_3$comparaison <- 3
    sigtab_3$vs <- vs[3]
    
    res_subset4 <- results(deseq, contrast = list(c(resultsNames(deseq)[2], resultsNames(deseq)[4]))) #putrescine wt vs il22 ko
    sigtab_4 <- cbind(as(res_subset4, "data.frame"), as(tax_table(ps)[rownames(res_subset4), ], "matrix"))
    sigtab_4$comparaison <- 4
    sigtab_4$vs <- vs[4]
    
  }else{ # If you performed deseq with 2 variables, each with 2 groups
    
    #Partition results for specific pairwise comparaisons
    res_subset1 <- results(deseq, contrast = list(c(resultsNames(deseq)[3]))) #wt putrescine vs vehicle
    sigtab_1 <- cbind(as(res_subset1, "data.frame"), as(tax_table(ps)[rownames(res_subset1), ], "matrix"))
    sigtab_1$comparaison <- 1
    sigtab_1$vs <- vs[1]
    
    res_subset2 <- results(deseq, contrast = list(c(resultsNames(deseq)[3], resultsNames(deseq)[4]))) #il22 ko putrescine vs vehicle
    sigtab_2 <- cbind(as(res_subset2, "data.frame"), as(tax_table(ps)[rownames(res_subset2), ], "matrix"))
    sigtab_2$comparaison <- 2
    sigtab_2$vs <- vs[2]
    
    res_subset3 <- results(deseq, contrast = list(c(resultsNames(deseq)[2]))) #vehicle wt vs il22 ko
    sigtab_3 <- cbind(as(res_subset3, "data.frame"), as(tax_table(ps)[rownames(res_subset3), ], "matrix"))
    sigtab_3$comparaison <- 3
    sigtab_3$vs <- vs[3]
    
    res_subset4 <- results(deseq, contrast = list(c(resultsNames(deseq)[2], resultsNames(deseq)[4]))) #putrescine wt vs il22 ko
    sigtab_4 <- cbind(as(res_subset4, "data.frame"), as(tax_table(ps)[rownames(res_subset4), ], "matrix"))
    sigtab_4$comparaison <- 4
    sigtab_4$vs <- vs[4]
    
  }
  
  #Append the sigtabs together
  sigtab <- bind_rows(sigtab_1, sigtab_2, sigtab_3, sigtab_4)
  
  #Add ASV variable col to sigtab (enables to store asv names, not only as rownames, because they will be changed when using rowbind)
  sigtab["asv"] <- gsub("\\..*", "", rownames(sigtab))
  
  #Keeping only ASVs for which they were taxa found at the taxonomical level of interest
  sigtab <- sigtab[!is.na(sigtab[[taxa]]),]
  
  if(FDR){
    pvalue <- as.character("padj")
  }else{
    pvalue <- as.character("pvalue")
  }
  
  #Replacing NA padj by 1 (they correspond to this anyways)
  sigtab$padj[is.na(sigtab[pvalue])] <- 1
  
  #Add column that adds symbols for the significance 
  # Define significance levels
  sigtab$significance <- as.character(cut(sigtab[,pvalue],
                                          breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                          labels = c("***", "**", "*", "NS")))
  
  sigtab <- sigtab[(sigtab[pvalue])<threshold,] # Reduce sigtab to only rows with significant value
  sigtab$Species <- paste(sigtab$Genus, " ", sigtab$Species, " (", sigtab$asv, ")", sep = "")
  sigtab$upper <- sigtab$log2FoldChange+sigtab$lfcSE
  sigtab$lower <- sigtab$log2FoldChange-sigtab$lfcSE
  sigtab$change_dir <- ifelse(sigtab$log2FoldChange < 0, "-", "+")
  
  for(cmp in unique(sigtab$vs)){
    
    sigtab_sub <- sigtab[sigtab$vs == cmp,] # Sub of sigtab for comparison of interest
    
    sigtab_sub <- sigtab_sub %>%
      arrange(log2FoldChange)
    sigtab_sub[[taxa]] <- factor(sigtab_sub[[taxa]], levels = rev(sigtab_sub[[taxa]]))
    
    max_val <- max(abs(c(sigtab_sub$lower, sigtab_sub$upper)))
    
    p <- ggplot(data = sigtab_sub, aes(x = log2FoldChange, y = .data[[taxa]], color = change_dir)) +
      geom_point(size = 3) +
      geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_color_manual(values = customColors) +
      scale_x_continuous(limits = c(-max_val, max_val), position = "top") +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.y = element_text(size = 10),
        legend.position = "right",
      ) +
      guides(color = "none")+
      labs(
        x = "Log Fold Change",
        y = NULL,
        color = ""
      )+annotate(
        "text",
        x = 0,
        y = -Inf,
        vjust = -1,      # Push below the panel edge
        hjust = 0.5,
        label = cmp
      )
    
    ggsave(filename = paste(path, cmp, "_lfc.png"), height = dim[1], width = dim[2], dpi = 300, bg = "white")
  }
}

log2foldChangeGraphSingleTimepoint <- function(ps, deseq, gg_group, timePoint, taxa = "Species", threshold = 0.05, FDR = TRUE,
                                              customColors, path, single_factor_design = FALSE,additionnalAes = NULL,
                                              dim = c(6,6)){
  
  #Creates directory for taxonomic level
  dir <- paste(path, taxa, sep = "")
  existingDirCheck(path = dir)
  
  #Save results at single timePoint
  res <- results(deseq, name = resultsNames(deseq)[2])
  
  #Save significance table
  sigtab <- cbind(as(res, "data.frame"), as(tax_table(ps)[rownames(res),], "matrix"))
  
  #Add ASV variable col to sigtab (enables to store asv names, not only as rownames, because they will be changed when using rowbind)
  sigtab["asv"] <- rownames(sigtab)
  
  #Keeping only ASVs for which they were taxa found at the taxonomical level of interest
  sigtab <- sigtab[!is.na(sigtab[[taxa]]),]
  
  if(FDR){
    pvalue <- as.character("padj")
  }else{
    pvalue <- as.character("pvalue")
  }
  
  #Replacing NA padj by 1 (they correspond to this anyways)
  sigtab[pvalue][is.na(sigtab[pvalue])] <- 1
  
  #Add column that adds symbols for the significance 
  # Define significance levels
  sigtab$significance <- as.character(cut(sigtab[,pvalue],
                                          breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                          labels = c("***", "**", "*", "NS")))
  
  sigtab <- sigtab[(sigtab[pvalue])<threshold,] # Reduce sigtab to only rows with significant value
  if(nrow(sigtab) == 0){ # Exit in case that there is nothing significant
    return(NULL)
  } 
  sigtab$Species <- paste(sigtab$Genus, " ", sigtab$Species, " (", sigtab$asv, ")", sep = "")
  sigtab$upper <- sigtab$log2FoldChange+sigtab$lfcSE
  sigtab$lower <- sigtab$log2FoldChange-sigtab$lfcSE
  sigtab$change_dir <- ifelse(sigtab$log2FoldChange < 0, "-", "+")
    
  sigtab <- sigtab %>% # Arrange taxa in increasing order of log2fc values
    arrange(log2FoldChange)
  sigtab[[taxa]] <- factor(sigtab[[taxa]], levels = rev(sigtab[[taxa]]))
    
  max_val <- max(abs(c(sigtab$lower, sigtab$upper)))
  
  p <- ggplot(data = sigtab, aes(x = log2FoldChange, y = .data[[taxa]], color = change_dir)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = customColors) +
    scale_x_continuous(limits = c(-max_val, max_val), position = "top") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.y = element_text(size = 10),
      legend.position = "right",
    ) +
    guides(color = "none")+
    labs(
      x = "Log Fold Change",
      y = NULL,
      color = ""
    )+annotate(
      "text",
      x = 0,
      y = -Inf,
      vjust = -1,      # Push below the panel edge
      hjust = 0.5,
      label = timePoint
    )
  
  ggsave(filename = paste(path, timePoint, "_lfc.png"), height = dim[1], width = dim[2], dpi = 300, bg = "white")
}