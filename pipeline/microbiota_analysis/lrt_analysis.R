# Function that identify bugs that are differentially abundant between timepoints and groups
relabTimepoints <- function(ps, deseq, measure = "log2fold", timeVariable, varToCompare, taxa = "Species", threshold = 0.01, customColors, path){
  
  # Creates directory for taxonomical level of interest
  dir <- paste0(path, taxa)
  existingDirCheck(dir)
  
  # Subset for taxonomical level of interest
  if(taxa != "Species"){
    ps <- tax_glom(ps, taxrank = taxa)
  }
  
  # Subsets for each group
  ps_sub1 <- prune_samples(sample_data(ps)$diet == "50", ps) # For 50 ppm 
  ps_sub2 <- prune_samples(sample_data(ps)$diet == "500", ps) # For 50 ppm 
  
  deseq_subset1 <- phyloseq_to_deseq2(ps_sub1, ~ week)
  deseq_subset1 <- DESeq(deseq_subset1, test="LRT", reduced = ~1)
  deseq_subset2 <- phyloseq_to_deseq2(ps_sub2, ~ week)
  deseq_subset2 <- DESeq(deseq_subset2, test="LRT", reduced = ~1)
  
  #Get counts
  counts <- t(otu_table(ps))
  
  #Save results for 3w vs 8w // 8w vs 10w
  res_1 <- results(deseq_subset1, contrast = c(timeVariable, "3","8"), test = "Wald")
  res_2 <- results(deseq_subset1, contrast = c(timeVariable, "8","10"), test = "Wald")
  res_3 <- results(deseq_subset2, contrast = c(timeVariable, "3","8"), test = "Wald")
  res_4 <- results(deseq_subset2, contrast = c(timeVariable, "8","10"), test = "Wald")
  
  #Save significance table
  sigtab_1 <- cbind(as(res_1, "data.frame"), as(tax_table(ps)[rownames(res_1),], "matrix"))
  sigtab_2 <- cbind(as(res_2, "data.frame"), as(tax_table(ps)[rownames(res_2),], "matrix"))
  sigtab_3 <- cbind(as(res_3, "data.frame"), as(tax_table(ps)[rownames(res_3),], "matrix"))
  sigtab_4 <- cbind(as(res_4, "data.frame"), as(tax_table(ps)[rownames(res_4),], "matrix"))
  
  # Add comparison column to the sigtab
  sigtab_1$comparison <- "3v8_50"
  sigtab_2$comparison <- "8v10_50"
  sigtab_3$comparison <- "3v8_500"
  sigtab_4$comparison <- "8v10_500"
  sigtab <- bind_rows(sigtab_1, sigtab_2, sigtab_3, sigtab_4)
  
  # #Create combined sigtab with stats at each timepoint
  # for(i in seq_along(resTimeVariable)){
  #   
  #   print(i)
  #   
  #   #Results subset for each timepoint
  #   res_subset <- results(deseq, name = resultsNames(deseq)[5+i])
  #   
  #   #Save significance table
  #   sigtab_subset <- cbind(as(res_subset, "data.frame"), as(tax_table(ps)[rownames(res_subset),], "matrix"))
  #   
  #   #Add timeVariable column to the sigtab
  #   sigtab_subset[[timeVariable]] <- levels(sample_data(ps)[[timeVariable]])[i+1]
  #   
  #   #Append the sigtabs together
  #   sigtab <- bind_rows(sigtab, sigtab_subset)
  #   
  # }
  
  # Add ASV variable col to sigtab (enables to store asv names, not only as rownames, because they will be changed when using rowbind)
  sigtab["asv"] <- gsub("\\..*", "", rownames(sigtab))
  
  # Keeping only ASVs for which they were taxa found at the taxonomical level of interest
  if(taxa == "Species"){
    sigtab <- sigtab[!is.na(sigtab[[taxa]]),]
  }
  
  #Replacing NA padj by 1 (they correspond to this anyways)
  sigtab$padj[is.na(sigtab$padj)] <- 1
  
  #Add column that adds symbols for the significance 
  # Define significance levels
  sigtab$significance <- cut(sigtab$padj,
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", "NS"))
  
  # Find asvs that have at least one significant value at some timepoint
  asvList <- unique(sigtab[(sigtab$padj)<threshold,"asv"])
  
  #Loop along ASVs (single taxons analyzed)
  for(i in seq_along(asvList)){
    
    #Save asv value
    asv = asvList[i]
    
    #Save speciesName for dir creation and for graph title
    if(taxa == "Species"){
      taxonName <- paste(unique(sigtab[sigtab$asv == asv, "Genus"]),unique(sigtab[sigtab$asv == asv, "Species"]), sep = " ")
    }else{
      taxonName <- unique(sigtab[sigtab$asv == asv, taxa])
    }
    
    #Sigtab specific to the taxon analyzed
    sigtab_taxon <- sigtab[sigtab$asv == asv,]
    
    #Creates directory with taxon name
    dir_taxon <- paste0(dir, "/", i, "-",taxonName)
    existingDirCheck(path = dir_taxon)
    
    # For each timepoint, calculate relative abundance data 
    relative_abundance <- data.frame()
    for(t in levels(sample_data(ps)[[timeVariable]])){
      
      #Select counts for ASV of interest and at tp of interest
      tp_counts <- counts[,sample_data(ps)[sample_data(ps)[[timeVariable]] == t,]$sample_id]
      asv_counts <- tp_counts[asv, ]
      df <- asv_counts * 100 / colSums(tp_counts)
      df <- data.frame(
        sample_id = colnames(df),
        rel_ab = as.numeric(df)
      )
      df$timepoint <- t
      relative_abundance <- rbind(relative_abundance, df)
    }
    
    #Merge relative abundance with sample metadata
    relative_abundance <- merge(relative_abundance, as(sample_data(ps), "data.frame"), by = "sample_id")
    relative_abundance[[timeVariable]] <- as.numeric(as.character(relative_abundance[[timeVariable]])) * 7
    relative_abundance[[varToCompare]] <- as.character(relative_abundance[[varToCompare]])
    
    # # Compute mean relative abundance by week and diet
    # means_df <- relative_abundance %>%
    #   group_by(week, diet) %>%
    #   summarise(mean_abundance = mean(rel_ab)) %>%
    #   group_by(week) %>%
    #   summarise(upper_limit = max(mean_abundance), lower_limit = min(mean_abundance))
    
    
    # #Before merging sigtab, ensure timeVariable is numeric and *7
    # sigtab_taxon[[timeVariable]] <- as.numeric(as.character(sigtab_taxon[[timeVariable]])) * 7
    
    # #Merge sigtab_taxon with means_df
    # means_df <- merge(sigtab_taxon, means_df, by = timeVariable)
    # 
    # print(means_df)
    
    
    p <- ggplot(data = relative_abundance, aes(x = !!sym(timeVariable), y = rel_ab, color = !!sym(varToCompare)), group = !!sym(varToCompare)) +
      
      geom_point(size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) + 
      
      # Error bars
      stat_summary(fun.data = "mean_se", geom = "errorbar",
                   aes(color = !!sym(varToCompare)),
                   width = 5, size = 1,
                   alpha = 0.5,
                   position = "identity")+
      
      #Mean lines
      stat_summary(fun.data = "mean_se", geom = "errorbar",
                   aes(ymin = ..y.., ymax = ..y.., group = !!sym(varToCompare)),
                   color = "black", linewidth = 0.5, width = 0.5,
                   position = "identity")+
      
      #Connecting mean points with lines
      stat_summary(fun = mean, geom = "line", size = 1.2) +  # Connecting means with lines
      
      
      
      labs(title = taxonName, y = "Relative abundance (%)", color = "Diet", x = "Time (days)") +
      scale_color_manual(values = customColors)+
      
      # #Add vertical line segments for significance at each timepoint
      # geom_segment(data = means_df, aes(x = !!sym(timeVariable)+1.6, xend = !!sym(timeVariable)+1.6,
      #                                   y = lower_limit, yend = upper_limit),
      #              color = "black", linetype = "dashed", ) +
      # 
      # #Complete with short horizontal segments to make the bar look nicer
      # geom_segment(data = means_df, aes(x = !!sym(timeVariable)+1.6, xend = !!sym(timeVariable),
      #                                   y = upper_limit, yend = upper_limit),
      #              color = "black", linetype = "dashed", ) +
      # 
      # geom_segment(data = means_df, aes(x = !!sym(timeVariable)+1.6, xend = !!sym(timeVariable),
      #                                   y = lower_limit, yend = lower_limit),
      #              color = "black", linetype = "dashed", ) +
      # 
      # # Add significance text at each timepoint
      # geom_text(data = means_df, aes(x = !!sym(timeVariable)+1.6, 
      #                                y = (upper_limit + lower_limit) / 2, 
      #                                label = significance),
      #           color = "black", size = 5, vjust = 0.5) +
      
      
     my_theme()
    ggsave(plot = p, filename = paste(dir_taxon,"/",gsub(" ", "_", taxonName),"_lrt.png", sep = ""), dpi = 300, height = 6, width = 6, bg = 'white')
    write.xlsx(sigtab_taxon, paste(dir_taxon,"/",gsub(" ", "_", taxonName),"_stats.xlsx", sep = ""))
    
  }
}