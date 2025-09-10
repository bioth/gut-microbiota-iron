#Produces the timeline with stats calculated with same fashion as for relabSingleTimepoint, works with design "~ factor + timeFactor + factor:timeFactor"
relabTimelineFullAnalysis <- function(ps, deseq, measure = "log2fold", timeVariable, varToCompare, taxa = "Species", threshold = 0.01, returnSigAsvs = FALSE, customColors, path){
  
  #Creates directory for taxonomic level
  dir <- paste(path, taxa, sep = "")
  existingDirCheck(path = dir)
  
  #Create empty sigtabCombined dataframe to merge sigtabs together
  sigtabCombined <- data.frame()
  
  #Define empty list that will contain pairs comparaisons names
  vs <- c()
  
  #Get counts from DESeq2 object
  counts <- t(otu_table(ps))
  
  #Save results at single timepoint
  res1 <- results(deseq, contrast = list(c(resultsNames(deseq)[2]))) # 50 vs 500 at week 3
  res2 <- results(deseq, contrast = list(c(resultsNames(deseq)[2],resultsNames(deseq)[6]))) # 50 vs 500 at week 8
  res3 <- results(deseq, contrast = list(c(resultsNames(deseq)[2],resultsNames(deseq)[7]))) # 50 vs 500 at week 10
  res4 <- results(deseq, contrast = list(c(resultsNames(deseq)[2],resultsNames(deseq)[8]))) # 50 vs 500 at week 14
  
  res5 <- results(deseq, contrast = list(c(resultsNames(deseq)[3]))) # week3 vs week8 for 50 ppm
  res6 <- results(deseq, contrast = list(c(resultsNames(deseq)[4]))) # week3 vs week10 for 50 ppm
  res7 <- results(deseq, contrast = list(c(resultsNames(deseq)[5]))) # week3 vs week14 for 50 ppm
  
  
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
    
  
  #Keeping only ASVs for which they were taxa found at the taxonomical level of interest
  sigtabCombined <- sigtabCombined[!is.na(sigtabCombined[[taxa]]),]
  
  #Find asvs that have at least one significant value at some timepoint
  asvList <- unique(sigtabCombined[(sigtabCombined$padj)<threshold,"asv"])
  
  #Loop along ASVs (single taxons analyzed)
  for(i in seq_along(asvList)){
    
    #Save asv value
    asv = asvList[i]
    
    #Save speciesName for dir creation and for graph title
    if(taxa == "Species"){
      taxonName <- paste(unique(sigtab[sigtab$asv == asv, "Genus"]),unique(sigtab[sigtab$asv == asv, "Species"]), paste("(", asv, ")", sep =""), sep = " ")
    }else{
      taxonName <- unique(sigtab[sigtab$asv == asv, taxa])
    }
    
    #Sigtab specific to the taxon analyzed
    sigtab_taxon <- sigtabCombined[sigtabCombined$asv == asv,]
    
    #Add the 16s sequence to the sigtab specific to ASV of interest (only for species, higher taxonomical levels have multiple ASVs associated with them)
    if(taxa == "Species"){
      sigtab_taxon$dna_sequence <- as.character(refseq(ps)[asv])
    }
    
    #Creates directory with taxon name
    dir_taxon <- paste(dir, "/", i, "-",taxonName, sep = "")
    existingDirCheck(path = dir_taxon)
    
    #Creates empty relative_abundance dataframe
    relative_abundance <- data.frame()
    
    #Iterate through timepoints and create an overall matrix of relative abundances
    for(j in seq_along(levels(sample_data(ps)[[timeVariable]]))){
      
      #Save timepoint value
      timePoint = levels(sample_data(ps)[[timeVariable]])[j]
      
      #Select counts for ASV of interest and at timepoint i
      nrm_counts <- nrm_counts_all[[j]]
      asv_counts <- nrm_counts[asv, ]
      
      #Calculate relative abundance as a percentage for each sample
      rel_ab <- asv_counts * 100 / colSums(nrm_counts)
      
      #Transform to dataframe and add timeVariable column
      rel_ab <- as.data.frame(rel_ab)
      
      #Append df to df with all
      relative_abundance <- rbind(relative_abundance, rel_ab)
      
    }
    
    #Add sample_id column for merge
    relative_abundance$sample_id <- rownames(relative_abundance)
    
    #Merge relative abundance with sample metadata
    relative_abundance <- merge(relative_abundance, as(sample_data(ps), "data.frame"), by = "sample_id")
    
    relative_abundance[[timeVariable]] <- as.numeric(as.character(relative_abundance[[timeVariable]])) * 7
    relative_abundance[[varToCompare]] <- as.character(relative_abundance[[varToCompare]])
    
    #Rearrange df in increasing order for timepoints
    relative_abundance <- relative_abundance %>%
      arrange(!!sym(varToCompare)) %>%
      arrange(!!sym(timeVariable))
    
    # Compute mean relative abundance by week and diet
    means_df <- relative_abundance %>%
      group_by(week, diet) %>%
      summarise(mean_abundance = mean(rel_ab)) %>%
      group_by(week) %>%
      summarise(upper_limit = max(mean_abundance), lower_limit = min(mean_abundance))
    
    #Before merging sigtab, ensure timeVariable is numeric and *7
    sigtab_taxon[[timeVariable]] <- as.numeric(as.character(sigtab_taxon[[timeVariable]])) * 7
    
    #Merge sigtab_taxon with means_df
    means_df <- merge(sigtab_taxon, means_df, by = timeVariable)
    
    p <- ggplot(data = relative_abundance, aes(x = !!sym(timeVariable), y = rel_ab, color = !!sym(varToCompare)), group = !!sym(varToCompare)) +
      
      geom_point(size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) + 
      
      # Error bars
      stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
                   aes(color = !!sym(varToCompare)),
                   width = 5, size = 1,
                   alpha = 0.5,
                   position = "identity")+
      
      #Mean lines
      stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
                   aes(ymin = ..y.., ymax = ..y.., group = !!sym(varToCompare)),
                   color = "black", linewidth = 0.5, width = 0.5,
                   position = "identity")+
      
      #Connecting mean points with lines
      stat_summary(fun = mean, geom = "line", size = 1.2) +  # Connecting means with lines
      
      
      
      labs(title = taxonName, y = "Relative abundance (%)", color = "Diet", x = "Time (days)") +
      scale_color_manual(values = customColors)+
      
      #Add vertical line segments for significance at each timepoint
      geom_segment(data = means_df, aes(x = !!sym(timeVariable)+1.6, xend = !!sym(timeVariable)+1.6,
                                        y = lower_limit, yend = upper_limit),
                   color = "black", linetype = "dashed", ) +
      
      #Complete with short horizontal segments to make the bar look nicer
      geom_segment(data = means_df, aes(x = !!sym(timeVariable)+1.6, xend = !!sym(timeVariable),
                                        y = upper_limit, yend = upper_limit),
                   color = "black", linetype = "dashed", ) +
      
      geom_segment(data = means_df, aes(x = !!sym(timeVariable)+1.6, xend = !!sym(timeVariable),
                                        y = lower_limit, yend = lower_limit),
                   color = "black", linetype = "dashed", ) +
      
      # Add significance text at each timepoint
      geom_text(data = means_df, aes(x = !!sym(timeVariable)+1.6, 
                                     y = (upper_limit + lower_limit) / 2, 
                                     label = significance),
                color = "black", size = 5, vjust = 0.5) +
      
      
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
    ggsave(plot = p, filename = paste(dir_taxon,"/",gsub(" ", "_", taxonName),"_relab.png", sep = ""), dpi = 300, height = 6, width = 6, bg = 'white')
    
    #Write as excel file the significance table specific to an ASV
    write.xlsx(sigtab_taxon, paste(dir_taxon,"/",gsub(" ", "_", taxonName),"_stats.xlsx", sep = ""))
    
    #Write as excel file the relative abundance data specific  to an ASV
    write.xlsx(relative_abundance, paste(dir_taxon,"/",gsub(" ", "_", taxonName),"_relab.xlsx", sep = ""))
    
    #Returns list of significant asvs
    if(returnSigAsvs){
      return(asvList)
    }
    
  }
}
