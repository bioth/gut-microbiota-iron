
relabSpeciesPairwise <- function(ps, deseq, measure = "log2fold", gg_group, pairs, threshold = 0.01, customColors, path){

  #Get normalized counts from DESeq2 object
  normalized_counts <- counts(deseq, normalized = TRUE)
  
  for(i in seq_along(pairs)){
    
    #Save pair
    pair <- unlist(pairs[i])
    
    #Save color pair
    colorPair <- unlist(customColors[i])
    
    #Partition results for a specific pair
    res <- results(deseq, contrast=c(gg_group, pair))
    
    #Save significance table
    sigtab <- cbind(as(res, "data.frame"), as(tax_table(ps)[rownames(res), ], "matrix"))
    
    #Keeping only ASVs for which species were identified 
    sigtab <- sigtab[!is.na(sigtab$Species),]
    
    #Replacing NA padj by 1 (they correspond to this anyways)
    sigtab$padj[is.na(sigtab$padj)] <- 1
    
    #Keeping only padj < threshold (default = 0.01)
    sigtab <- sigtab[(sigtab$padj)<threshold,]
    
    #Save list of ASVs that are significantly different between pair groups
    asvList <- rownames(sigtab)
    
    # results()$log2FoldChange
    
    #define name for folders 
    vs <- paste(clean_string(pair[1]), "vs", clean_string(pair[2]), sep="_")
    
    #Define dir path where graphs and stats are gonna be saved
    dir <- paste(path, vs, sep = "")
    
    #Checking if dir already exists, otherwise creates it
    existingDirCheck(path = dir)
    
    for(asv in asvList){
      
      #Save speciesName for dir creation and for graph title
      speciesName <- paste(sigtab[asv,"Genus"], sigtab[asv,"Species"], sep = " ")
      
      #Creates directory with species name of ASV
      dir_species <- paste(dir, "/", gsub(" ", "_", speciesName), sep = "")
      existingDirCheck(path = dir_species)
      
      #Select counts for ASV of interest
      asv_counts <- normalized_counts[asv, ]
      
      #Calculate relative abundance as a percentage for each sample
      relative_abundance <- asv_counts * 100 / colSums(normalized_counts)
      
      #Remove sample IDs that are not part of the gg_group pair
      relative_abundance <- relative_abundance[sample_data(ps)$sample_id[sample_data(ps)[[gg_group]] %in% pair]]
      
      #Convert relative abundance to a data frame
      relative_abundance <- data.frame(
        sample_id = names(relative_abundance),
        rel_ab = as.numeric(relative_abundance)
      )
      
      # Merge relative abundance with sample metadata
      relative_abundance <- merge(relative_abundance, as(sample_data(ps), "data.frame"), by = "sample_id")
      
      #Save associated p-value
      p_value <- sigtab[asv, "padj"]

      p <- ggplot(data = relative_abundance, aes(x = gg_group, y = rel_ab, color = gg_group)) +
        geom_point(size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) + 
        
        # Error bars
        stat_summary(fun.data = "mean_se", geom = "errorbar",
                     aes(color = gg_group),
                     width = 0.2, size = 0.7,
                     position = position_dodge(-0.75)) +
        
        #Mean lines
        stat_summary(fun.data = "mean_se", geom = "errorbar",
                     aes(ymin = ..y.., ymax = ..y.., group = gg_group),
                     color = "black", linewidth = 0.5, width = 0.5,
                     position = position_dodge(-0.75))+
        
        
        labs(title = speciesName,
             y = "Relative abundance (%)", color = "Groups", x = "Groups") +
        scale_color_manual(values = colorPair)+
        
        
        #Add significance bar
        geom_signif(comparisons = list(pair),
                    annotations = paste0("p = ", format(p_value, digits = 2, scientific = TRUE)),  # Display p-value
                    y_position = max(relative_abundance$rel_ab, na.rm = TRUE) + 0.1,
                    tip_length = 0.02,
                    size = 1.2,  # Make the bar wider
                    color = "black") +

        theme_minimal()+
        theme(
          plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
          axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust x-axis tick label font size
          axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
          legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
          legend.text = element_text(size = 12),  # Adjust legend font size
          panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bar
      ggsave(plot = p, filename = paste(dir_species,"/",measure,"_",vs,".png", sep = ""), dpi = 300, height = 6, width = 6, bg = 'white')
      
    }
    

    
    # #Ordination plot
    # p <- plot_ordination(ps_subset, pcoa_results, type = "samples", 
    #                      color = gg_group) + 
    #   theme_classic() +
    #   theme(strip.background = element_blank())+
    #   stat_ellipse(aes(group = !!sym(gg_group)),      # Add ellipses grouping points by genotype
    #                type = "t",  # t-distribution for better fit
    #                level = 0.95,  # Confidence level for the ellipse                     
    #                geom = "polygon", alpha = 0)+
    #   labs(title = paste("PCoA of", distCharacter, "distance matrix.", sep = " ")) +
    #   scale_color_manual(values = colorPair)+
    #   labs(color = "Groups")+
    #   
    #   theme(
    #     plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    #     axis.title.x = element_text(size = 12),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    #     axis.text.x = element_text(size = 12),  # Adjust x-axis tick label font size
    #     axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
    #     legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    #     legend.text = element_text(size = 12),  # Adjust legend font size
    #     panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
    #     panel.grid.minor = element_blank(),  # Remove minor grid lines
    #     axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
    # 
    # 
    # test.adonis <- adonis(as.formula(paste("dist_subset ~", gg_group)), data = data.frame(sample_data(ps_subset)))
    # test.adonis <- as.data.frame(test.adonis$aov.tab)
    # print(test.adonis)
    # write.table(test.adonis, file = paste(dir,"/",measure,"_",vs,".tsv", sep = ""), col.names = NA, row.names = TRUE, sep = '\t', quote = FALSE)
    # ggsave(plot = p, filename = paste(dir,"/",measure,"_",vs,".png", sep = ""), dpi = 300, height = 8, width = 8, bg = 'white')
    
    
    # cbn <- combn(x=unique(metadata$body.site), m = 2)
    # p <- c()
    # 
    # for(i in 1:ncol(cbn)){
    #   ps.subs <- subset_samples(ps.rarefied, body.site %in% cbn[,i])
    #   metadata_sub <- data.frame(sample_data(ps.subs))
    #   permanova_pairwise <- adonis(phyloseq::distance(ps.subs, method = "bray") ~ body.site, 
    #                                data = metadata_sub)
    #   p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
    # }
    # 
    # p.adj <- p.adjust(p, method = "BH")
    # p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
    # p.table
    
  }
}

relabSpeciesTimePoint <- function(ps, deseq, measure = "log2fold", timeVariable, varToCompare, threshold = 0.01, customColors, path){
  
  #Get normalized counts from DESeq2 object
  normalized_counts <- counts(deseq, normalized = TRUE)
  
  for(timepoint in levels(sample_data(ps)[[timeVariable]])){
    
    #Partition results for a specific date
    res <- results(deseq, contrast=c(timeVariable, timepoint))
    
    #Save significance table
    sigtab <- cbind(as(res, "data.frame"), as(tax_table(ps)[rownames(res), ], "matrix"))
    
    #Keeping only ASVs for which species were identified 
    sigtab <- sigtab[!is.na(sigtab$Species),]
    
    #Replacing NA padj by 1 (they correspond to this anyways)
    sigtab$padj[is.na(sigtab$padj)] <- 1
    
    #Keeping only padj < threshold (default = 0.01)
    sigtab <- sigtab[(sigtab$padj)<threshold,]
    
    #Save list of ASVs that are significantly different between pair groups
    asvList <- rownames(sigtab)
    
    # results()$log2FoldChange
    
    #Define dir path where graphs and stats are gonna be saved
    dir <- paste(path,"week_",as.character(timepoint), sep = "")
    
    #Checking if dir already exists, otherwise creates it
    existingDirCheck(path = dir)
    
    for(asv in asvList){
      
      #Save speciesName for dir creation and for graph title
      speciesName <- paste(sigtab[asv,"Genus"], sigtab[asv,"Species"], sep = " ")
      
      #Creates directory with species name of ASV
      dir_species <- paste(dir, "/", gsub(" ", "_", speciesName), sep = "")
      existingDirCheck(path = dir_species)
      
      #Select counts for ASV of interest
      asv_counts <- normalized_counts[asv, ]
      
      #Calculate relative abundance as a percentage for each sample
      relative_abundance <- asv_counts * 100 / colSums(normalized_counts)
      
      #Remove sample IDs that are not part of the gg_group pair
      relative_abundance <- relative_abundance[sample_data(ps)$sample_id[sample_data(ps)[[gg_group]] %in% pair]]
      
      #Convert relative abundance to a data frame
      relative_abundance <- data.frame(
        sample_id = names(relative_abundance),
        rel_ab = as.numeric(relative_abundance)
      )
      
      # Merge relative abundance with sample metadata
      relative_abundance <- merge(relative_abundance, as(sample_data(ps), "data.frame"), by = "sample_id")
      
      #Save associated p-value
      p_value <- sigtab[asv, "padj"]
      
      p <- ggplot(data = relative_abundance, aes(x = gg_group, y = rel_ab, color = gg_group)) +
        geom_point(size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) + 
        
        # Error bars
        stat_summary(fun.data = "mean_se", geom = "errorbar",
                     aes(color = gg_group),
                     width = 0.2, size = 0.7,
                     position = position_dodge(-0.75)) +
        
        #Mean lines
        stat_summary(fun.data = "mean_se", geom = "errorbar",
                     aes(ymin = ..y.., ymax = ..y.., group = gg_group),
                     color = "black", linewidth = 0.5, width = 0.5,
                     position = position_dodge(-0.75))+
        
        
        labs(title = speciesName,
             y = "Relative abundance (%)", color = "Groups", x = "Groups") +
        scale_color_manual(values = colorPair)+
        
        
        #Add significance bar
        geom_signif(comparisons = list(pair),
                    annotations = paste0("p = ", format(p_value, digits = 2, scientific = TRUE)),  # Display p-value
                    y_position = max(relative_abundance$rel_ab, na.rm = TRUE) + 0.1,
                    tip_length = 0.02,
                    size = 1.2,  # Make the bar wider
                    color = "black") +
        
        
        theme(
          plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
          axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust x-axis tick label font size
          axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
          legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
          legend.text = element_text(size = 12),  # Adjust legend font size
          panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bar
      ggsave(plot = p, filename = paste(dir_species,"/",measure,"_",vs,".png", sep = ""), dpi = 300, height = 8, width = 8, bg = 'white')
      
    }
    
    
    
    # #Ordination plot
    # p <- plot_ordination(ps_subset, pcoa_results, type = "samples", 
    #                      color = gg_group) + 
    #   theme_classic() +
    #   theme(strip.background = element_blank())+
    #   stat_ellipse(aes(group = !!sym(gg_group)),      # Add ellipses grouping points by genotype
    #                type = "t",  # t-distribution for better fit
    #                level = 0.95,  # Confidence level for the ellipse                     
    #                geom = "polygon", alpha = 0)+
    #   labs(title = paste("PCoA of", distCharacter, "distance matrix.", sep = " ")) +
    #   scale_color_manual(values = colorPair)+
    #   labs(color = "Groups")+
    #   
    #   theme(
    #     plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    #     axis.title.x = element_text(size = 12),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    #     axis.text.x = element_text(size = 12),  # Adjust x-axis tick label font size
    #     axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
    #     legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    #     legend.text = element_text(size = 12),  # Adjust legend font size
    #     panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
    #     panel.grid.minor = element_blank(),  # Remove minor grid lines
    #     axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
    # 
    # 
    # test.adonis <- adonis(as.formula(paste("dist_subset ~", gg_group)), data = data.frame(sample_data(ps_subset)))
    # test.adonis <- as.data.frame(test.adonis$aov.tab)
    # print(test.adonis)
    # write.table(test.adonis, file = paste(dir,"/",measure,"_",vs,".tsv", sep = ""), col.names = NA, row.names = TRUE, sep = '\t', quote = FALSE)
    # ggsave(plot = p, filename = paste(dir,"/",measure,"_",vs,".png", sep = ""), dpi = 300, height = 8, width = 8, bg = 'white')
    
    
    # cbn <- combn(x=unique(metadata$body.site), m = 2)
    # p <- c()
    # 
    # for(i in 1:ncol(cbn)){
    #   ps.subs <- subset_samples(ps.rarefied, body.site %in% cbn[,i])
    #   metadata_sub <- data.frame(sample_data(ps.subs))
    #   permanova_pairwise <- adonis(phyloseq::distance(ps.subs, method = "bray") ~ body.site, 
    #                                data = metadata_sub)
    #   p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
    # }
    # 
    # p.adj <- p.adjust(p, method = "BH")
    # p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
    # p.table
    
  }
}

relabSpeciesTimeline <- function(ps, deseq, measure = "log2fold", timeVariable, varToCompare, threshold = 0.01, customColors, path){
  
  #Get normalized counts from DESeq2 object
  normalized_counts <- counts(deseq, normalized = TRUE)
    
    #Save results
    res <- results(deseq)
    
    #Save significance table
    sigtab <- cbind(as(res, "data.frame"), as(tax_table(ps)[rownames(res), ], "matrix"))
    
    #Keeping only ASVs for which species were identified 
    sigtab <- sigtab[!is.na(sigtab$Species),]
    
    #Replacing NA padj by 1 (they correspond to this anyways)
    sigtab$padj[is.na(sigtab$padj)] <- 1
    
    #Keeping only padj < threshold (default = 0.01)
    sigtab <- sigtab[(sigtab$padj)<threshold,]
    
    #Save list of ASVs that are significantly different
    asvList <- rownames(sigtab)
    
    #Write statistics into a txt file in the path folder
    write_xlsx(sigtab, paste(path,"/statistics_species.xlsx", sep =""), col_names = TRUE)
    
    # results()$log2FoldChange
    
    #Using seq along enables to add index in case same species was found for an ASV
    for(i in seq_along(asvList)){
      
      
      asv = asvList[i]
      
      #Save speciesName for dir creation and for graph title
      speciesName <- paste(sigtab[asv,"Genus"], sigtab[asv,"Species"], sep = " ")
      
      #Creates directory with species name of ASV
      dir_species <- paste(path, i, "-", gsub(" ", "_", speciesName), sep = "")
      existingDirCheck(path = dir_species)
      
      #Select counts for ASV of interest
      asv_counts <- normalized_counts[asv, ]
      
      #Calculate relative abundance as a percentage for each sample
      relative_abundance <- asv_counts * 100 / colSums(normalized_counts)
      
      #Convert relative abundance to a data frame
      relative_abundance <- data.frame(
        sample_id = names(relative_abundance),
        rel_ab = as.numeric(relative_abundance)
      )
      
      # Merge relative abundance with sample metadata
      relative_abundance <- merge(relative_abundance, as(sample_data(ps), "data.frame"), by = "sample_id")
      
      relative_abundance[[timeVariable]] <- as.numeric(as.character(relative_abundance[[timeVariable]])) * 7
      relative_abundance[[varToCompare]] <- as.character(relative_abundance[[varToCompare]])
      
      #Save associated p-value
      p_value <- sigtab[asv, "padj"]
      
      p <- ggplot(data = relative_abundance, aes(x = !!sym(timeVariable), y = rel_ab, color = !!sym(varToCompare)), group = !!sym(varToCompare)) +
        
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
        
     
        
        labs(title = speciesName,
             y = "Relative abundance (%)", color = "Diet", x = "Time (days)") +
        scale_color_manual(values = customColors)+
        # scale_x_continuous(limits = c(0, NA))+
        
        
        #Add significance bar
        # geom_signif(comparisons = list(pair),
        #             annotations = paste0("p = ", format(p_value, digits = 2, scientific = TRUE)),  # Display p-value
        #             y_position = max(relative_abundance$rel_ab, na.rm = TRUE) + 0.1,
        #             tip_length = 0.02,
        #             size = 1.2,  # Make the bar wider
        #             color = "black") +
        
        theme_minimal()+
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
          axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bar
      ggsave(plot = p, filename = paste(dir_species,"/",gsub(" ", "_", speciesName),"_relab.png", sep = ""), dpi = 300, height = 6, width = 6, bg = 'white')
      
    }
    
    
    
    # #Ordination plot
    # p <- plot_ordination(ps_subset, pcoa_results, type = "samples", 
    #                      color = gg_group) + 
    #   theme_classic() +
    #   theme(strip.background = element_blank())+
    #   stat_ellipse(aes(group = !!sym(gg_group)),      # Add ellipses grouping points by genotype
    #                type = "t",  # t-distribution for better fit
    #                level = 0.95,  # Confidence level for the ellipse                     
    #                geom = "polygon", alpha = 0)+
    #   labs(title = paste("PCoA of", distCharacter, "distance matrix.", sep = " ")) +
    #   scale_color_manual(values = colorPair)+
    #   labs(color = "Groups")+
    #   
    #   theme(
    #     plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    #     axis.title.x = element_text(size = 12),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    #     axis.text.x = element_text(size = 12),  # Adjust x-axis tick label font size
    #     axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
    #     legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    #     legend.text = element_text(size = 12),  # Adjust legend font size
    #     panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
    #     panel.grid.minor = element_blank(),  # Remove minor grid lines
    #     axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
    # 
    # 
    # test.adonis <- adonis(as.formula(paste("dist_subset ~", gg_group)), data = data.frame(sample_data(ps_subset)))
    # test.adonis <- as.data.frame(test.adonis$aov.tab)
    # print(test.adonis)
    # write.table(test.adonis, file = paste(dir,"/",measure,"_",vs,".tsv", sep = ""), col.names = NA, row.names = TRUE, sep = '\t', quote = FALSE)
    # ggsave(plot = p, filename = paste(dir,"/",measure,"_",vs,".png", sep = ""), dpi = 300, height = 8, width = 8, bg = 'white')
    
    
    # cbn <- combn(x=unique(metadata$body.site), m = 2)
    # p <- c()
    # 
    # for(i in 1:ncol(cbn)){
    #   ps.subs <- subset_samples(ps.rarefied, body.site %in% cbn[,i])
    #   metadata_sub <- data.frame(sample_data(ps.subs))
    #   permanova_pairwise <- adonis(phyloseq::distance(ps.subs, method = "bray") ~ body.site, 
    #                                data = metadata_sub)
    #   p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
    # }
    # 
    # p.adj <- p.adjust(p, method = "BH")
    # p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
    # p.table
    
}

#Can creates relative abundance graph for a single varToCompare and for several time points, works at any taxonomic level
relabTimeline <- function(ps, deseq, measure = "log2fold", timeVariable, varToCompare, taxa = "Species", threshold = 0.01, customColors, path){
  
  #Creates directory for taxonomic level
  dir <- paste(path, taxa, sep = "")
  existingDirCheck(path = dir)
  
  #Get normalized counts from DESeq2 object
  normalized_counts <- counts(deseq, normalized = TRUE)
  
  #Variable similar to timeVariable but without the first term (because first time point is reference time point)
  resTimeVariable <- levels(sample_data(ps)[[timeVariable]])[-1]
  
  resultsNames(deseq)
  
  #Save results at the first timepoint
  res <- results(deseq, contrast = list(resultsNames(deseq)[5]),resultsNames(deseq)[1])
  
  #Save significance table
  sigtab <- cbind(as(res, "data.frame"), as(tax_table(ps)[rownames(res),], "matrix"))
  
  #Add timeVariable column to the sigtab
  sigtab[[timeVariable]] <- levels(sample_data(ps)[[timeVariable]])[1]
  
  print(sigtab)
  
  #Create combined sigtab with stats at each timepoint
  for(i in seq_along(resTimeVariable)){
    
    print(i)
    
    #Results subset for each timepoint
    res_subset <- results(deseq, name = resultsNames(deseq)[5+i])
    
    #Save significance table
    sigtab_subset <- cbind(as(res_subset, "data.frame"), as(tax_table(ps)[rownames(res_subset),], "matrix"))
    
    #Add timeVariable column to the sigtab
    sigtab_subset[[timeVariable]] <- levels(sample_data(ps)[[timeVariable]])[i+1]
    
    #Append the sigtabs together
    sigtab <- bind_rows(sigtab, sigtab_subset)
    
  }
  
  #Add ASV variable col to sigtab (enables to store asv names, not only as rownames, because they will be changed when using rowbind)
  sigtab["asv"] <- gsub("\\..*", "", rownames(sigtab))

  #Keeping only ASVs for which they were taxa found at the taxonomical level of interest
  sigtab <- sigtab[!is.na(sigtab[[taxa]]),]
  
  #Replacing NA padj by 1 (they correspond to this anyways)
  sigtab$padj[is.na(sigtab$padj)] <- 1
  
  #Add column that adds symbols for the significance 
  # Define significance levels
  sigtab$significance <- cut(sigtab$padj,
                                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                   labels = c("***", "**", "*", "NS"))

  #Find asvs that have at least one significant value at some timepoint
  asvList <- unique(sigtab[(sigtab$padj)<threshold,"asv"])

  # results()$log2FoldChange
  
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
    dir_taxon <- paste(dir, "/", i, "-",taxonName, sep = "")
    existingDirCheck(path = dir_taxon)
    
    #Select counts for ASV of interest
    asv_counts <- normalized_counts[asv, ]
    
    #Calculate relative abundance as a percentage for each sample
    relative_abundance <- asv_counts * 100 / colSums(normalized_counts)
    
    #Convert relative abundance to a data frame
    relative_abundance <- data.frame(
      sample_id = names(relative_abundance),
      rel_ab = as.numeric(relative_abundance)
    )
    
    #Merge relative abundance with sample metadata
    relative_abundance <- merge(relative_abundance, as(sample_data(ps), "data.frame"), by = "sample_id")
    
    relative_abundance[[timeVariable]] <- as.numeric(as.character(relative_abundance[[timeVariable]])) * 7
    relative_abundance[[varToCompare]] <- as.character(relative_abundance[[varToCompare]])

    
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
    
    print(means_df)

    
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
    
  }
}

#For design with 4 groups based on 2 conditions - this the latest version used
# gg_group must be order with correct order prior to that (as a factor)
relabGroups <- function(ps, deseq, measure = "log2fold", gg_group, taxa = "Species", threshold = 0.01, FDR = TRUE,
                        returnSigAsvs = FALSE, normalizeCounts = FALSE, customColors, pairs, path, single_factor_design = FALSE,
                        additionnalAes = NULL, dim = c(6,6), displayPvalue = TRUE, displaySignificance = TRUE, includeUnknownSpecies = FALSE){
  
  #Creates directory for taxonomic level
  dir <- paste(path, taxa, sep = "")
  existingDirCheck(path = dir)
  
  if(normalizeCounts){
    #Get normalized counts from DESeq2 object
    normalized_counts <- counts(deseq, normalized = TRUE)
    print(normalized_counts)
  }else{
    # ps <- transformCounts(ps, transformation = "rel_ab")
    normalized_counts <- t(otu_table(ps))
    # # Ensure consistent structure: convert to a matrix
    # normalized_counts <- as.data.frame(normalized_counts)
  }
  
  #Define empty list that will contain pairs comparaisons names
  vs <- c()
  
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
    # res_subset1 <- results(deseq, contrast = list(c(resultsNames(deseq)[2]))) #wt putrescine vs vehicle / # 50 vs 500 ctrl
    # sigtab_1 <- cbind(as(res_subset1, "data.frame"), as(tax_table(ps)[rownames(res_subset1), ], "matrix"))
    # sigtab_1$comparaison <- 1
    # sigtab_1$vs <- vs[1]
    # 
    # res_subset2 <- results(deseq, contrast = list(c(resultsNames(deseq)[3], resultsNames(deseq)[4]))) #il22 ko putrescine vs vehicle / 50 dss vs 500 dss
    # sigtab_2 <- cbind(as(res_subset2, "data.frame"), as(tax_table(ps)[rownames(res_subset2), ], "matrix"))
    # sigtab_2$comparaison <- 2
    # sigtab_2$vs <- vs[2]
    # 
    # res_subset3 <- results(deseq, contrast = list(c(resultsNames(deseq)[3]))) #vehicle wt vs il22 ko / 50 water vs 50 dss
    # sigtab_3 <- cbind(as(res_subset3, "data.frame"), as(tax_table(ps)[rownames(res_subset3), ], "matrix"))
    # sigtab_3$comparaison <- 3
    # sigtab_3$vs <- vs[3]
    # 
    # res_subset4 <- results(deseq, contrast = list(c(resultsNames(deseq)[2], resultsNames(deseq)[4]))) #putrescine wt vs il22 ko / 500 water vs 500 dss
    # sigtab_4 <- cbind(as(res_subset4, "data.frame"), as(tax_table(ps)[rownames(res_subset4), ], "matrix"))
    # sigtab_4$comparaison <- 4
    # sigtab_4$vs <- vs[4]
    # 
    # interaction <- results(deseq, contrast= list(c(resultsNames(deseq)[1]))) # Do genotypes respond differently to treatment, comparisons of comparisons
    # sigtab_interaction <- cbind(as(interaction, "data.frame"), as(tax_table(ps)[rownames(interaction), ], "matrix"))
    # sigtab_interaction$comparaison <- 5
    # sigtab_interaction$vs <- "interaction" 
    
    res_subset1 <- results(deseq, contrast = list(c(resultsNames(deseq)[2]))) #wt putrescine vs vehicle / # 50 vs 500 ctrl
    sigtab_1 <- cbind(as(res_subset1, "data.frame"), as(tax_table(ps)[rownames(res_subset1), ], "matrix"))
    sigtab_1$comparaison <- 1
    sigtab_1$vs <- vs[1]
    
    res_subset2 <- results(deseq, contrast = list(c(resultsNames(deseq)[2], resultsNames(deseq)[4]))) #il22 ko putrescine vs vehicle / 50 dss vs 500 dss
    sigtab_2 <- cbind(as(res_subset2, "data.frame"), as(tax_table(ps)[rownames(res_subset2), ], "matrix"))
    sigtab_2$comparaison <- 2
    sigtab_2$vs <- vs[2]
    
    res_subset3 <- results(deseq, contrast = list(c(resultsNames(deseq)[3]))) #vehicle wt vs il22 ko / 50 water vs 50 dss
    sigtab_3 <- cbind(as(res_subset3, "data.frame"), as(tax_table(ps)[rownames(res_subset3), ], "matrix"))
    sigtab_3$comparaison <- 3
    sigtab_3$vs <- vs[3]
    
    res_subset4 <- results(deseq, contrast = list(c(resultsNames(deseq)[3], resultsNames(deseq)[4]))) #putrescine wt vs il22 ko / 500 water vs 500 dss
    sigtab_4 <- cbind(as(res_subset4, "data.frame"), as(tax_table(ps)[rownames(res_subset4), ], "matrix"))
    sigtab_4$comparaison <- 4
    sigtab_4$vs <- vs[4]
    
    # interaction <- results(deseq, contrast= list(c(resultsNames(deseq)[1]))) # Do genotypes respond differently to treatment, comparisons of comparisons
    # sigtab_interaction <- cbind(as(interaction, "data.frame"), as(tax_table(ps)[rownames(interaction), ], "matrix"))
    # sigtab_interaction$comparaison <- 5
    # sigtab_interaction$vs <- "interaction" 
    
    
  }else{ # If you performed deseq with 2 variables, each with 2 groups
    
    #Partition results for specific pairwise comparaisons
    res_subset1 <- results(deseq, contrast = list(resultsNames(deseq)[3])) #wt putrescine vs vehicle // 50vs500 ctrl
    sigtab_1 <- cbind(as(res_subset1, "data.frame"), as(tax_table(ps)[rownames(res_subset1), ], "matrix"))
    sigtab_1$comparaison <- 1
    sigtab_1$vs <- vs[1]

    res_subset2 <- results(deseq, contrast = list(c(resultsNames(deseq)[3], resultsNames(deseq)[4]))) #il22 ko putrescine vs vehicle // 50vs500 dss
    sigtab_2 <- cbind(as(res_subset2, "data.frame"), as(tax_table(ps)[rownames(res_subset2), ], "matrix"))
    sigtab_2$comparaison <- 2
    sigtab_2$vs <- vs[2]

    res_subset3 <- results(deseq, contrast = list(resultsNames(deseq)[2])) #vehicle wt vs il22 ko // 50 ctrl vs 50 dss
    sigtab_3 <- cbind(as(res_subset3, "data.frame"), as(tax_table(ps)[rownames(res_subset3), ], "matrix"))
    sigtab_3$comparaison <- 3
    sigtab_3$vs <- vs[3]

    res_subset4 <- results(deseq, contrast = list(c(resultsNames(deseq)[2], resultsNames(deseq)[4]))) #putrescine wt vs il22 ko // 500 ctrl vs 500 dss
    sigtab_4 <- cbind(as(res_subset4, "data.frame"), as(tax_table(ps)[rownames(res_subset4), ], "matrix"))
    sigtab_4$comparaison <- 4
    sigtab_4$vs <- vs[4]
    # 
    # interaction <- results(deseq, contrast= list(c(resultsNames(deseq)[4]))) # Do genotypes respond differently to treatment, comparisons of comparisons
    # sigtab_interaction <- cbind(as(interaction, "data.frame"), as(tax_table(ps)[rownames(interaction), ], "matrix"))
    # sigtab_interaction$comparaison <- 5
    # sigtab_interaction$vs <- "interaction" 
    
    # [1] "Intercept"              "treatment_dss_vs_water" "diet_500_vs_50"        
    # [4] "cage"                   "treatmentdss.diet500" 
    
    # res_subset1 <- results(deseq, contrast = list(resultsNames(deseq)[3])) # 50 vs 500 ctrl
    # sigtab_1 <- cbind(as(res_subset1, "data.frame"), as(tax_table(ps)[rownames(res_subset1), ], "matrix"))
    # sigtab_1$comparaison <- 1
    # sigtab_1$vs <- vs[1]
    # 
    # res_subset2 <- results(deseq, contrast = list(c(resultsNames(deseq)[3], resultsNames(deseq)[5]))) # 50 dss vs 500 dss
    # sigtab_2 <- cbind(as(res_subset2, "data.frame"), as(tax_table(ps)[rownames(res_subset2), ], "matrix"))
    # sigtab_2$comparaison <- 2
    # sigtab_2$vs <- vs[2]
    # 
    # res_subset3 <- results(deseq, contrast = list(resultsNames(deseq)[2])) # 50 water vs 50 dss
    # sigtab_3 <- cbind(as(res_subset3, "data.frame"), as(tax_table(ps)[rownames(res_subset3), ], "matrix"))
    # sigtab_3$comparaison <- 3
    # sigtab_3$vs <- vs[3]
    # 
    # res_subset4 <- results(deseq, contrast = list(c(resultsNames(deseq)[2], resultsNames(deseq)[5]))) # 500 water vs 500 dss
    # sigtab_4 <- cbind(as(res_subset4, "data.frame"), as(tax_table(ps)[rownames(res_subset4), ], "matrix"))
    # sigtab_4$comparaison <- 4
    # sigtab_4$vs <- vs[4]

    # interaction <- results(deseq, contrast= list(c(resultsNames(deseq)[4]))) # cage?
    # sigtab_interaction <- cbind(as(interaction, "data.frame"), as(tax_table(ps)[rownames(interaction), ], "matrix"))
    # sigtab_interaction$comparaison <- 5
    # sigtab_interaction$vs <- "interaction"
    
  }
    
  #Append the sigtabs together
  sigtab <- bind_rows(sigtab_1, sigtab_2, sigtab_3, sigtab_4) #sigtab_interaction
  
  #Add ASV variable col to sigtab (enables to store asv names, not only as rownames, because they will be changed when using rowbind)
  sigtab["asv"] <- gsub("\\..*", "", rownames(sigtab))
  
  #Keeping only ASVs for which they were taxa found at the taxonomical level of interest
  if(taxa == "Species" & includeUnknownSpecies){
    sigtab[[taxa]][is.na(sigtab[[taxa]])] <- "Unknown" # If includeUnknownSpecies, keep ASVs for which there was Genus tax annotation even if no species annotation
    sigtab <- sigtab[!is.na(sigtab[["Genus"]]),]
  }else{
    sigtab <- sigtab[!is.na(sigtab[[taxa]]),]
  }
  
  
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
  
  
  #Find asvs that have at least one significant value at some timepoint
  asvList <- unique(sigtab[(sigtab[pvalue])<threshold,"asv"])

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
      sigtab_taxon <- sigtab[sigtab$asv == asv,]
      
      #Add the 16s sequence to the sigtab specific to ASV of interest (only for species, higher taxonomical levels have multiple ASVs associated with them)
      if(taxa == "Species"){
        sigtab_taxon$dna_sequence <- as.character(refseq(ps)[asv])
      }
      
      #Creates directory with taxon name
      dir_taxon <- paste(dir, "/", i, "-",taxonName, sep = "")
      existingDirCheck(path = dir_taxon)
      
      ##Calculate relative abundance as a percentage for each sample
      relative_abundance <- apply(normalized_counts, 2, prop.table) * 100
        #normalized_counts[asv, ]
      
      #Keep for taxa of interest
      relative_abundance <- relative_abundance[asv, ]
      # relative_abundance <- asv_counts * 100 / colSums(normalized_counts)
      
      #Convert relative abundance to a data frame
      relative_abundance <- data.frame(
        sample_id = names(relative_abundance),
        rel_ab = as.numeric(relative_abundance)
      )
      
      #Merge relative abundance with sample metadata
      relative_abundance <- merge(relative_abundance, as(sample_data(ps), "data.frame"), by = "sample_id")
      
      #Save group names for using them in the ggsignif
      groups <- levels(sample_data(ps)[[gg_group]])
      
      p <- ggplot(data = relative_abundance, aes(x = .data[[gg_group]], y = rel_ab, color = .data[[gg_group]])) +
        geom_point(size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) + 
        
        #Error bars
        stat_summary(fun.data = "mean_se", geom = "errorbar",
                     aes(color = .data[[gg_group]]),
                     width = 0.2, size = 0.7,
                     position = position_dodge(-0.75)) +
        
        #Mean lines
        stat_summary(fun.data = "mean_se", geom = "errorbar",
                     aes(ymin = ..y.., ymax = ..y.., group = .data[[gg_group]]),
                     color = "black", linewidth = 0.5, width = 0.5,
                     position = position_dodge(-0.75))+
        
        
        labs(title = taxonName,
             y = "Relative abundance (%)", color = "Groups", x = "Groups") +
        scale_color_manual(values = customColors)+
        theme_minimal()+
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
          axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bar
        
        if(displaySignificance){
          #Add significance bars
          p <- p + geom_signif(comparisons = list(c(groups[1],groups[2])),
                      annotations = ifelse(displayPvalue, paste("p = ", 
                                                                format(sigtab_taxon[sigtab_taxon$comparaison == 1, pvalue], digits = 2, scientific = TRUE)),
                                           sigtab_taxon[sigtab_taxon$comparaison == 1, "significance"]),
                      tip_length = 0.02,
                      y_position =  max(relative_abundance$rel_ab)+1/12*max(relative_abundance$rel_ab),
                      size = 1.2,  # Make the bar wider
                      color = "black") +
            
            geom_signif(comparisons = list(c(groups[3],groups[4])),
                        annotations = ifelse(displayPvalue, paste("p = ", 
                                                                  format(sigtab_taxon[sigtab_taxon$comparaison == 2, pvalue], digits = 2, scientific = TRUE)),
                                             sigtab_taxon[sigtab_taxon$comparaison == 2, "significance"]),
                        tip_length = 0.02,
                        y_position =  max(relative_abundance$rel_ab)+1/12*max(relative_abundance$rel_ab),
                        size = 1.2,  # Make the bar wider
                        color = "black") +
            
            geom_signif(comparisons = list(c(groups[1],groups[3])),
                        annotations = ifelse(displayPvalue, paste("p = ", 
                                                                  format(sigtab_taxon[sigtab_taxon$comparaison == 3, pvalue], digits = 2, scientific = TRUE)),
                                             sigtab_taxon[sigtab_taxon$comparaison == 3, "significance"]),
                        tip_length = 0.02,
                        y_position =  max(relative_abundance$rel_ab)+2/12*max(relative_abundance$rel_ab),
                        size = 1.2,  # Make the bar wider
                        color = "black") +
            
            geom_signif(comparisons = list(c(groups[2],groups[4])),
                        annotations = ifelse(displayPvalue, paste("p = ", 
                                                                  format(sigtab_taxon[sigtab_taxon$comparaison == 4, pvalue], digits = 2, scientific = TRUE)),
                                             sigtab_taxon[sigtab_taxon$comparaison == 4, "significance"]),
                        tip_length = 0.02,
                        y_position =  max(relative_abundance$rel_ab)+3/12*max(relative_abundance$rel_ab),
                        size = 1.2,  # Make the bar wider
                        color = "black")
        }
        
      if(isFALSE(is.null(additionnalAes))){
        # p <- do.call("+", c(list(p), additionnalAes))
        
        p <- Reduce("+", c(list(p), additionnalAes))
      }
      
      ggsave(plot = p, filename = paste(dir_taxon,"/",taxonName,"_relab.png", sep = ""), dpi = 300, height = dim[1], width = dim[2], bg = 'white')
      
      #Write as excel file the significance table specific to an ASV
      write.xlsx(sigtab_taxon, paste(dir_taxon,"/",gsub(" ", "_", taxonName),"_stats.xlsx", sep = ""))

      #Write as excel file the relative abundance data specific  to an ASV
      write.xlsx(relative_abundance, paste(dir_taxon,"/",gsub(" ", "_", taxonName),"_relab.xlsx", sep = ""))
      
    }

  #Returns list of significant asvs
  if(returnSigAsvs){
    return(asvList)
  }
    
}

#For design with 4 groups based on 2 conditions - this the latest version used
#gg_group must be order with correct order prior to that (as a factor)
relabSingleGroup <- function(ps, deseq, measure = "log2fold", gg_group, taxa = "Species", threshold = 0.01, displayPvalue = FALSE, returnSigAsvs = FALSE, normalizeCounts = FALSE, customColors, pairs, path){
  
  #Creates directory for taxonomic level
  dir <- paste(path, taxa, sep = "")
  existingDirCheck(path = dir)
  
  if(normalizeCounts){
    #Get normalized counts from DESeq2 object
    normalized_counts <- counts(deseq, normalized = TRUE)
    print(normalized_counts)
  }else{
    # ps <- transformCounts(ps, transformation = "rel_ab")
    normalized_counts <- t(otu_table(ps))
    # # Ensure consistent structure: convert to a matrix
    # normalized_counts <- as.data.frame(normalized_counts)
    print(normalized_counts)
  }
  
  
  
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
  res_subset1 <- results(deseq, contrast = list(resultsNames(deseq)[3])) #wt putrescine vs vehicle
  sigtab_1 <- cbind(as(res_subset1, "data.frame"), as(tax_table(ps)[rownames(res_subset1), ], "matrix"))
  sigtab_1$comparaison <- 1
  sigtab_1$vs <- vs[1]
  
  res_subset2 <- results(deseq, contrast=list(c(resultsNames(deseq)[3], resultsNames(deseq)[4]))) #il22 ko putrescine vs vehicle
  sigtab_2 <- cbind(as(res_subset2, "data.frame"), as(tax_table(ps)[rownames(res_subset2), ], "matrix"))
  sigtab_2$comparaison <- 2
  sigtab_2$vs <- vs[2]
  
  res_subset3 <- results(deseq, contrast=list(resultsNames(deseq)[2])) #vehicle wt vs il22 ko
  sigtab_3 <- cbind(as(res_subset3, "data.frame"), as(tax_table(ps)[rownames(res_subset3), ], "matrix"))
  sigtab_3$comparaison <- 3
  sigtab_3$vs <- vs[3]
  
  res_subset4 <- results(deseq, contrast=list(c(resultsNames(deseq)[2], resultsNames(deseq)[4]))) #putrescine wt vs il22 ko
  sigtab_4 <- cbind(as(res_subset4, "data.frame"), as(tax_table(ps)[rownames(res_subset4), ], "matrix"))
  sigtab_4$comparaison <- 4
  sigtab_4$vs <- vs[4]
  
  interaction <- results(deseq, contrast=list(resultsNames(deseq)[4])) # Do genotypes respond differently to treatment, comparisons of comparisons
  sigtab_interaction <- cbind(as(interaction, "data.frame"), as(tax_table(ps)[rownames(interaction), ], "matrix"))
  sigtab_interaction$comparaison <- 5
  sigtab_interaction$vs <- "interaction" 
  
  #Append the sigtabs together
  sigtab <- bind_rows(sigtab_1, sigtab_2, sigtab_3, sigtab_4, sigtab_interaction)
  
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
    sigtab_taxon <- sigtab[sigtab$asv == asv,]
    
    #Add the 16s sequence to the sigtab specific to ASV of interest (only for species, higher taxonomical levels have multiple ASVs associated with them)
    if(taxa == "Species"){
      sigtab_taxon$dna_sequence <- as.character(refseq(ps)[asv])
    }
    
    #Creates directory with taxon name
    dir_taxon <- paste(dir, "/", i, "-",taxonName, sep = "")
    existingDirCheck(path = dir_taxon)
    
    ##Calculate relative abundance as a percentage for each sample
    relative_abundance <- apply(normalized_counts, 2, prop.table) * 100
    #normalized_counts[asv, ]
    
    #Keep for taxa of interest
    relative_abundance <- relative_abundance[asv, ]
    # relative_abundance <- asv_counts * 100 / colSums(normalized_counts)
    
    #Convert relative abundance to a data frame
    relative_abundance <- data.frame(
      sample_id = names(relative_abundance),
      rel_ab = as.numeric(relative_abundance)
    )
    
    #Merge relative abundance with sample metadata
    relative_abundance <- merge(relative_abundance, as(sample_data(ps), "data.frame"), by = "sample_id")
    
    #Save group names for using them in the ggsignif
    groups <- levels(sample_data(ps)[[gg_group]])
    
    p <- ggplot(data = relative_abundance, aes(x = gg_group, y = rel_ab, color = gg_group)) +
      geom_point(size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) + 
      
      #Error bars
      stat_summary(fun.data = "mean_se", geom = "errorbar",
                   aes(color = gg_group),
                   width = 0.2, size = 0.7,
                   position = position_dodge(-0.75)) +
      
      #Mean lines
      stat_summary(fun.data = "mean_se", geom = "errorbar",
                   aes(ymin = ..y.., ymax = ..y.., group = gg_group),
                   color = "black", linewidth = 0.5, width = 0.5,
                   position = position_dodge(-0.75))+
      
      
      labs(title = taxonName,
           y = "Relative abundance (%)", color = "Groups", x = "Groups") +
      scale_color_manual(values = customColors)+
      
      #Add significance bars
      geom_signif(comparisons = list(c(groups[1],groups[2])),
                  annotations = ifelse(displayPvalue, paste("p = ", 
                                                            format(sigtab_taxon[sigtab_taxon$comparaison == 1, "padj"], digits = 2, scientific = TRUE)),
                                       sigtab_taxon[sigtab_taxon$comparaison == 1, "significance"]),
                  tip_length = 0.02,
                  y_position =  max(relative_abundance$rel_ab)+1/12*max(relative_abundance$rel_ab),
                  size = 1.2,  # Make the bar wider
                  color = "black") +
      
      geom_signif(comparisons = list(c(groups[3],groups[4])),
                  annotations = ifelse(displayPvalue, paste("p = ", 
                                                            format(sigtab_taxon[sigtab_taxon$comparaison == 2, "padj"], digits = 2, scientific = TRUE)),
                                       sigtab_taxon[sigtab_taxon$comparaison == 2, "significance"]),
                  tip_length = 0.02,
                  y_position =  max(relative_abundance$rel_ab)+1/12*max(relative_abundance$rel_ab),
                  size = 1.2,  # Make the bar wider
                  color = "black") +
      
      geom_signif(comparisons = list(c(groups[1],groups[3])),
                  annotations = ifelse(displayPvalue, paste("p = ", 
                                                            format(sigtab_taxon[sigtab_taxon$comparaison == 3, "padj"], digits = 2, scientific = TRUE)),
                                       sigtab_taxon[sigtab_taxon$comparaison == 3, "significance"]),
                  tip_length = 0.02,
                  y_position =  max(relative_abundance$rel_ab)+2/12*max(relative_abundance$rel_ab),
                  size = 1.2,  # Make the bar wider
                  color = "black") +
      
      geom_signif(comparisons = list(c(groups[2],groups[4])),
                  annotations = ifelse(displayPvalue, paste("p = ", 
                                                            format(sigtab_taxon[sigtab_taxon$comparaison == 4, "padj"], digits = 2, scientific = TRUE)),
                                       sigtab_taxon[sigtab_taxon$comparaison == 4, "significance"]),
                  tip_length = 0.02,
                  y_position =  max(relative_abundance$rel_ab)+3/12*max(relative_abundance$rel_ab),
                  size = 1.2,  # Make the bar wider
                  color = "black") +
      
      theme_minimal()+
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
        axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bar
    ggsave(plot = p, filename = paste(dir_taxon,"/",taxonName,"_relab.png", sep = ""), dpi = 300, height = 6, width = 6, bg = 'white')
    
    #Write as excel file the significance table specific to an ASV
    write.xlsx(sigtab_taxon, paste(dir_taxon,"/",gsub(" ", "_", taxonName),"_stats.xlsx", sep = ""))
    
    #Write as excel file the relative abundance data specific  to an ASV
    write.xlsx(relative_abundance, paste(dir_taxon,"/",gsub(" ", "_", taxonName),"_relab.xlsx", sep = ""))
    
  }
  
  #Returns list of significant asvs
  if(returnSigAsvs){
    return(asvList)
  }
  
}

#Revised function that does deseq analysis but only for one factor with two groups (ex: diet 50 vs 500)
relabSingleTimepoint <- function(ps, deseq, measure = "log2fold", varToCompare, timePoint, taxa = "Species", threshold = 0.01, FDR = TRUE,
                                 LDA = FALSE, customColors, path, additionnalAes = NULL, dim = c(6,6), displayPvalue = TRUE, blockFactor = FALSE, displaySampleID = FALSE){
  
  #Creates directory for taxonomic level
  dir <- paste(path, taxa, sep = "")
  existingDirCheck(path = dir)
  
  # Get counts from DESeq2 object
  normalized_counts <- counts(deseq, normalized = FALSE)
  
  #Save results at single timepoint
  if(blockFactor){ # Block factor is used to explain some variability, but overall results are used
    res <- results(deseq, name = resultsNames(deseq)[3])
  }else{
    res <- results(deseq, name = resultsNames(deseq)[2])
  }
  
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
  sigtab$significance <- cut(sigtab[,pvalue],
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", "NS"))
  
  #Find asvs that have at least one significant value at some timepoint
  asvList <- unique(sigtab[(sigtab[pvalue])<threshold,"asv"])

  if(length(asvList)==0){
    return(NULL)
  }
  
  #Loop along ASVs (single taxons analyzed)
  for(i in seq_along(asvList)){
    
    # Save asv value
    asv = asvList[i]

    # Save speciesName for dir creation and for graph title
    if(taxa == "Species"){
      taxonName <- paste(unique(sigtab[sigtab$asv == asv, "Genus"]),unique(sigtab[sigtab$asv == asv, "Species"]), paste("(", asv, ")", sep =""), sep = " ")
    }else{
      taxonName <- unique(sigtab[sigtab$asv == asv, taxa])
    }

    # Sigtab specific to the taxon analyzed
    sigtab_taxon <- sigtab[sigtab$asv == asv,]

    # Add the 16s sequence to the sigtab specific to ASV of interest (only for species, higher taxonomical levels have multiple ASVs associated with them)
    if(taxa == "Species"){
      sigtab_taxon$dna_sequence <- as.character(refseq(ps)[asv])
    }

    #Save p-value
    p_value <- sigtab_taxon[,pvalue]
    significance <- sigtab_taxon["significance"]

    #Creates directory with taxon name
    dir_taxon <- paste(dir, "/", i, "-",taxonName, sep = "")
    existingDirCheck(path = dir_taxon)

    #Select counts for ASV of interest
    asv_counts <- normalized_counts[asv, ]
    
    #Calculate relative abundance as a percentage for each sample
    relative_abundance <- asv_counts * 100 / colSums(normalized_counts)

    #Convert relative abundance to a data frame
    relative_abundance <- data.frame(
      sample_id = names(relative_abundance),
      rel_ab = as.numeric(relative_abundance)
    )

    #Merge relative abundance with sample metadata
    relative_abundance <- merge(relative_abundance, as(sample_data(ps), "data.frame"), by = "sample_id")

    # relative_abundance[[timeVariable]] <- as.numeric(as.character(relative_abundance[[timeVariable]])) * 7
    relative_abundance[[varToCompare]] <- as.character(relative_abundance[[varToCompare]])
    relative_abundance$week <- timePoint

    #Saving groups variables as list of levels for varToCompare
    groups <- levels(sample_data(ps)[[varToCompare]])
    


    p <- ggplot(data = relative_abundance, aes(x = !!sym(varToCompare), y = rel_ab, color = !!sym(varToCompare)), group = !!sym(varToCompare)) +

      geom_point(size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) +

      # Error bars
      stat_summary(fun.data = "mean_se", geom = "errorbar",
                   aes(color = !!sym(varToCompare)),
                   width = 0.5, size = 1,
                   alpha = 0.5,
                   position = "identity")+

      #Mean lines
      stat_summary(fun.data = "mean_se", geom = "errorbar",
                   aes(ymin = ..y.., ymax = ..y.., group = !!sym(varToCompare)),
                   color = "black", linewidth = 0.5, width = 1,
                   position = "identity")+


      labs(title = paste(taxa ,taxonName, "\nat last timepoint", sep = " "), y = "Relative abundance (%)", color = "Diet", x = "") +
      scale_color_manual(values = customColors)+
      # scale_x_discrete(labels = c("50 ppm\nDSS","500 ppm\nDSS"))+
      scale_x_discrete(labels = c(levels(sample_data(ps)[[varToCompare]])))+

      #Add stat bar
      geom_signif(comparisons = list(c(groups[1],groups[2])),
                  annotations = ifelse(displayPvalue, paste("p = ", format(p_value, digits = 2, scientific = TRUE)), significance),
                  tip_length = 0.02,
                  y_position =  max(relative_abundance$rel_ab),
                  size = 1.2,  # Make the bar wider
                  color = "black") +
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
    
    # Add sample IDs as text labels if displaySample is TRUE
    if (displaySampleID) {
      p <- p+geom_text(aes(label = sample_id), position = position_jitter(width = 0.1, height = 0), size = 3, vjust = -0.5)
    }
    
    if(isFALSE(is.null(additionnalAes))){
      p <- p + additionnalAes
    }
    
    ggsave(plot = p, filename = paste(dir_taxon,"/",gsub(" ", "_", taxonName),"_relab.png", sep = ""), dpi = 300, height = dim[1], width = dim[2], bg = 'white')

    #Write as excel file the significance table specific to an ASV
    write.xlsx(sigtab_taxon, paste(dir_taxon,"/",gsub(" ", "_", taxonName),"_stats.xlsx", sep = ""))

    #Write as excel file the relative abundance data specific  to an ASV
    write.xlsx(relative_abundance, paste(dir_taxon,"/",gsub(" ", "_", taxonName),"_relab.xlsx", sep = ""))

  }
  
  if(LDA){
    
    # Apply rlog transformation to the count data
    rlog_data <- rlog(deseq, blind = TRUE)
    
    # Extract the log-transformed counts for the significant features
    lda_data <- assay(rlog_data)[asvList, ]
    
    # Combine the data with the group labels (e.g., condition)
    lda_input <- merge(t(lda_data), sample_data(ps)[,varToCompare], by = "row.names")
    row.names(lda_input) <- lda_input$Row.names
    lda_input <- lda_input[-1]
  
    
    train_control <- trainControl(method = "cv", number = 10)
    # Train the LDA model using caret (this returns a train object)
    lda_model <- train(as.formula(paste(varToCompare, "~ .")),
                       data = lda_input,
                       method = "lda",
                       trControl = train_control)  # ensure you have defined train_control
    
    # Extract ASV-level LDA coefficients (effect sizes)
    lda_effects <- data.frame(
      ASV = rownames(lda_model$finalModel$scaling),
      LDA_Score = lda_model$finalModel$scaling[,1],  # Use first discriminant
      taxa = paste0(substring(tax_table(ps)[rownames(lda_model$finalModel$scaling),"Genus"], 1, 1), ". ",
                                    tax_table(ps)[rownames(lda_model$finalModel$scaling),"Species"])
    )
    
    lda_summary <- lda_effects %>%
      group_by(taxa) %>%
      summarise(mean_LDA = mean(LDA_Score)) %>%
      mutate(Group = ifelse(mean_LDA > 0, levels(sample_data(ps)[[varToCompare]])[2], levels(sample_data(ps)[[varToCompare]])[1]))  # Separate step for clarity
    
    print(lda_summary)
    
    # Plot species-level mean LDA scores
    p <- ggplot(lda_summary, aes(x = reorder(taxa, mean_LDA), y = mean_LDA, fill = Group)) +
      geom_col() +
      coord_flip() +
      geom_hline(yintercept = 0, color = "black") +
      scale_fill_manual(values = customColors) +
      labs(x = "", y = "LDA Score") +
      theme_minimal() +
      theme(legend.position = "top",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.x = element_text(face = "bold"),
            axis.text.x = element_text(face = "bold"),
            axis.text.y = element_text(face = "bold.italic")
            )
    
    ggsave(plot = p, filename = paste0(dir, "/lda.png"), dpi = 300, bg = "white", width = 4, height = 6)
  }
}

#Produces the timeline with stats calculated with same fashion as for relabSingleTimepoint, works with design "~ factor", and is done separately for each week
relabTimelineRevised <- function(ps, measure = "log2fold", timeVariable, varToCompare, taxa = "Species", threshold = 0.01, returnSigAsvs = FALSE, customColors, path, displayIndivValues = FALSE, dim = c(5,5)){

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
    normalized_counts <- counts(deseq_subset, normalized = FALSE)
    
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
    
    #Before merging sigtab, ensure timeVariable is numeric
    sigtab_taxon[[timeVariable]] <- as.numeric(as.character(sigtab_taxon[[timeVariable]])) *7
    
    #Merge sigtab_taxon with means_df
    means_df <- merge(sigtab_taxon, means_df, by = timeVariable)
    
    p <- ggplot(data = relative_abundance, aes(x = !!sym(timeVariable), y = rel_ab, color = !!sym(varToCompare)), group = !!sym(varToCompare)) +
      scale_x_continuous(n.breaks = 8)+
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
      # theme(
      #   plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
      #   axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
      #   axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
      #   axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust x-axis tick label font size
      #   axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
      #   legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
      #   legend.text = element_text(size = 12),  # Adjust legend font size
      #   panel.grid.major = element_blank(),  # Add major grid lines
      #   panel.grid.minor = element_blank(),  # Remove minor grid lines
      #   axis.line = element_line(color = "black", size = 1),
      #   panel.background = element_blank()) # Include axis lines  # Include axis bar
      my_theme()
    
    # Show individual values
    if(displayIndivValues){
      p <- p + geom_point(size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) 
    }
    
    ggsave(plot = p, filename = paste(dir_taxon,"/",gsub(" ", "_", taxonName),"_relab.png", sep = ""), dpi = 300, height = dim[1], width = dim[2], bg = 'white')
    
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

volcanoPlot2Groups <- function(ps, deseq, varToCompare, taxa = "Species", threshold = 0.05, customColors){
  
  #Save results at single timepoint
  res <- results(deseq, name = resultsNames(deseq)[2])

  #Save significance table
  sigtab <- cbind(as(res, "data.frame"), as(tax_table(ps)[rownames(res),], "matrix"))
  
  #Add ASV variable col to sigtab (enables to store asv names, not only as rownames, because they will be changed when using rowbind)
  sigtab["asv"] <- rownames(sigtab)
  
  #Keeping only ASVs for which they were taxa found at the taxonomical level of interest
  sigtab <- sigtab[!is.na(sigtab[[taxa]]),]
  
  #Replacing NA padj by 1 (they correspond to this anyways)
  sigtab$padj[is.na(sigtab$padj)] <- 1
  
  res_df <- sigtab %>%
    mutate(significance = case_when(
      padj < 0.05 & log2FoldChange >= 1  ~ "Up",
      padj < 0.05 & log2FoldChange <= -1 ~ "Down",
      TRUE                               ~ "NotSig"
    ))
  
  EnhancedVolcano(
    res_df,
    lab = res_df[[taxa]],
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    FCcutoff = 1,
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~ 'adjusted p'),
    legendPosition = 'right',
    title = "50 ppm vs 500 ppm at T35",
    subtitle = taxa
  )
  
}

volcanoPlot2GroupsMultifactorDesign <- function(ps, deseq, varToCompare, taxa = "Species", threshold = 0.05, FCcutoff = 1,
                                                customColors, FDR = TRUE, includeUnknownSpecies = TRUE, selectedComparison = 1,
                                                title = ""){
  
  #Partition results for specific pairwise comparaisons
  res_subset1 <- results(deseq, contrast = list(resultsNames(deseq)[3])) #wt putrescine vs vehicle // 50vs500 ctrl
  sigtab_1 <- cbind(as(res_subset1, "data.frame"), as(tax_table(ps)[rownames(res_subset1), ], "matrix"))
  sigtab_1$comparaison <- 1
  
  res_subset2 <- results(deseq, contrast = list(c(resultsNames(deseq)[3], resultsNames(deseq)[4]))) #il22 ko putrescine vs vehicle // 50vs500 dss
  sigtab_2 <- cbind(as(res_subset2, "data.frame"), as(tax_table(ps)[rownames(res_subset2), ], "matrix"))
  sigtab_2$comparaison <- 2
  
  res_subset3 <- results(deseq, contrast = list(resultsNames(deseq)[2])) #vehicle wt vs il22 ko // 50 ctrl vs 50 dss
  sigtab_3 <- cbind(as(res_subset3, "data.frame"), as(tax_table(ps)[rownames(res_subset3), ], "matrix"))
  sigtab_3$comparaison <- 3
  
  res_subset4 <- results(deseq, contrast = list(c(resultsNames(deseq)[2], resultsNames(deseq)[4]))) #putrescine wt vs il22 ko // 500 ctrl vs 500 dss
  sigtab_4 <- cbind(as(res_subset4, "data.frame"), as(tax_table(ps)[rownames(res_subset4), ], "matrix"))
  sigtab_4$comparaison <- 4
  
  # Append the sigtabs together
  sigtab <- bind_rows(sigtab_1, sigtab_2, sigtab_3, sigtab_4) 
  
  # Add ASV variable col to sigtab (enables to store asv names, not only as rownames, because they will be changed when using rowbind)
  sigtab["asv"] <- gsub("\\..*", "", rownames(sigtab))
  
  # Keeping only ASVs for which they were taxa found at the taxonomical level of interest
  if(taxa == "Species" & includeUnknownSpecies){
    sigtab[[taxa]][is.na(sigtab[[taxa]])] <- "unknown" # If includeUnknownSpecies, keep ASVs for which there was Genus tax annotation even if no species annotation
    sigtab <- sigtab[!is.na(sigtab[["Genus"]]),]
  }else{
    sigtab <- sigtab[!is.na(sigtab[[taxa]]),]
  }
  
  # If taxa == species prepare names that are displayed 
  if(taxa == "Species"){
    sigtab[[taxa]] <- paste(sigtab[["Genus"]],sigtab[[taxa]], sep = "\n")
  }
  
  # FDR corrected pvalues or uncorrected pvalues
  if(FDR){
    pvalue <- as.character("padj")
  }else{
    pvalue <- as.character("pvalue")
  }
  
  # Keeping only ASVs for which they were taxa found at the taxonomical level of interest
  sigtab <- sigtab[!is.na(sigtab[[taxa]]),]
  
  # Replacing NA padj by 1 (they correspond to this anyways)
  sigtab$padj[is.na(sigtab$padj)] <- 1
  
  # Select values of interest
  sigtab <- sigtab[sigtab$comparaison == selectedComparison,]
  
  res_df <- sigtab %>%
    mutate(significance = case_when(
      padj < threshold & log2FoldChange >= FCcutoff  ~ "Up",
      padj < threshold & log2FoldChange <= -FCcutoff ~ "Down",
      TRUE                               ~ "NotSig"
    ))
  
  EnhancedVolcano(
    res_df,
    lab = res_df[[taxa]],
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = threshold,
    FCcutoff = FCcutoff,
    xlab = bquote(bold(~Log[2]~ 'fold change')),
    ylab = bquote(bold(~-Log[10]~ 'adjusted p')),
    legendPosition = 'bottom',
    title = title,
    subtitle = NULL,
    ylim = c(NA,max(-log10(res_df$padj))),
    xlim = c(-max(abs(res_df$log2FoldChange)),max(abs(res_df$log2FoldChange))),
    caption = NULL,
    labSize = 2,
    legendLabels = c("n.s.",bquote(bold(Log[2]~ 'FC')),"C",bquote(bold('p-value and'~Log[2]~ 'FC'))),
    legendIconSize = 2

    
    
  )
  
}