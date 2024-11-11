#Revised version, not doing separate distance calculations and adding path option

#Beta diversity analysis for different timepoints. You must provide a filtered ps object, the timeVariable and the varToCompare (present in sample_data)
betaDiversityTimepoint <- function(ps, timeVariable, varToCompare, distMethod, customColors, path){
  
  #calculating distance matrix
  dist <- phyloseq::distance(ps, method = distMethod)
  
  for(timepoint in levels(sample_data(ps)[[timeVariable]])){
    
    #calculating distance matrix
    dist_subset <- as.dist(as.matrix(dist)[sample_data(ps)$sample_id[sample_data(ps)[[timeVariable]] %in% timepoint],sample_data(ps)$sample_id[sample_data(ps)[[timeVariable]] %in% timepoint]])
    
    #Creates subset for each different timepoint
    ps_subset <- prune_samples(sample_data(ps)[[timeVariable]] == timepoint, ps)
    
    #Define dir path where graphs and stats are gonna be saved
    dir <- paste(path,"week_",as.character(timepoint), sep = "")
    
    #Checking if dir already exists, otherwise creates it
    existingDirCheck(path = dir)
    
    #Save distance method variable as string for titles
    if(distMethod == "bray"){
      distCharacter = "Bray-Curtis"
    }
    else if (distMethod == "wunifrac"){
      distCharacter = "Weighted Unifrac"
    }
    
    #Perform PCoA
    pcoa_results <- ordinate(ps_subset, method = "PCoA", distance = dist_subset)
    
    #Ordination plot
    p <- plot_ordination(ps_subset, pcoa_results, type = "samples", 
                         color = varToCompare) + 
      theme_classic() +
      theme(strip.background = element_blank())+
      stat_ellipse(aes(group = !!sym(varToCompare)),      # Add ellipses grouping points by genotype
                   type = "t",  # t-distribution for better fit
                   level = 0.95,  # Confidence level for the ellipse                     
                   geom = "polygon", alpha = 0)+
      labs(title = paste("PCoA of", distCharacter, "distance matrix at", timepoint, "weeks.", sep = " ")) +
      scale_color_manual(values = customColors)+
      labs(color = "Diet")+
      
      theme(
        plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
        axis.title.x = element_text(size = 12),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
        axis.text.x = element_text(size = 12),  # Adjust x-axis tick label font size
        axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
        legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
        legend.text = element_text(size = 12),  # Adjust legend font size
        panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
    
    
    test.adonis <- adonis(as.formula(paste("dist_subset ~", varToCompare)), data = data.frame(sample_data(ps_subset)))
    test.adonis <- as.data.frame(test.adonis$aov.tab)
    print(test.adonis)
    write.table(test.adonis, file = paste(dir,"/",distMethod,"_","week_",timepoint,".tsv", sep = ""), col.names = NA, row.names = TRUE, sep = '\t', quote = FALSE)
    ggsave(plot = p, filename = paste(dir,"/",distMethod,"_","week_",timepoint,".png", sep = ""), dpi = 600, height = 6, width = 6, bg = 'white')
    
    
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


#Beta diversity analysis for pairwise comparaisons, provide filtered ps object and pairs to compare
#WARNING: Formula assumes that the id of samples is named "sample_id" and is not flexible contrary to gg_group variable
betaDiversityPairwise <- function(ps, gg_group, pairs, distMethod, customColors, path){
  
  #calculating distance matrix
  dist <- phyloseq::distance(ps, method = distMethod)
  
  for(i in seq_along(pairs)){
    
    #Save pair
    pair <- unlist(pairs[i])
    
    #Save color pair
    colorPair <- unlist(customColors[i])
    
    #Creates subset of distance matrix for each selected pairs (transform dist object into matrix for which we select particular cols and rows and then transform it back to a dist object so that it can be used in the PCoA)
    dist_subset <- as.dist(as.matrix(dist)[sample_data(ps)$sample_id[sample_data(ps)[[gg_group]] %in% pair & !is.na(sample_data(ps)[[gg_group]])],sample_data(ps_samuel)$sample_id[sample_data(ps)[[gg_group]] %in% pair & !is.na(sample_data(ps)[[gg_group]])]])
    
    #Creates subset of phyloseq object for selected pairs
    ps_subset <- prune_samples(sample_data(ps)[[gg_group]] %in% pair, ps)
    
    #define name for folders 
    vs <- paste(clean_string(pair[1]), "vs", clean_string(pair[2]), sep="_")
    
    #Define dir path where graphs and stats are gonna be saved
    dir <- paste(path ,vs, sep = "")
    
    #Checking if dir already exists, otherwise creates it
    existingDirCheck(path = dir)
    
    #Save distance method variable as string for titles
    if(distMethod == "bray"){
      distCharacter = "Bray-Curtis"
    }
    else if (distMethod == "wunifrac"){
      distCharacter = "Weighted Unifrac"
    }
    
    #Perform PCoA
    pcoa_results <- ordinate(ps_subset, method = "PCoA", distance = dist_subset)
    
    #Ordination plot
    p <- plot_ordination(ps_subset, pcoa_results, type = "samples", 
                         color = gg_group) + 
      theme_classic() +
      theme(strip.background = element_blank())+
      stat_ellipse(aes(group = !!sym(gg_group)),      # Add ellipses grouping points by genotype
                   type = "t",  # t-distribution for better fit
                   level = 0.95,  # Confidence level for the ellipse                     
                   geom = "polygon", alpha = 0)+
      labs(title = paste("PCoA of", distCharacter, "distance matrix.", sep = " ")) +
      scale_color_manual(values = colorPair)+
      labs(color = "Groups")+
      
      theme(
        plot.title = element_text(size = 6, face = "bold"),  # Adjust title font size and style
        axis.title.x = element_text(size = 12),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
        axis.text.x = element_text(size = 12),  # Adjust x-axis tick label font size
        axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
        legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
        legend.text = element_text(size = 12),  # Adjust legend font size
        panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
    
    
    test.adonis <- adonis(as.formula(paste("dist_subset ~", gg_group)), data = data.frame(sample_data(ps_subset)))
    test.adonis <- as.data.frame(test.adonis$aov.tab)
    print(test.adonis)
    write.table(test.adonis, file = paste(dir,"/",distMethod,"_",vs,".tsv", sep = ""), col.names = NA, row.names = TRUE, sep = '\t', quote = FALSE)
    ggsave(plot = p, filename = paste(dir,"/",distMethod,"_",vs,".png", sep = ""), dpi = 300, height = 6, width = 6, bg = 'white')
    
    
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
