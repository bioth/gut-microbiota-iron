#Revised version, not doing separate distance calculations and adding path option

#Beta diversity analysis for different timepoints. You must provide a filtered ps object, the timeVariable and the varToCompare (present in sample_data)
betaDiversityTimepoint <- function(ps, timeVariable, varToCompare, distMethod, customColors, font, path){
  
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
        plot.title = element_text(size = 16, face = "bold", family = font),  # Adjust title font size and style
        axis.title.x = element_text(size = 14, face = "bold", family = font), # Adjust x-axis label font size and style   
        axis.title.y = element_text(size = 14, face = "bold", family = font), # Adjust y-axis label font size and style
        axis.text.x = element_text(size = 14, face = "bold", family = font),  # Adjust x-axis tick label font size
        axis.text.y = element_text(size = 14, face = "bold", family = font),  # Adjust y-axis tick label font size
        legend.title = element_text(size = 14, face = "bold", family = font),  # Remove legend title
        legend.text = element_text(size = 14, family = font),  # Adjust legend font size
        panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
    
    
    test.adonis <- adonis(as.formula(paste("dist_subset ~", varToCompare)), data = data.frame(sample_data(ps_subset)))
    test.adonis <- as.data.frame(test.adonis$aov.tab)
    print(test.adonis)
    
    #Save stats as an excel file
    write.xlsx(test.adonis, paste(dir,"/",distMethod,"_","week_",timepoint,".xlsx", sep = ""))
    
    #Save figure
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
betaDiversityPairwise <- function(ps, gg_group, pairs, distMethod, customColors, font, displayPValue = TRUE, dim = c(6,6), transform = FALSE, path){
  
  #Transform abundance into relative abundances or log_transformed values
  if(transform == "rel_ab"){
    ps <- transformCounts(ps, transformation = "rel_ab")
  } else if(transform == "log"){
    ps <- transformCounts(ps, transformation = "log", log_base = 10)
  }
  
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
    colnames(pcoa_results$vectors) <- gsub("Axis.", "PC", colnames(pcoa_results$vectors)) #Replace colnames "Axis.n" by "PCn"
    
    #Statistics, and save table in folder
    test.adonis <- adonis(as.formula(paste("dist_subset ~", gg_group)), data = data.frame(sample_data(ps_subset)))
    test.adonis <- as.data.frame(test.adonis$aov.tab)
    p_value <- test.adonis[["Pr(>F)"]][1]  #Extract p-value from adonis result
    print(test.adonis)
    
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
      theme(aspect.ratio = 1)+ #Scale the x and y axis the same 
      
      theme(
        plot.title = element_text(size = 12, face = "bold", family = font),  # Adjust title font size and style
        axis.title.x = element_text(size = 12, face = "bold", family = font),
        axis.title.y = element_text(size = 12, face = "bold", family = font),# Adjust y-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
        axis.text.x = element_text(size = 12, face = "bold", family = font, color = "black"),  # Adjust x-axis tick label font size
        axis.text.y = element_text(size = 12, face = "bold", family = font, color = "black"),  # Adjust y-axis tick label font size
        legend.title = element_text(size = 10, face = "bold", family = font),  # Remove legend title
        legend.text = element_text(size = 9, face = "bold", family = font),  # Adjust legend font size
        panel.grid.major = element_blank(),  # Add major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black", size = 1))+ # Include axis lines  # Include axis bars
      
      # Add the p-value under the legend
      if(displayPValue){
        annotate("text", 
                 x = Inf, y = -Inf,  # Position at bottom right, under the legend
                 label = paste("p =", round(p_value, 3)), 
                 hjust = 1, vjust = -0.5, 
                 size = 4, color = "black", 
                 fontface = "bold",
                 family = font)  # Make the entire text bold
      }

    ggsave(plot = p, filename = paste(dir,"/",distMethod,"_",vs,".png", sep = ""), dpi = 300, height = dim[1], width = dim[2], bg = 'white')
    
    #Save stats as an excel file
    write.xlsx(test.adonis, paste(dir,"/",distMethod,"_",vs,".xlsx", sep = ""))
    
    
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


#Beta diversity analysis for all groups
#WARNING: Formula assumes that the id of samples is named "sample_id" and is not flexible contrary to gg_group variable
# betaDiversityAll <- function(ps, gg_group, distMethod, customColors, font, displayPValue = TRUE, dim = c(6,6), transform = FALSE, variancePlot = FALSE, path){
#   
#   #Transform abundance into relative abundances or log_transformed values
#   if(transform == "rel_ab"){
#     ps <- transformCounts(ps, transformation = "rel_ab")
#   } else if(transform == "log"){
#     ps <- transformCounts(ps, transformation = "log", log_base = 10)
#   }
#   
#   #calculating distance matrix
#   dist <- phyloseq::distance(ps, method = distMethod)
#   
#   #Save distance method variable as string for titles
#   if(distMethod == "bray"){
#     distCharacter = "Bray-Curtis"
#   }
#   else if (distMethod == "wunifrac"){
#     distCharacter = "Weighted Unifrac"
#   }
#   
#   #Perform PCoA
#   pcoa_results <- ordinate(ps, method = "PCoA", distance = dist)
#   colnames(pcoa_results$vectors) <- gsub("Axis.", "PC", colnames(pcoa_results$vectors)) #Replace colnames "Axis.n" by "PCn"
#   
#   #Plot portion of variance explained by each PCs
#   if(variancePlot){
#     
#     # Extract eigenvalues from PCoA results
#     variance_explained <- pcoa_results$values$Relative_eig*100
#     
#     # Create a data frame for plotting
#     variance_df <- data.frame(PC = paste0("PC", 1:length(variance_explained)),
#                               VarianceExplained = variance_explained)
#     
#     #Display only for the first 5
#     variance_df <- variance_df[1:5,]
#     
#     
#     # Plot the variance explained by each principal component
#     variance_plot <- ggplot(variance_df, aes(x = PC, y = VarianceExplained)) +
#       geom_point() +
#       geom_line(group = 1)+
#       theme_minimal() +
#       labs(title = paste("Variance Explained by Principal Components\n(", distCharacter, " distance)", sep = ""),
#            x = "Principal Component", y = "Variance Explained (%)") +
#       ylim(0,NA)+
#       theme(plot.title = element_text(size = 14, face = "bold"),
#             axis.title.x = element_text(size = 12, face = "bold"),
#             axis.title.y = element_text(size = 12, face = "bold"),
#             axis.text.x = element_text(size = 12, face = "bold"),
#             axis.text.y = element_text(size = 12, face = "bold"))
#     
#     # Save the variance explained plot
#     ggsave(plot = variance_plot, filename = paste(path, distMethod, "_variance_explained.png", sep = ""), dpi = 300, height = 6, width = 6, bg = 'white')
#     
#   }
#   
#   #Statistics, and save table in folder
#   test.adonis <- adonis(as.formula(paste("dist ~", gg_group)), data = data.frame(sample_data(ps)))
#   test.adonis <- as.data.frame(test.adonis$aov.tab)
#   p_value <- test.adonis[["Pr(>F)"]][1]  #Extract p-value from adonis result
#   print(test.adonis)
#   write.table(test.adonis, file = paste(path,distMethod,"_all_groups.tsv", sep = ""), col.names = NA, row.names = TRUE, sep = '\t', quote = FALSE)
#   
#   #Ordination plot
#   p <- plot_ordination(ps, pcoa_results, type = "samples", 
#                        color = gg_group) + 
#     theme_classic() +
#     theme(strip.background = element_blank())+
#     stat_ellipse(aes(group = !!sym(gg_group)),      # Add ellipses grouping points by genotype
#                  type = "t",  # t-distribution for better fit
#                  level = 0.95,  # Confidence level for the ellipse                     
#                  geom = "polygon", alpha = 0)+
#     labs(title = paste("PCoA of", distCharacter, "distance matrix.", sep = " ")) +
#     scale_color_manual(values = customColors)+
#     labs(color = "Groups")+
#     theme(aspect.ratio = 1)+ #Scale the x and y axis the same 
#     
#     theme(
#       plot.title = element_text(size = 12, face = "bold", family = font),  # Adjust title font size and style
#       axis.title.x = element_text(size = 12, face = "bold", family = font),
#       axis.title.y = element_text(size = 12, face = "bold", family = font),# Adjust y-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
#       axis.text.x = element_text(size = 12, face = "bold", family = font, color = "black"),  # Adjust x-axis tick label font size
#       axis.text.y = element_text(size = 12, face = "bold", family = font, color = "black"),  # Adjust y-axis tick label font size
#       legend.title = element_text(size = 10, face = "bold", family = font),  # Remove legend title
#       legend.text = element_text(size = 9, face = "bold", family = font),  # Adjust legend font size
#       panel.grid.major = element_blank(),  # Add major grid lines
#       panel.grid.minor = element_blank(),  # Remove minor grid lines
#       axis.line = element_line(color = "black", size = 1))+ # Include axis lines  # Include axis bars
#     
#     # Add the p-value under the legend
#     if(displayPValue){
#       annotate("text", 
#                x = Inf, y = -Inf,  # Position at bottom right, under the legend
#                label = paste("p =", round(p_value, 3)), 
#                hjust = 1, vjust = -0.5, 
#                size = 4, color = "black", 
#                fontface = "bold",
#                family = font)  # Make the entire text bold
#     }
#   
#   ggsave(plot = p, filename = paste(path,distMethod,"_all_groups.png", sep = ""), dpi = 300, height = dim[1], width = dim[2], bg = 'white')
#   
# }


betaDiversityAll <- function(ps, gg_group, distMethod, customColors, font, displayPValue = TRUE, dim = c(6,6), transform = FALSE, variancePlot = FALSE, display3PCs = FALSE, path){
  
  # Transform abundance into relative abundances or log-transformed values
  if(transform == "rel_ab"){
    ps <- transformCounts(ps, transformation = "rel_ab")
  } else if(transform == "log"){
    ps <- transformCounts(ps, transformation = "log", log_base = 10)
  }
  
  # Calculating distance matrix
  dist <- phyloseq::distance(ps, method = distMethod)
  
  # Save distance method variable as string for titles
  if(distMethod == "bray"){
    distCharacter = "Bray-Curtis"
  }
  else if (distMethod == "wunifrac"){
    distCharacter = "Weighted Unifrac"
  }
  
  # Perform PCoA
  pcoa_results <- ordinate(ps, method = "PCoA", distance = dist)
  colnames(pcoa_results$vectors) <- gsub("Axis.", "PC", colnames(pcoa_results$vectors)) # Replace colnames "Axis.n" by "PCn"
  
  # Plot portion of variance explained by each PC
  if(variancePlot){
    
    # Extract eigenvalues from PCoA results
    variance_explained <- pcoa_results$values$Relative_eig*100
    
    # Create a data frame for plotting
    variance_df <- data.frame(PC = paste0("PC", 1:length(variance_explained)),
                              VarianceExplained = variance_explained)
    
    # Display only for the first 5
    variance_df <- variance_df[1:5,]
    
    # Plot the variance explained by each principal component
    variance_plot <- ggplot(variance_df, aes(x = PC, y = VarianceExplained)) +
      geom_point() +
      geom_line(group = 1) +
      theme_minimal() +
      labs(title = paste("Variance Explained by Principal Components\n(", distCharacter, " distance)", sep = ""),
           x = "Principal Component", y = "Variance Explained (%)") +
      ylim(0, NA) +
      theme(plot.title = element_text(size = 14, face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 12, face = "bold"))
    
    # Save the variance explained plot
    ggsave(plot = variance_plot, filename = paste(path, distMethod, "_variance_explained.png", sep = ""), dpi = 300, height = 6, width = 6, bg = 'white')
  }
  
  # Statistics, and save table in folder
  test.adonis <- adonis(as.formula(paste("dist ~", gg_group)), data = data.frame(sample_data(ps)))
  test.adonis <- as.data.frame(test.adonis$aov.tab)
  p_value <- test.adonis[["Pr(>F)"]][1]  # Extract p-value from adonis result
  print(test.adonis)
  
  # If display3PCs is TRUE, create a 3D PCoA plot using plotly
  if(display3PCs){
    # Extract the first three PCs for 3D plotting
    pcoa_df <- data.frame(pcoa_results$vectors[, c("PC1", "PC2", "PC3")])
    pcoa_df$Group <- sample_data(ps)[[gg_group]]
    
    # Create the 3D PCoA plot
    pcoa_3d_plot <- plot_ly(pcoa_df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Group, 
                            type = "scatter3d", mode = "markers", marker = list(size = 5)) %>%
      layout(title = paste("PCoA of", distCharacter, "distance matrix (3D view)"),
             scene = list(xaxis = list(title = "PC1"), 
                          yaxis = list(title = "PC2"),
                          zaxis = list(title = "PC3")))
    
    # Save the 3D plot
    htmlwidgets::saveWidget(pcoa_3d_plot, file = paste(path, distMethod, "_3D_PCoA_plot.html", sep = ""))
  } else {
    # Ordination plot (2D)
    p <- plot_ordination(ps, pcoa_results, type = "samples", 
                         color = gg_group) + 
      theme_classic() +
      theme(strip.background = element_blank()) +
      stat_ellipse(aes(group = !!sym(gg_group)),      # Add ellipses grouping points by genotype
                   type = "t",  # t-distribution for better fit
                   level = 0.95,  # Confidence level for the ellipse                     
                   geom = "polygon", alpha = 0) +
      labs(title = paste("PCoA of", distCharacter, "distance matrix.", sep = " ")) +
      scale_color_manual(values = customColors) +
      labs(color = "Groups") +
      theme(aspect.ratio = 1) + # Scale the x and y axis the same 
      theme(
        plot.title = element_text(size = 12, face = "bold", family = font),  # Adjust title font size and style
        axis.title.x = element_text(size = 12, face = "bold", family = font),
        axis.title.y = element_text(size = 12, face = "bold", family = font),
        axis.text.x = element_text(size = 12, face = "bold", family = font, color = "black"),
        axis.text.y = element_text(size = 12, face = "bold", family = font, color = "black"),
        legend.title = element_text(size = 10, face = "bold", family = font),
        legend.text = element_text(size = 9, face = "bold", family = font),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 1)
      ) +
      # Add the p-value under the legend
      if(displayPValue){
        annotate("text", 
                 x = Inf, y = -Inf,  # Position at bottom right, under the legend
                 label = paste("p =", round(p_value, 3)), 
                 hjust = 1, vjust = -0.5, 
                 size = 4, color = "black", 
                 fontface = "bold",
                 family = font)  # Make the entire text bold
      }
    
    # Save the 2D PCoA plot
    ggsave(plot = p, filename = paste(path, distMethod, "_all_groups.png", sep = ""), dpi = 300, height = dim[1], width = dim[2], bg = 'white')
    
    #Save stats as an excel file
    write.xlsx(test.adonis, paste(path, distMethod,"_all_groups.xlsx", sep = ""))
  }
}
