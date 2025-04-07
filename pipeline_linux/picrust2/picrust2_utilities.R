# Function that iterates through list of compounds/pathways, subtract KO_abundance dataframe for KOs associated, and returns multiple graphs
# ko_abundance needs to have kos as rownames
KOsToCompoundAbundanceGraphs <- function(ko_abundance, ko_annotations, metadata, group, customColors, sample_id_col, path, dim = c(6,6), additionalAes){
  
  existingDirCheck(path)
  
  # Iterate through the compounds/pathways
  for(i in 1:nrow(ko_annotations)){
    compound <- ko_annotations$`compound/pathway`[i]
    kos <- unlist(strsplit(ko_annotations$Kos[i], ";"))
    abundance <- ko_abundance[rownames(ko_abundance) %in% kos, ]
    total_ab <- as.data.frame(colSums(abundance))
    
    # # Associate gg_group with each sample_id using metadata information
    # for(i in 1:nrow(total_ab)){
    # total_ab[[group]][i] <- metadata[metadata[[sample_id_col]]==rownames(total_ab)[i],group]
    # }
    
    total_ab <- merge(total_ab, meta, by = "row.names")
    
    # Replacing a colname and adding a sample_id column
    colnames(total_ab)[2] <- compound
    
    p <- ggplot(total_ab, aes(x = .data[[group]], y = .data[[compound]], color = .data[[group]])) +
      geom_point(size = 1, 
                 position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) + 
      
      #Error bars
      stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
                   aes(color = .data[[group]]),
                   width = 0.2, size = 0.7,
                   position = position_dodge(-0.75)) +
      
      #Mean lines
      stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
                   aes(ymin = ..y.., ymax = ..y.., group = .data[[group]]),
                   color = "black", linewidth = 0.5, width = 0.5,
                   position = position_dodge(-0.75))+
      
      #Connecting mean points with lines
      # stat_summary(fun = mean, geom = "line", size = 1.2) +  # Connecting means with lines
      
      
      
      labs(title = compound, y = paste("Normalized Abundance"), color = "Group", x = "") +
      scale_color_manual(values = customColors)+
      scale_y_continuous(limits = c(0, NA))+
      
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
    
    if(isFALSE(is.null(additionalAes))){
      # p <- do.call("+", c(list(p), additionnalAes))
      
      p <- Reduce("+", c(list(p), additionalAes))
    }
    
    # Save plot
    ggsave(plot = p, filename = paste0(path, "/", clean_string(compound), "_predicted.png"), dpi = 300, width = dim[1], height = dim[2], bg = "white")
  }
}