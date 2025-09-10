library(ggplot2)
library(phyloseq)
library(dplyr)

# requires a ps object, with metadata present, gg_group annotated as a factor with right group order. Requires path where graph is saved
alphaDiversityGgGroup <- function(ps, path, group){
  
  #Save name of path where graphs will be saved
  dir_path <- paste(path,"alpha_diversity", sep = "")
  
  #Check if directory exists, creates it if not
  existingDirCheck(dir_path)
  
  #Save alpha diversity measures variable to iterate through it
  measures <- c("Shannon", "Simpson","InvSimpson","Chao1","Observed")
  
  #Estinate richness measures for dataset
  richness_data <- estimate_richness(ps, measures = c("Shannon", "Simpson","InvSimpson","Chao1","Observed"))
  
  #Add sample metadata to richness dataframe
  richness_data <- cbind(as.data.frame(sample_data(ps)), richness_data)
  
  #Save richness data as csv
  write.csv(richness_data, paste(dir_path,"/richness_data.csv", sep = ""), col.names = TRUE, row.names = FALSE)

  #Save group names to create loop of pairwise comparaisons while keeping order (levels)
  groups <- levels(sample_data(ps)[[group]])
  
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      
      #Save pair variable
      pair <- paste(groups[i],"_vs_",groups[j], sep = "")
      
      #Creates path for pairwise comparaisons
      pair_path <- paste(dir_path,"/",pair, sep = "")
      
      #Check if directory exists, creates it if not
      existingDirCheck(pair_path)
      
      #Create subset of richness_data specific to pairwise comparaison
      subset <- richness_data[richness_data[[group]] %in% c(groups[i], groups[j]),]
      
      #Iterate through each alpha diversity value
      for(measure in measures){
      
        p <- ggplot(subset, aes(x = .data[[group]], y = .data[[measure]], color = .data[[group]])) +
          geom_boxplot()+
          
          # Error bars
          stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
                       aes(color = .data[[group]]),
                       width = 0.2, size = 0.7) +
          
          labs(title = paste(pair,measure), y = measure, color = "Group") +
          
          theme_minimal()+
          
          theme(
            plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
            axis.title.x = element_blank(),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
            axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust x-axis tick label font size
            axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
            legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
            legend.text = element_text(size = 12),  # Adjust legend font size
            panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
            panel.grid.minor = element_blank(),  # Remove minor grid lines
            axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
        
        
        #Save graph in appropriate folder
        ggsave(file = paste(pair_path,"/",pair,"_",measure,".png", sep = ""), plot = p, dpi = 300, bg = "white")
      }
      
    }}
}

# requires a ps object, with metadata present, gg_group annotated as a factor with right group order. Requires path where graph is saved
alphaDiversityTimeSeries<- function(ps, path, time, group){
  
  #Save name of path where graphs will be saved
  dir_path <- paste(path,"alpha_diversity", sep = "")
  
  #Check if directory exists, creates it if not
  existingDirCheck(dir_path)
  
  #Save alpha diversity measures variable to iterate through it
  measures <- c("Shannon", "Simpson","InvSimpson","Chao1","Observed")
  
  #Estinate richness measures for dataset
  richness_data <- estimate_richness(ps, measures = c("Shannon", "Simpson","InvSimpson","Chao1","Observed"))
  
  #Add sample metadata to richness dataframe
  richness_data <- cbind(as.data.frame(sample_data(ps)), richness_data)
  
  #Save richness data as csv
  write.csv(richness_data, paste(dir_path,"/richness_data.csv", sep = ""), col.names = TRUE, row.names = FALSE)
  
  #Save group names to create loop of pairwise comparaisons while keeping order (levels)
  groups <- levels(sample_data(ps)$gg_group)
  
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      
      #Save pair variable
      pair <- paste(groups[i],"_vs_",groups[j], sep = "")
      
      #Creates path for pairwise comparaisons
      pair_path <- paste(dir_path,"/",pair, sep = "")
      
      #Check if directory exists, creates it if not
      existingDirCheck(pair_path)
      
      #Create subset of richness_data specific to pairwise comparaison
      subset <- richness_data[richness_data$gg_group %in% c(groups[i], groups[j]),]
      
      #Iterate through each alpha diversity value
      for(measure in measures){
        
        p <- ggplot(subset, aes(x = gg_group, y = subset[[measure]], color = gg_group)) +
          geom_boxplot()+
          
          # Error bars
          stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
                       aes(color = gg_group),
                       width = 0.2, size = 0.7) +
          
          labs(title = paste(pair,measure), y = measure, color = "Group") +
          
          theme_minimal()+
          
          theme(
            plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
            axis.title.x = element_blank(),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
            axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust x-axis tick label font size
            axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
            legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
            legend.text = element_text(size = 12),  # Adjust legend font size
            panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add major grid lines
            panel.grid.minor = element_blank(),  # Remove minor grid lines
            axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
        
        
        #Save graph in appropriate folder
        ggsave(file = paste(pair_path,"/",pair,"_",measure,".png", sep = ""), plot = p, dpi = 300, bg = "white")
      }
      
    }}
}

# Plots alpha diversity but for all groups combinaisons
alphaDiversityGgGroup2 <- function(ps, path, gg_group, statPairs = NA, customColors){ #var1, var2
  
  #Save name of path where graphs will be saved
  dir <- paste(path,"alpha_diversity", sep = "")
  
  #Check if directory exists, creates it if not
  existingDirCheck(dir)
  
  #Save alpha diversity measures variable to iterate through it
  measures <- c("Shannon", "Simpson","InvSimpson","Chao1","Observed")
  
  #Estinate richness measures for dataset
  richness_data <- estimate_richness(ps, measures = c("Shannon", "Simpson","InvSimpson","Chao1","Observed"))
  
  #Add sample metadata to richness dataframe
  richness_data <- cbind(as.data.frame(sample_data(ps)), richness_data)
  
  #Save richness data as xlsx
  write.xlsx(richness_data, paste(dir,"/richness_data.xlsx", sep = ""))
  
  # #Save group names to create loop of pairwise comparaisons while keeping order (levels)
  # groups <- levels(sample_data(ps)[[group]])
  # 
  # for (i in 1:(length(groups)-1)) {
  #   for (j in (i+1):length(groups)) {
  #     
  #     #Save pair variable
  #     pair <- paste(groups[i],"_vs_",groups[j], sep = "")
  #     
  #     #Creates path for pairwise comparaisons
  #     pair_path <- paste(dir_path,"/",pair, sep = "")
  #     
  #     #Check if directory exists, creates it if not
  #     existingDirCheck(pair_path)
  #     
  #     #Create subset of richness_data specific to pairwise comparaison
  #     subset <- richness_data[richness_data[[group]] %in% c(groups[i], groups[j]),]
  
  
  # #Stats
  # for(i in seq_along(statPairs)){
  #   pair = unlist(statPairs[i])
  #   
  #   #Create subset of richness_data specific to pairwise comparaison
  #   subset <- richness_data[richness_data[[gg_group]] %in% c(pair[1], groups[2]),]
  #   
  #   #anova test specific to subset
  #   
  #   
  #   
  # }
  #     
  
  # Ensure that 'genotype' and 'treatment' are factors
  richness_data$genotype <- as.factor(richness_data$genotype)
  richness_data$treatment <- as.factor(richness_data$treatment)

  
  
      #Iterate through each alpha diversity value
      for(measure in measures){
        
        #stats
        # Fit the ANOVA model
        aov_model <- aov(richness_data[[measure]] ~ genotype + treatment + genotype:treatment, data = richness_data)
        
        # Summarize the results
        # print(summary(aov_model))
        
        # Perform post-hoc test (Tukey's Honest Significant Difference)
        post_hoc <- TukeyHSD(aov_model, which = "genotype:treatment")
        
        # View the results of the post-hoc test
        print(post_hoc[["genotype:treatment"]][,"p adj"])
        
        
        # plot(post_hoc)
        # 
        
        
        p <- ggplot(richness_data, aes(x = .data[[gg_group]], y = .data[[measure]], fill = .data[[gg_group]])) +
          # Error bars
          stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
                       aes(color = .data[[gg_group]]),
                       width = 0.2, size = 0.7) +
          
          geom_boxplot(
            color = customColors,
            fill = scales::alpha(customColors, 0.3)
          )+
          

          
          labs(title = paste(measure, "diversity", sep = " "), y = measure, color = "Group", fill = "Group") +
          scale_color_manual(values = customColors)+
          scale_fill_manual(values = customColors)+
          ylim(0,NA)+
          theme_minimal()+
          
          theme(
            plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
            axis.title.x = element_blank(),  # Adjust x-axis label font size and style          axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
            axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust x-axis tick label font size
            axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
            legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
            legend.text = element_text(size = 12),  # Adjust legend font size
            panel.grid.major = element_blank(),  # Add major grid lines
            panel.grid.minor = element_blank(),  # Remove minor grid lines
            axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bars
        
        
        #Save graph in appropriate folder
        ggsave(file = paste(dir,"/",measure,".png", sep = ""), plot = p, dpi = 300, bg = "white", height = 6, width =6)
      }
      
}

# requires a ps object, with metadata present, gg_group annotated as a factor with right group order. Requires path where graph is saved
alphaDiversityTimeSeries2<- function(ps, path, time, group, writeData = TRUE, saveFig = FALSE){
  
  # Empty list for graphs and data
  graphs <- list()
  
  #Save name of path where graphs will be saved
  dir_path <- paste(path,"alpha_diversity", sep = "")
  
  #Check if directory exists, creates it if not
  existingDirCheck(dir_path)
  
  #Save alpha diversity measures variable to iterate through it
  # measures <- c("Shannon", "Simpson","InvSimpson","Chao1","Observed")
  measures <- c("Chao1", "Shannon", "InvSimpson")
  
  #Estinate richness measures for dataset
  richness_data <- estimate_richness(ps, measures = c("Shannon", "Simpson","InvSimpson","Chao1","Observed"))
  
  #Add sample metadata to richness dataframe
  richness_data <- cbind(as.data.frame(sample_data(ps)), richness_data)
  
  # Write data as excel file
  if(writeData){
    write_xlsx(richness_data, path = paste0(dir_path, "/alpha_diversity_data.xlsx"))
  }
  
  # Loop through alpha diversity measures
  for(measure in measures){
    p = ggplot(richness_data, aes(x = .data[[time]], y = .data[[measure]], fill = .data[[group]])) + # , pattern = .data[["treatment"]]
      
      geom_boxplot(position = position_dodge(width = 0.8), width = 0.6,
                           color = "black") +
      labs(title = measure, y = measure, fill = "Group") # , pattern = NA

    
    # Append graph to graph list
    graphs <- append(graphs, list(p))

    # Save graph in appropriate folder
    if(saveFig){
      ggsave(file = paste(dir_path,"/",measure,".png", sep = ""), plot = p, dpi = 300, bg = "white")
      
    }
  }
  return(graphs)
}

# Plots alpha diversity as timeline, time must numeric and group must be a factor
alphaDiversityTimeline <- function(ps, time, group, shape, custom_colors, semRibbons = TRUE){
  
  # measures <- c("Shannon", "Simpson","InvSimpson","Chao1","Observed")
  measures <- c("Chao1", "Shannon", "InvSimpson", "Observed")
  
  #Estinate richness measures for dataset
  alpha_data <- estimate_richness(ps, measures = measures)
  
  #Add sample metadata to richness dataframe
  alpha_data <- cbind(as.data.frame(sample_data(ps)), alpha_data)
  
  # Empty list for graphs and data
  graphs <- list()
  
  for(measure in measures){
    
    p <- ggplot(data = alpha_data, aes(x = .data[[time]], y = .data[[measure]], group = .data[[group]], fill = .data[[group]], color = .data[[group]], shape = .data[[shape]])) +
      stat_summary(fun = mean, geom = "line", linewidth = 0.6) +
      stat_summary(fun = mean, geom = "point", size = 1.5, color = "black" ) +
      scale_color_manual(values = custom_colors)+
      scale_fill_manual(values = custom_colors)+
      scale_shape_manual(values = c(21,22))+
      labs(x = "Time", y = measure, title = measure)+
      my_theme()
    
    if(semRibbons){
      p <- p + stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.3, color = NA, aes(fill = .data[[group]]))
    }else{
      p <- p+ stat_summary(
        fun.data = mean_se,
        geom = "errorbar",
        width = 0.2,
        aes(color = .data[[group]])
      )
    }
    
    # Append graph to graph list
    graphs <- append(graphs, list(p))
    
  }
  return(graphs)
}

# Plots alpha diversity as boxplots for each group, for a single timepoint
alphaDiversityBoxplot <- function(ps, group, shape, custom_colors){
  
  # measures <- c("Shannon", "Simpson","InvSimpson","Chao1","Observed")
  measures <- c("Chao1", "Shannon", "InvSimpson")
  
  #Estinate richness measures for dataset
  alpha_data <- estimate_richness(ps, measures = measures)
  
  #Add sample metadata to richness dataframe
  alpha_data <- cbind(as.data.frame(sample_data(ps)), alpha_data)
  
  # Empty list for graphs and data
  graphs <- list()
  
  for(measure in measures){
    
    p <- ggplot(data = alpha_data, aes(x = !!sym(group), y = !!sym(measure), fill = !!sym(group), color = !!sym(group), shape = !!sym(shape))) +
      geom_boxplot(inherit.aes = FALSE,
                   mapping = aes(fill = !!sym(group), x = !!sym(group), y = !!sym(measure)),
                   colour = "black", # crans
                   alpha = 0.7)+ # transparence
      
      # stat_summary(fun = mean, geom = "", linewidth = 0.6) +
      # stat_summary(fun = mean, geom = "point", size = 1.5, color = "black" ) +
      scale_color_manual(values = custom_colors)+
      scale_fill_manual(values = custom_colors)+
      # scale_shape_manual(values = c(21,22))+
      labs(x = "", y = "Index", title = measure, fill = NULL)+
      guides(fill = "none", shape = "none")+
      my_theme()
    
    # Append graph to graph list
    graphs <- append(graphs, list(p))
    
  }
  return(graphs)
}
