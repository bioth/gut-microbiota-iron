library(ggplot2)
library(phyloseq)
library(dplyr)

#function to add SEM (1.96 for 95% confidence interval)
mean_cl_normal <- function(x, mult = 1.96) { #mult is 1.96 for a 95% confidence interval
  # Calculate the mean of the input vector x
  mean_val <- mean(x, na.rm = TRUE)
  
  # Calculate the standard error of the mean
  se_val <- sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
  
  # Return a data frame with the mean (y), and the lower (ymin) and upper (ymax) bounds
  data.frame(y = mean_val, ymin = mean_val - mult * se_val, ymax = mean_val + mult * se_val)
}

#Function checking if a dir exists and creating it otherwise
existingDirCheck <- function(path){
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Directory created: ", path)
  } else {
    message("Directory already exists: ", path)
  }
  
}

#requires a ps object, with metadata present, gg_group annotated as a factor with right group order. Requires path where graph is saved
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

#requires a ps object, with metadata present, gg_group annotated as a factor with right group order. Requires path where graph is saved
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