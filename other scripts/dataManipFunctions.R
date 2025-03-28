#Loading functions for data manipulation (functions can be used separately for different data format)
{
  #Function removing rows with empty values (dead mice and useless stuff)
  emptyRow <- function(df){
    rows_to_remove <- c()  # Initialize an empty vector to store row indices to remove
    for(i in 1:nrow(df)){
      for(k in 1:ncol(df)){
        if(is.na(df[i,k]) || (is.character(df[i,k]) && df[i,k] == "")){ #check if they are empty cells or empty chars
          rows_to_remove <- c(rows_to_remove, i)  # Add the index of the row to remove
          break  # Exit the inner loop once an empty cell is found in the row
        }
      }
    }
    if(is.null(rows_to_remove)==FALSE){ #rows_to_remove might be empty, need to take care of that
      df <- df[-rows_to_remove, ]  # Remove the rows with empty cells
    } 
    return(df)  # Return the modified dataframe
  }
  
  #changing colnames with dates
  changeDateCols <- function(df){
    startDateCol <- grep("202",colnames(df))[1] #startDateCol corresponds to the position where cols with dates as names start
    for(i in startDateCol:ncol(df)){
      colnames(df)[i] <- gsub("\\.", "-",colnames(df)[i]) #replacing . with -
      colnames(df)[i]<- substring(colnames(df)[i], first = 2) #getting rid of X at the start
      #colnames(weight_measure_file)[i] <- substr(colnames(weight_measure_file)[i], start = 1, stop = nchar(colnames(weight_measure_file)[i]) - 1)
      
    }
    return(df)
  }
  
  #function that does pivot_longer from cols corresponding to weight measurement on a particular date
  #long tidy format
  pivotLongerWeightDate <- function(df){
    startDate <- grep("202",colnames(df))[1] #finding at which position the weight measurements date cols start appearing
    return(as.data.frame(pivot_longer(df, c(startDate:ncol(df)), cols_vary = "slowest", names_to = "date", values_to = "weight")))
  }
  
  #setting the weight column data into numeric values
  charToNum <- function(df){
    for(i in 1:ncol(df)){
      if(any(grepl(",",df[,i]))){
        df[,i] <- as.numeric(gsub("\\,", ".", df[,i]))} #checking if character is a number
    }
    return(df)
  }
  
  #convert cols containing only decimal, negative values etc into numeric
  charToNumRegEx <- function(df, colsToIgnore) {
    for (i in (colsToIgnore+1):ncol(df)) {
    # Check if the column is character
      if (is.character(df[, i])) {
        # Convert values to numeric if they match the regular expression pattern
        if (all(grepl("^-?\\d*\\.?\\d+$", df[, i]))) {
          df[, i] <- as.numeric(df[, i])
        }
      }
    }
    return(df)
  }
  
  #specific to DSS format
  commaInsteadPoint <- function(df){
    for(i in 1:ncol(df)){
      if(any(grepl(",",df[,i]))){
        for(n in 1:nrow(df)){
          if((grepl("[A-Za-z]", df[n,i]))==FALSE){
            df[n,i] <- gsub("\\,", ".", df[n,i])} #checking if character is a number
        }
        
      } 
    }
    return(df)
  }

  #Transforming dates into numeric format for statistical measurements
  dateToNum<- function(df){
    for(i in 1:ncol(df)){
      if(class(df[,i])=="Date"){ #check if dataframe col is of type "Date"
        referenceDate <- min(df[,i]) #choose a the earliest date as ref
        df$time_numeric <- as.numeric(df[,i] - referenceDate) #substraction by reference and as_numeric for every date
      } 
    } 
    return(df)
  }
  
  #function producing col with combination of treatment and diet
  #colPlacement is where we want to place the new column 
  ggGrouping <- function(df, colPlacement) {
    gg_group <- paste(df$diet, df$treatment)
    df <- cbind(df[, 1:(colPlacement - 1)], gg_group, df[, colPlacement:ncol(df)])
    return(df)
  }
  
  #defining a function that takes a dataframe, the string we are looking for in the first row
  #and the columns we wanna ignore and add at the end (such as diet, cage etc)
  #probably could be optimized by finding the patterns inside the function but works for now
  colStringSelect <- function(df,string,cols_ignored){
    selected_cols <- grep(string, df[1,(cols_ignored+1):ncol(df)])
    selected_cols <- selected_cols + cols_ignored
    df_select <- df[, selected_cols, drop = FALSE]
    result_df <- cbind(df[,1:cols_ignored],df_select)
    return(as.data.frame(result_df))
  }
  
  #Function that replaces colnames with dates, according to DSS length and day starting
  dayColnames <- function(df,dssLength,dssStart,cols_ignored){
    dates <- c() #creating vector
    dates <- c(dates, as.numeric(as.Date(dssStart))) #putting date of start of DSS in numeric
    for(i in 1:dssLength+1){
      dates <-c(dates, dates[i-1]+1)
    }
    colnames(df)[(cols_ignored + 1):(cols_ignored + dssLength + 1)] <- format(as.Date(dates),"%Y-%m-%d")
    return(df)
  }
  
  #Calculates the percentage change in weight according to first weight measure in the dataframe
  #((Day X)/(Day 1))×100 = important for disease index calcuation
  percentageWeightChange <- function(df) {
    weightColPositions <- grep("weight", df[1,])
    #adjusting weightColPositons values according to shift caused by adding new cols
    adjustedColPositions <- c()
    for(i in 1:length(weightColPositions)){
      adjustedColPositions <- c(adjustedColPositions, (weightColPositions[i]+i-1))
    }
    
    # Loop over weight columns
    for (i in adjustedColPositions) {
      # Calculate percentage weight change based on day 0 weight measure
      percentage_change <- (((as.numeric(df[, i])/as.numeric(df[, weightColPositions[1]])) * 100)-100)
      
      # Add the calculated percentage weight change as a new column next to the current weight column
      new_col_index <- i
      new_col_name <- paste0(colnames(df)[i], "_w_change")
      df <- cbind(df[, 1:i], new_column = percentage_change, df[, (i+1):ncol(df)])
      colnames(df)[new_col_index+1] <- new_col_name
      df[1,i+1] <- "w_change" #first row filled with this information = useful for creating the indiv dataframes and merging them
    }
    
    return(df)
  }
  
  # Function to calculate weight percentage change
  calculate_weight_percentage <- function(df, fromDay0) {
    df$weight_pct_change <- 100  # Initialize PercentageChange column
    
    nmice <- length(unique(df$id)) #Count number of mice
    
    if(fromDay0) {
      initial_weights <- df$weight[1:nmice]
      k = 1
      for (i in (nmice + 1):nrow(df)) {
        current_weight <- df$weight[i]
        
        if(i %% 48 != 0){
          k <- k+1
          }
        
        # Calculate percentage change and update the weight_pct_change column
        df$weight_pct_change[i] <- (((current_weight - initial_weights[i-nmice*k]) / initial_weights[i-nmice*k]) * 100)+100
      }
    }
    else {
      # Calculate percentage change for each time point relative to the previous time point
      for (i in (nmice + 1):nrow(df)) {
        current_weight <- df$weight[i]
        prev_weight <- df$weight[i - nmice]
        
        # Calculate percentage change and update the weight_pct_change column
        df$weight_pct_change[i] <- (((current_weight - prev_weight) / prev_weight) * 100)+100
      }
    }
    
    return(df)
  }
  
  read_excel_allsheets <- function(filename, tibble = FALSE) {
    # I prefer straight data.frames
    # but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
  }

  
  
  
  
  #Function that calculates disease index for DSS
  #	Weight Variation: 0 - None, 1 - 1%-5%, 2 - 5%-10%, 3 - 10%-20%, 4 - >20%
  #Stool Consistency: 0 - Normal, 2 - Loose, 4 – Diarrhea
  #Fecal Bleeding: Fecal Bleeding: 0 - None, 2 - Hemoccult Positive, 4 - Gross rectal bleeding
  dssDiseaseIndexCalculation <- function(df,negativeOnly){ #negativeOnly TRUE if wanna test only negative weight variation, FALSE otherwise
    df$index <- 0 #initializing disease index col
    for(i in 1:nrow(df)){
      #consistency, NL are considered as L and LD are considered as D
      if(grepl("D",df$consistency[i])){
        df$index[i] <- df$index[i]+4
      }else if(grepl("L",df$consistency[i])){
        df$index[i] <- df$index[i]+2
        
      }else if(grepl("NL",df$consistency[i])){
        df$index[i] <- df$index[i]+1
        
      }
      else{
        df$index[i] <- df$index[i]+0
      }
      #hemoccult
      if(grepl("G",df$hemo[i])){
        df$index[i] <- df$index[i]+4
      }else if(grepl("Yes",df$hemo[i])){
        df$index[i] <- df$index[i]+2
      }else{
        df$index[i] <- df$index[i]+0
      }
      if(isFALSE(negativeOnly)){
        #weight variation, we use absolute percent variation here
        if (abs(df$weight_change_percent[i]) >= 0 & abs(df$weight_change_percent[i]) < 1) {
          df$index[i] <- df$index[i] + 0
        } else if (abs(df$weight_change_percent[i]) >= 1 & abs(df$weight_change_percent[i]) < 5) {
          df$index[i] <- df$index[i] + 1
        } else if (abs(df$weight_change_percent[i]) >= 5 & abs(df$weight_change_percent[i]) < 10) {
          df$index[i] <- df$index[i] + 2
        } else if (abs(df$weight_change_percent[i]) >= 10 & abs(df$weight_change_percent[i]) < 20) {
          df$index[i] <- df$index[i] + 3
        } else if (abs(df$weight_change_percent[i]) >= 20) {
          df$index[i] <- df$index[i] + 4
        } else if(negativeOnly){
          # Weight variation, considering only negative percent variation
          if (df$weight_change_percent[i] < 0 & abs(df$weight_change_percent[i]) >= 0 & abs(df$weight_change_percent[i]) < 1) {
            df$index[i] <- df$index[i] + 0
          } else if (df$weight_change_percent[i] < 0 & abs(df$weight_change_percent[i]) >= 1 & abs(df$weight_change_percent[i]) < 5) {
            df$index[i] <- df$index[i] + 1
          } else if (df$weight_change_percent[i] < 0 & abs(df$weight_change_percent[i]) >= 5 & abs(df$weight_change_percent[i]) < 10) {
            df$index[i] <- df$index[i] + 2
          } else if (df$weight_change_percent[i] < 0 & abs(df$weight_change_percent[i]) >= 10 & abs(df$weight_change_percent[i]) < 20) {
            df$index[i] <- df$index[i] + 3
          } else if (df$weight_change_percent[i] < 0 & abs(df$weight_change_percent[i]) >= 20) {
            df$index[i] <- df$index[i] + 4
          }
          
        }
      }
      
    }
    return(df)
  }
}




# Custom function to calculate standard error with optional scaling
mean_cl_normal <- function(x, mult = 2.201) {
  mean_val <- mean(x)
  se_val <- sd(x) / sqrt(length(x))
  ymin_val <- pmax(mean_val - mult * se_val, 0)
  data.frame(y = mean_val, ymin = ymin_val, ymax = mean_val + mult * se_val)
}

#function that is combination of functions above, startDateCol is position at which date cols start appearing
weightDataManipulation<- function(df,groupInfoCols, fromDay0){
  df <- emptyRow(df)
  df <- changeDateCols(df)
  df <- pivotLongerWeightDate(df)
  
  #putting "date" col as a date variable type (not character)
  df$date <- as.Date(df$date)
  
  #setting the diet column data into a string format so that it can be used into ggplot (if its numeric)
  if(is.numeric(df$diet)){
    df$diet <- as.character(df$diet)
  }
  
  df <- charToNum(df)
  
  df <- dateToNum(df)
  df <- ggGrouping(df, colPlacement = (groupInfoCols+1))
  df <- calculate_weight_percentage(df, if(fromDay0 == TRUE){fromDay0 = TRUE}else{fromDay0 = FALSE})
  
  return(df)
}











#function for manipulating dss follow up sheet data
dssFollowupManipulation <- function(df, groupInfoCols, dateStart, nbrDays, negativeOnly){
  
  #getting rid of empty rows if they are (dead mice or issues)
  df <- emptyRow(df)
  
  #replacing "," by "."
  df <- commaInsteadPoint(df)
  
  #adding weight_change percentages cols
  df <- percentageWeightChange(df)

  
  #creating 3 different dataframes for each metric (weight, hemoccult, stool consistency)
  dss_weight <- colStringSelect(df,"weight",cols_ignored = groupInfoCols) #everytime we skip the first 4 cols (diet, treatment etc)
  dss_weight_change <- colStringSelect(df,"w_change",cols_ignored = groupInfoCols)
  dss_hemo <- colStringSelect(df,"bleeding",cols_ignored = groupInfoCols)
  dss_consistency <- colStringSelect(df,"consistency",cols_ignored = groupInfoCols)
  
  print("3 df created")
  
  #replacing date colnames for each individual df (Day0 replaced by "2024-xx-xx" and so on)
  dss_weight <- dayColnames(dss_weight,nbrDays,dateStart,cols_ignored = groupInfoCols)
  dss_weight_change <- dayColnames(dss_weight_change,nbrDays,dateStart,cols_ignored = groupInfoCols)
  dss_hemo <- dayColnames(dss_hemo,nbrDays,dateStart,cols_ignored = groupInfoCols)
  dss_consistency <- dayColnames(dss_consistency,nbrDays,dateStart,cols_ignored = groupInfoCols)
  
  print("colnames replacement")
  
  #removing first row (useless information)
  dss_weight <- dss_weight[-1,]
  dss_weight_change <- dss_weight_change[-1,]
  dss_hemo <- dss_hemo[-1,]
  dss_consistency <- dss_consistency[-1,]
  
  print("1st row removal")
  
  #using pivot_longer to put data into long tidy format (makes it usable by ggplot)
  groupInfoCols <- groupInfoCols+1 #date col is added as information col, groupInfoCols increases then by one
  dss_weight <- pivot_longer(dss_weight, c(groupInfoCols:length(dss_weight)), cols_vary = "slowest", names_to = "date", values_to = "weight")
  dss_weight_change <- pivot_longer(dss_weight_change, c(groupInfoCols:length(dss_weight_change)), cols_vary = "slowest", names_to = "date", values_to = "weight_change_percent")
  dss_hemo <- pivot_longer(dss_hemo, c(groupInfoCols:length(dss_hemo)), cols_vary = "slowest", names_to = "date", values_to = "hemo")
  dss_consistency <- pivot_longer(dss_consistency, c(groupInfoCols:length(dss_consistency)), cols_vary = "slowest", names_to = "date", values_to = "consistency")
  
  print("pivot longer")
  
  #putting inside this variable colnames from a individual dataframe for later
  df_colnames <- colnames(dss_weight)[0:groupInfoCols]
  
  #merging the data
  combined_df <- bind_cols(dss_weight, dss_weight_change, dss_hemo, dss_consistency)
  
  print("combined df created")
  
  #Identify duplicated columns
  duplicated_cols <- duplicated(names(combined_df)) | duplicated(t(combined_df))
  
  #Remove duplicated columns
  combined_df <- combined_df[, !duplicated_cols, drop = FALSE]
  
  print("removal of duplicates")
  
  #changing colnames because they are messy after the merge
  colnames(combined_df)[0:4] <- df_colnames
  colnames(combined_df)[groupInfoCols] <- "date"
  
  print("colnames modification")
  
  #adding a gg_grouping variable (combination of treatment and diet)
  combined_df <- ggGrouping(combined_df,colPlacement = (groupInfoCols+1))
  groupInfoCols <- groupInfoCols+1 #adding gg_group, new info col

  #transforming measure such as weight measures in numeric variables
  combined_df <- charToNumRegEx(combined_df,groupInfoCols)

  #putting "date" col as a date variable type (not character)
  combined_df$date <- as.Date(combined_df$date)

  #adding col with date as numerical values
  combined_df <- dateToNum(combined_df)

  #calculating disease index
  combined_df <- dssDiseaseIndexCalculation(combined_df, negativeOnly = negativeOnly)
  
  return(combined_df)
}


#function for data manipulation of the dissection measures (body weight, liver, spleen and colon length)
dissectionDataManipulation <- function(df,groupInfoCols, numerical = TRUE){
  
  #getting rid of empty rows (dead mice)
  df <- emptyRow(df)
  
  #measures cols into numerical variables
  if(isFALSE(numerical)){
    df <- charToNum(df)
  }
  

  
  
  #setting the diet column data into a string format so that it can be used into ggplot (if its numeric)
  if(is.numeric(df$diet)){
    df$diet <- as.character(df$diet)
  }
  
  #creating groups which are combinations of treatment and diet
  df <- ggGrouping(df, colPlacement = (groupInfoCols+1))
  
  #dividing these measures by final weight of the mice (normalization)
  df$std_spleen_weigth <- df$spleen_weight/df$body_weight
  df$std_liver_weigth <- df$liver_weight/df$body_weight
  df$std_colon_len <- df$colon_length/df$body_weight
  
  return(df)
}


#order in which categorical groups are displaced
desired_order <- c("50 water", "500 water", "50 dss", "500 dss")

# Define your custom color palette
custom_colors <- c("blue","red")




#function for plotting DSS disease index data
dssDiseaseIndexPlot <- function(df){
 plot <- df %>%
    ggplot(aes(x = time_numeric, y = index, color = diet))+
    stat_summary(aes(group = gg_group, shape = treatment),fun ="mean", geom = "point", size = 3)+
    stat_summary(aes(group = gg_group), fun = "mean", geom = "line", size = 1)+
    stat_summary(aes(color = diet), fun.data = mean_cl_normal, geom="errorbar", width=0.1, size = 1) + #adding SEM error bars
   
    labs(title = "Disease activity index (DAI) during DSS exposure",
         x = "Day",
         y = "DAI",
         color = "Diet")+
   scale_color_manual(values = custom_colors)+
    guides(shape = 'none')+
    theme_minimal()+
    theme(
      plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
      axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
      axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
      axis.text.x = element_text(size = 12),  # Adjust x-axis tick label font size
      axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
      legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
      legend.text = element_text(size = 12),  # Adjust legend font size
      panel.grid.major = element_blank(),  # Add major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black", size = 1)  # Include axis lines  # Include axis bars
    )+
   ylim(0, max(df$index))
  return(plot)
}


#plotting DSI at final day of DSS period
dssDsiFinalDay <- function(df){
  
  #order in which categorical groups are displaced
  desired_order1 <- c("50 dss", "500 dss")
  
  #subsetting for only the final day 
  df <- subset(df, time_numeric == 5)
  
  #Statistical measurements
  #Subsetting the dataframe to include only the DSS groups
  data <- df[df$treatment == "dss",]
  group1 <- data[data$diet == "50", ]
  group2 <- data[data$diet == "500", ]
  
  # Shapiro-Wilk test for normality
  print(shapiro.test(group1$index))
  print(shapiro.test(group2$index))
  
  # Levene's test for homogeneity of variance
  print(leveneTest(group1$index, group2$index))
  
  good_test <- wilcox.test(index ~ gg_group, data = data, exact = FALSE)
  print(good_test)
  p_value <- good_test$p.value
  print(str(p_value))
  
  plot <- df %>%
    ggplot(aes(x = factor(gg_group, levels = desired_order1),  y = index, color = as.character(diet))) +
    
    stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..-0.25, yend=..y.., color = diet), size =1)+ #adding horizontal bars representing means
    stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..+0.25, yend=..y.., color = diet), size =1)+ 
    stat_summary(aes(color = diet), fun.data="mean_cl_normal", geom="errorbar", width=0.2, size = 1) + #adding SEM error bars
    #geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), alpha = 0.5) +
    
    geom_signif(comparisons =   desired_order1,
                map_signif_level = TRUE, # Pour afficher l'étoile
                y_position = 9, # Ajuste cette valeur selon ton graphique
                annotations = ifelse(p_value < 0.05, "*", "ns")) + # Affiche * si significatif +
    
    
    labs(title = "DSI at day 5 of DSS administration",
         x = "Treatment",
         y = "DSI",
         color = "Diet")+ 
    scale_color_discrete(labels = c("50 ppm FeSO4", "500 ppm FeSO4"))+
    scale_color_manual(values = custom_colors)+
    scale_x_discrete(labels = c("50 ppm + DSS", "500 ppm + DSS"))+
    theme_minimal() +  
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12),
      panel.grid.major =element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 1)
    ) +
    ylim(0, max(df$index))
  
  

  # # Perform t-test
  # print(t.test(index ~ diet, data = data))
  # 
  # #perform Kruskal-Wallis test
  # print(kruskal.test(index ~ diet, data = data))
  

  
  return(plot)
}

#plotting the weight measures scatter plot
weightPlot <- function(df, percentage = FALSE, diet_only = FALSE, stats = TRUE){
  
  # Custom function to calculate standard error with optional scaling
  mean_cl_normal <- function(x, mult = 1) {
    mean_val <- mean(x)
    se_val <- sd(x) / sqrt(length(x))
    ymin_val <- pmax(mean_val - mult * se_val, 0)
    data.frame(y = mean_val, ymin = ymin_val, ymax = mean_val + mult * se_val)
  }
  
  
  desired_order <- c("50 water", "500 water", "50 DSS", "500 DSS")
  
  plot <- df %>%
    ggplot(aes(x = time_numeric, y = if(percentage) { weight_pct_change } else { weight }, color = diet)) +
    stat_summary(aes(group = if(diet_only) {diet} else {gg_group}, shape = if(diet_only) {NULL} else {treatment}), fun = "mean", geom = "point", size = 3) +
    stat_summary(fun = "mean", geom = "line", aes(group = gg_group, linetype = ifelse(grepl("dss", gg_group, ignore.case = TRUE), "DSS", "Water")), size = 1) +
    stat_summary(fun.data="mean_cl_normal", geom="errorbar", color="red", width=0.5, aes(group = if(diet_only) {diet} else {gg_group})) + #adding SEM error bars
    labs(title = "Body weight change over time",
         x = "Time (days)",
         y = if(percentage) {"Weight % change"} else {"Weight (g)"},
         color = "Diet") +
    scale_linetype_manual(name = "Treatment", 
                          values = c("DSS" = "dashed", "Water" = "solid")) +
    scale_color_manual(values = custom_colors)+
    guides(shape = 'none')+
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
      axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
      axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
      axis.text.x = element_text(size = 12),  # Adjust x-axis tick label font size
      axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
      legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
      legend.text = element_text(size = 12),  # Adjust legend font size
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black", size = 1)  # Include axis lines  # Include axis bars
    )
  if(stats){
    # Initialize an empty list to store Tukey's HSD results for each time point
    tukey_results <- list()
    
    # Loop through each unique time point
    for (time_point in unique(df$time_numeric)) {
      # Subset the data for the current time point
      df_time <- subset(df, time_numeric == time_point)
      
      # Fit a simple ANOVA model for the current time point
      simple_model <- aov(weight_pct_change ~ treatment * diet, data = df_time)
      
      # Perform Tukey's HSD test and store the results in the list
      tukey_results[[as.character(time_point)]] <- TukeyHSD(simple_model)
    }
    
    # Alternatively, print all results
    for (time_point in names(tukey_results)) {
      cat("\nTukey's HSD results for Time Point:", time_point, "\n")
      print(tukey_results[[time_point]])
    }
    # #Statistical measurements
    # dss50 <- df[df$gg_group == "50 DSS",]
    # dss500 <- df[df$gg_group == "500 DSS",]
    # water50 <- df[df$gg_group == "50 water",]
    # water500 <- df[df$gg_group == "500 water",]
    # model <- aov(weight_pct_change ~ treatment * diet * time_numeric + Error(id/time_numeric), data = df)
    # print(summary(model))
    # print(TukeyHSD(model))
  }
  return(plot)
}


#function for dissection measures box_plot
dissecBoxplot <- function(df,organ, display_significance_bars, permAnova = FALSE){
  
  custom_colors <- c("blue","red","darkblue","darkred")
  
  #titles for the boxplots and Y axis labels
  nrmString <- "Normalized"
  plotTitle1 <- "weight measures"
  yAxisLabel1 <- "weight ("
  yAxisLabel2 <- "mg/g)"
  plotTitle <- paste(nrmString,organ,plotTitle1,sep = " ")
  yAxisLabel <- paste(nrmString, organ, yAxisLabel1, yAxisLabel2)
  
  
  if(organ == "body"){
    plotTitle <- "Body weight measures at final day of the experiment"
    yAxisLabel <- "Body weight (g)"
    responseVariable <- "body_weight"
  }else if(organ == "liver"){
    plotTitle <- plotTitle
    yAxisLabel <- yAxisLabel
    responseVariable <- "std_liver_weigth"
  }else if(organ == "spleen"){
    plotTitle <- plotTitle
    yAxisLabel <- yAxisLabel
    responseVariable <- "std_spleen_weigth"
  }else if(organ == "colon"){
    plotTitle <- "Colon length measures"
    yAxisLabel <- "Colon length (cm)"
    responseVariable <- "colon_length"
  }
  else if(organ == "colon_nrm"){
    plotTitle <- "Normalized colon length measures"
    yAxisLabel <- "Normalized colon length (cm/g)"
    responseVariable <- "colon_length_nrm"
  }
  
  # Plot with mean points
  plot <- ggplot(data = df, aes(x = gg_group,  y = !!sym(responseVariable), color = gg_group)) +
          geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), alpha = 0.5) +
          stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..-0.25, yend=..y..), color = "black", size = 1)+ #adding horizontal bars representing means
          stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..+0.25, yend=..y..), color = "black", size = 1)+ 
          stat_summary(fun.data="mean_cl_normal", geom="errorbar", width=0.2, size = 0.7) + #adding SEM error bars
          labs(title = plotTitle,
               x = "",
               y = yAxisLabel,
               color = "")+ 
          scale_color_manual(values = custom_colors)+
          scale_x_discrete(labels = c("50 ppm\ncontrol", "500 ppm\ncontrol", "50 ppm\nDSS", "500 ppm\nDSS"))+
          theme_minimal() +  
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black", size = 1)
          ) +
          ylim(0, max(df[[responseVariable]]))
  
  #Statistical measurements
  dss50 <- df[df$gg_group == "50 dss",]
  dss500 <- df[df$gg_group == "500 dss",]
  water50 <- df[df$gg_group == "50 water",]
  water500 <- df[df$gg_group == "500 water",]
  
  # Shapiro-Wilk test for normality
  print(shapiro.test(dss50[[responseVariable]]))
  print(shapiro.test(dss500[[responseVariable]]))
  print(shapiro.test(water50[[responseVariable]]))
  print(shapiro.test(water500[[responseVariable]]))
  
  # Levene's test for homogeneity of variance
  print(leveneTest(df[[responseVariable]] ~ gg_group, data = df))
  
  # Perform anova
  if(permAnova){
    model <- lmp(df[[responseVariable]] ~ treatment * diet, data = df)
    print(summary(model))
  }else{
    result <- aov(df[[responseVariable]] ~ treatment * diet, data = df)
    print(summary(result))
    print(TukeyHSD(result))
  }

  
  
  # print("this is non parametric")
  # kruskal_test(df[[responseVariable]] ~ treatment * diet, data=df)
  # post_hoc <- data %>%
  #   dunn_test(body_weight ~ treatment, p.adjust.method = "bonferroni")
  # 
  # print(post_hoc)
  
  
  return(plot)
}


#Function for ior iron measurements, gg_group must be a factor and ordered, diet must be as character
ironBoxplot <- function(df, measure, display_significance_bars = TRUE, title, y_axis_title, custom_colors, path){
  
  #Plot with mean points
  plot <- df %>%
    ggplot(aes(x = gg_group,  y = !!sym(measure), color = diet)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), alpha = 0.5) +
    stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..-0.25, yend=..y..), size = 1)+ #adding horizontal bars representing means
    stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..+0.25, yend=..y..), size = 1)+ 
    stat_summary(fun.data="mean_cl_normal", geom="errorbar", width=0.2, size = 0.7) + #adding SEM error bars
    labs(title = title,
         y = y_axis_title,
         x = "",
         color = "Diet")+ 
    scale_color_manual(values = custom_colors)+
    ylim(0,max(df[[measure]])+1/4*max(df[[measure]]))+
    # ylim(0,250)+ #For comparing to Claire Balb/C measurements
    theme_minimal() +  
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 1)
    )
  
  #Saving the graph
  ggsave(plot = plot, filename = paste(path, "/", title, ".png", sep = ""), height = 6, width = 5, dpi = 300, bg = "white")
  
  # #Statistical measurements
  # dss50 <- df[df$gg_group == "50 dss",]
  # dss500 <- df[df$gg_group == "500 dss",]
  # water50 <- df[df$gg_group == "50 water",]
  # water500 <- df[df$gg_group == "500 water",]
  # 
  # # Shapiro-Wilk test for normality
  # print(shapiro.test(dss50[[responseVariable]]))
  # print(shapiro.test(dss500[[responseVariable]]))
  # print(shapiro.test(water50[[responseVariable]]))
  # print(shapiro.test(water500[[responseVariable]]))
  # 
  # # Levene's test for homogeneity of variance
  # print(leveneTest(df[[responseVariable]] ~ gg_group, data = df))
  # 
  # # Perform anova
  # if(permAnova){
  #   model <- lmp(df[[responseVariable]] ~ treatment * diet, data = df)
  #   print(summary(model))
  # }else{
  #   result <- aov(df[[responseVariable]] ~ treatment * diet, data = df)
  #   print(summary(result))
  #   print(TukeyHSD(result))
  # }
  
  
  
  # print("this is non parametric")
  # kruskal_test(df[[responseVariable]] ~ treatment * diet, data=df)
  # post_hoc <- data %>%
  #   dunn_test(body_weight ~ treatment, p.adjust.method = "bonferroni")
  # 
  # print(post_hoc)
  
  
  return(plot)
}
