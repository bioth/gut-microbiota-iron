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
        for(n in 1:nrow(df)){
          if((grepl("[A-Za-z]", df[n,i]))==FALSE){
            df[n,i] <- as.numeric(gsub("\\,", ".", df[n,i]))} #checking if character is a number
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
  ggGrouping<- function(df){
    df$group <- paste(df$diet,df$treatment)
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
      new_col_name <- paste0(colnames(df)[i], "_weight_change")
      df <- cbind(df[, 1:i], new_column = percentage_change, df[, (i+1):ncol(df)])
      colnames(df)[new_col_index+1] <- new_col_name
    }
    
    return(df)
  }
}









#function that is combination of functions above, startDateCol is position at which date cols start appearing
weightDataManipulation<- function(df){
  df <- emptyRow(df)
  df <- changeDateCols(df)
  df <- pivotLongerWeightDate(df)
  
  #putting "date" col as a date variable type (not character)
  df$date <- as.Date(df$date)
  
  #setting the diet column data into a string format so that it can be used into ggplot (if its numeric)
  if(class(df$diet)=="numeric"){
    df$diet <- as.character(df$diet)
  }
  
  df <- charToNum(df)
  
  df <- dateToNum(df)
  df <- ggGrouping(df)
  
  return(df)
}



#function for manipulating dss follow up sheet data
dssFollowupManipulation <- function(df, groupInfoCols, dateStart, nbrDays){
  
  #getting rid of empty rows if they are (dead mice or issues)
  df <- emptyRow(df)
  
  #transforming weight measure from char to num
  df <- charToNum(df)
  
  #adding weight_change percentages cols
  df <- percentageWeightChange(df)
  
  #creating 3 different dataframes for each metric (weight, hemoccult, stool consistency)
  dss_weight <- colStringSelect(df,"weight",cols_ignored = groupInfoCols) #everytime we skip the first 4 cols (diet, treatment etc)
  dss_weight_change <- colStringSelect(df,"weight_change",cols_ignored = groupInfoCols)
  dss_hemo <- colStringSelect(df,"hemoccult",cols_ignored = groupInfoCols)
  dss_consistency <- colStringSelect(df,"consistency",cols_ignored = groupInfoCols)
  
  #replacing date colnames for each individual df (Day0 replaced by "2024-xx-xx" and so on)
  dss_weight <- dayColnames(dss_weight,nbrDays,dateStart,cols_ignored = groupInfoCols)
  dss_weight_change <- dayColnames(dss_weight_change,nbrDays,dateStart,cols_ignored = groupInfoCols)
  dss_hemo <- dayColnames(dss_hemo,nbrDays,dateStart,cols_ignored = groupInfoCols)
  dss_consistency <- dayColnames(dss_consistency,nbrDays,dateStart,cols_ignored = groupInfoCols)
  
  #removing first row (useless information)
  dss_weight <- dss_weight[-1,]
  dss_weight_change <- dss_weight_change[-1,]
  dss_hemo <- dss_hemo[-1,]
  dss_consistency <- dss_consistency[-1,]
  
  #using pivot_longer to put data into long tidy format (makes it usable by ggplot)
  dss_weight <- pivot_longer(dss_weight, c(5:length(dss_weight)), cols_vary = "slowest", names_to = "date", values_to = "weight")
  dss_weight_change <- pivot_longer(dss_weight_change, c(5:length(dss_weight)), cols_vary = "slowest", names_to = "date", values_to = "weight_change_percent")
  dss_hemo <- pivot_longer(dss_hemo, c(5:length(dss_hemo)), cols_vary = "slowest", names_to = "date", values_to = "hemo")
  dss_consistency <- pivot_longer(dss_consistency, c(5:length(dss_consistency)), cols_vary = "slowest", names_to = "date", values_to = "consistency")
  
  #putting inside this variable colnames from a individual dataframe for later
  df_colnames <- colnames(dss_weight)[0:groupInfoCols]
  
  #merging the data
  combined_df <- bind_cols(dss_weight, dss_weight_change, dss_hemo, dss_consistency)
  
  #Identify duplicated columns
  duplicated_cols <- duplicated(names(combined_df)) | duplicated(t(combined_df))
  
  #Remove duplicated columns
  combined_df <- combined_df[, !duplicated_cols, drop = FALSE]
  
  #chqnging colnames because they are messy after the merge
  colnames(combined_df)[0:4] <- df_colnames
  colnames(combined_df)[groupInfoCols+1] <- "date"
  
  return(combined_df)
  
}


