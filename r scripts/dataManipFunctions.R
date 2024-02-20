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
        df[,i] <- as.numeric(gsub("\\,", ".", df[,i]))
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
dssFollowupManipulation <- function(df, groupInfoCols, dateStart){
  
  #getting rid of empty rows if they are (dead mice or issues)
  df <- emptyRow(df)
  
  #creating 3 different dataframes for each metric (weight, hemoccult, stool consistency)
  dss_weight <- colStringSelect(dss_followup,"weight",cols_ignored = groupInfoCols) #everytime we skip the first 4 cols (diet, treatment etc)
  dss_hemo <- colStringSelect(dss_followup,"hemoccult",cols_ignored = groupInfoCols)
  dss_consistency <- colStringSelect(dss_followup,"consistency",cols_ignored = groupInfoCols)
  
  
  
  
  
  
}


