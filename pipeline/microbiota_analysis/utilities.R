# Utilities functions used accross different higher level functions

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

#Function to remove "/" and "-" characters from a string + lowercase
clean_string <- function(input_string) {
  result <- tolower(gsub("[/-]", "", input_string))
  return(result)
}

# Function to make strings ready for publication (capitals letter and no underscores)
pretty_string <- function(input_string) {
  result <- gsub("[/_]", " ", input_string) # replace underscores by spaces
  result <- paste0(toupper(substring(result, 1, 1)), substring(result, 2, nchar(result)))
  return(result)
}

# Function to apply different type of counts transform for ASVs (rel_ab, log and CLR)
transformCounts <- function(ps, transformation = "rel_ab", log_base = 10) {
  
  # Check orientation of OTU table and transpose if necessary
  if (isFALSE(taxa_are_rows(otu_table(ps)))) {
    otu_table(ps) <- t(otu_table(ps))
    message("OTU table was transposed to have ASVs as rows.")
  } else {
    message("ASVs already as rows, no need to transpose OTU table.")
  }
  
  if (transformation == "rel_ab") {
    # Apply prop.table for relative abundance calculation and multiply by 100 for percentage
    otu_matrix <- apply(otu_table(ps), 2, prop.table) * 100
    
    # Ensure that the result is an otu_table object
    otu_table(ps) <- otu_table(otu_matrix, taxa_are_rows = TRUE)
    
  } else if (transformation == "log") {
    
    # Check if there are any zero counts to avoid log(0) issues
    if (any(otu_table(ps) == 0)) {
      message("Warning: Zero values detected in the OTU table. Adding a small constant to avoid log(0).")
      otu_table(ps) <- otu_table(ps) + 1e-6
    }
    
    # Apply log transformation (base 10 by default)
    otu_matrix <- log(otu_table(ps), base = log_base)
    
    # Ensure that the result is an otu_table object
    otu_table(ps) <- otu_table(otu_matrix, taxa_are_rows = TRUE)
  } else if (transformation == "CLR") {
    
    # Convert OTU table to matrix
    otu_matrix <- as.matrix(otu_table(ps))
    
    # Check if there are any zero counts and replace them to avoid log(0)
    if (any(otu_matrix == 0)) {
      message("Zero values detected, replacing with small constant to avoid log(0).")
      min_non_zero <- min(as.data.frame(otu_table(ps))[as.data.frame(otu_table(ps))>0])
      pseudocount <- min_non_zero / 2
      otu_matrix[otu_matrix == 0] <- pseudocount
    }
    # Compute geometric mean per sample (column)
    geometric_mean <- apply(otu_matrix, 2, function(x) exp(mean(log(x))))
    
    # Apply CLR transformation
    clr_matrix <- sapply(seq_len(ncol(otu_matrix)), function(j) {
      log(otu_matrix[, j] / geometric_mean[j])
    })
    
    row.names(clr_matrix) <- row.names(otu_matrix)
    colnames(clr_matrix) <- colnames(otu_matrix)
    
    # Ensure that the result is an otu_table object
    otu_table(ps) <- otu_table(clr_matrix, taxa_are_rows = TRUE)
    
  } else {
    stop("Transformation method not recognized. Use 'rel_ab' for relative abundance, 'log' for log transformation, or 'CLR' for centered log ratio.")
  }
  
  return(ps) # Return the transformed phyloseq object
}

# Function to produce picrust2 required inputs
producePicrust2Inputs <- function(ps, output_dir){
  
  # Creates picrust2 output folder if does not exist yet
  existingDirCheck(paste0(output_dir, "/picrust2/"))
  existingDirCheck(paste0(output_dir, "/picrust2/input"))
  
  # Define the output file path
  output_file <- file.path(paste0(output_dir, "/picrust2/input"), "seqs.fna")
  
  # Open a connection to write the file
  file_conn <- file(output_file, open = "w")
  
  for (i in seq_along(refseq(ps))){
    
    seq = as.character(refseq(ps)[i]) # Sequence
    id = as.character(names(refseq(ps)[i])) # ID (ASV number associated with sequence)
    
    # Write the FASTA entry to the file
    writeLines(paste0(">", id), file_conn) # FASTA header
    writeLines(seq, file_conn)            # Sequence
    
  }
  
  # Close the file connection
  close(file_conn)
  message("Produced fasta file successfully.")
  
  # Produce the biom file
  # Get otu_table
  otu_table <- as.data.frame(t(otu_table(ps)))
  
  # Convert table to desired format
  otu_table <- tibble::rownames_to_column(otu_table, var = "#OTU ID")
  
  # Add the header
  header <- "# Constructed from biom file"
  
  # Write the OTU table to a text file
  output_file <- file.path(paste0(output_dir, "/picrust2/input"), "otu_table.txt")
  write(header, file = output_file)
  write.table(otu_table, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)
  message("Produced otu_table file successfully.")
}

# Function to write and save stackbarExtended sig_table
writeStackbarExtendedSigTable <- function(sig_table_list,filepath){
  
  # Initialize empty dataframe to append tables to
  table_to_write <- data.frame()
  
  # Iterate over the list of tables
  for (i in seq_along(sig_table_list)) {
    
    # Extract the table
    table <- sig_table_list[[i]]
    
    # Add a column with the name of the current table
    table$comparaison <- names(sig_table_list)[i]
    
    # Append to the master table
    table_to_write <- rbind(table_to_write, table)
  }
  
  write_xlsx(x = table_to_write, path = filepath)
  
}

# Function to check if a variable is a factor
checkIfFactor <- function(var){
  if(isFALSE(is.factor(var))){
    stop(paste(deparse(substitute(var)), "is not a factor. Please put it as a factor for function to work." ))
  }
}

# Input a phyloseq object with species level taxa information, group species together as a singular ASV (careful, this function can produce deep changes regarding the analysis performed)
groupSameSpecies <- function(ps){
  ps <- tax_glom(ps, taxrank = "Genus_species")
}

# Ensure consistent themes across all graphs generated with ggplot2
my_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black", face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.line = element_line(color = "black", size = 1),
      legend.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

# Verifies assumptions of normality and variances for further statistical tests
verifyStatsAssumptions <- function(df, group, measure){
  
  # Levene's Test for homogeneity of variance
  print(leveneTest(as.formula(paste(measure, "~", group)), data = df))
  
  # Shapiro test per group for normality assumptions
  print(by(df[[measure]], df[[group]], shapiro.test))
  
}

# Graph for F/B ratio
fbRatioGraphTimeSeries <- function(df, group, measure, time, custom_colors, custom_theme = NULL){
  
  p <- ggplot(df, aes(x = .data[[time]], y = .data[[measure]], fill = .data[[group]])) +
    
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.6,
                 color = "black") +
    scale_fill_manual(values = custom_colors)+
    labs(title = "F/B ratio overtime\n according to diet exposure", y = "F/B ratio", fill = "Group", x = "Time") # , pattern = NA
  print(p+custom_theme)
  
}