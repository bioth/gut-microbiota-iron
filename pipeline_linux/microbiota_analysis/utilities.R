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

# Transform otu_table of a phyloseq object based on the chosen transformation
transformCounts <- function(ps, transformation = "rel_ab", log_base = 10) {
  if (transformation == "rel_ab") {
    
    # Check orientation of OTU table and transpose if necessary
    if (isFALSE(taxa_are_rows(otu_table(ps)))) {
      otu_table(ps) <- t(otu_table(ps))
      message("OTU table was transposed to have ASVs as rows.")
    } else {
      message("ASVs already as rows, no need to transpose OTU table.")
    }
    
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
  } else {
    stop("Transformation method not recognized. Use 'rel_ab' for relative abundance or 'log' for log transformation.")
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