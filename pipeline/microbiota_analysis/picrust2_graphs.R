# Provide kegg pathway ids - returns their names and their hierarchies
getBriteFromKeggPathID <- function(ids) {
  library(KEGGREST)
  library(stringr)
  
  final_df <- data.frame(
    id = character(),
    name = character(),
    cat1 = character(),
    cat2 = character(),
    stringsAsFactors = FALSE
  )
  
  for (id in ids) {
    message("Processing: ", id)
    
    # Try to fetch KEGG entry
    data <- tryCatch(keggGet(id), error = function(e) {
      warning("Failed to retrieve ", id, ": ", e$message)
      return(NULL)
    })
    
    # Skip if NULL or unexpected structure
    if (is.null(data) || length(data) == 0 || !is.list(data[[1]])) {
      warning("Skipping invalid or empty result for ID: ", id)
      next
    }
    
    entry <- data[[1]]
    
    # Safely extract name
    name <- if (!is.null(entry$NAME)) entry$NAME[1] else NA
    
    # Safely extract category levels from CLASS
    if (!is.null(entry$CLASS)) {
      class_split <- str_split(entry$CLASS, pattern = "; ")[[1]]
      cat1 <- class_split[1]
      cat2 <- ifelse(length(class_split) >= 2, class_split[2], NA)
    } else {
      cat1 <- NA
      cat2 <- NA
    }
    
    # Append to result
    final_df <- rbind(
      final_df,
      data.frame(
        id = id,
        name = name,
        cat1 = cat1,
        cat2 = cat2,
        stringsAsFactors = FALSE
      )
    )
  }
  
  return(final_df)
}

# Input is a kegg abundance table and a brite mapping table generated above to produce as output a heatmap with some grouping at some brite hierarchy level
KeggPathwayHmap <- function(kegg_ab, brite_mapping, metadata, group, custom_colors_group, custom_colors_cat, hierarchy = "none"){
  
  if(hierarchy == "1"){
    cat <- "cat1"
  }
  else if(hierarchy == "2"){
    cat <- "cat2"
  }
  
  # Subset kegg_ab df with for ids in mapping df
  kegg_ab <- kegg_ab[rownames(kegg_ab) %in% c(brite_mapping$id),]
  
  # Transform df into long format
  df_long <- kegg_ab %>%
    rownames_to_column(var = "id") %>%      # Move rownames to a column
    pivot_longer(
      cols = -id,                           # Keep KO column fixed
      names_to = "sample_id",                 
      values_to = "abundance"
    )
  

  
  # Bind metadata and kegg brite mapping
  df_long <- merge(df_long, metadata, by = "sample_id")
  df_long <- merge(df_long, brite_mapping, by = "id")
  
  # Calculate z_score
  df_long <- df_long %>%
    group_by(id) %>%
    mutate(zscore = (abundance - mean(abundance)) / sd(abundance)) %>%
    ungroup()
  
  # Reorder df for cat2 and cat1
  df_long <- df_long[order(df_long$cat2),]
  df_long <- df_long[order(df_long$cat1),]
  df_long$name <- factor(df_long$name, levels = rev(unique(df_long$name)))
  
  # Hmap of pathways
  p <- ggplot(df_long, aes(x = sample_id, y = name, fill = zscore)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, space = "Lab",
                         name = "Z score") + # limit = c(-1, 1)
    scale_y_discrete(position = "right")+
    labs(x = "", y = "")
  
  if(hierarchy != "none"){
    p <- p+
      facet_grid2(
        rows = vars(.data[[cat]]), cols = vars(.data[[group]]),
        scales = "free", space = "free",
        switch = "y",                       
        strip = strip_themed(background_x = elem_list_rect(fill = c(custom_colors_group[1],custom_colors_group[2])),
                             background_y = elem_list_rect(fill = custom_col_cat),
                             text_y = elem_list_text(angle = c(0, 0)),
                             by_layer_y   = FALSE))
  }else{
    p <- p+
      facet_wrap(~.data[[group]], scales = "free_x")
  }
  
  p <- p+
    theme_minimal()+
    theme(axis.text.y = element_text(face = "bold.italic"),
          axis.text.x = element_blank(),
          axis.title.x = element_text(size = 11, face = "bold"),
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.title.y.right = element_markdown(size = 14, face = "bold", hjust = 0.5, vjust = 0),
          strip.text.x = element_text(face = "bold", size = 11, color = ifelse(hierarchy == "none", "black", "white")),
          strip.text.y = element_text(face = "bold", size = 9),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())
  
  return(p)
  
}