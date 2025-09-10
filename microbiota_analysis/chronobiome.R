# Function that creates a timeline for gut microbiota relative abundance composition
library(ggstream)

plot_timeline_2_groups <- function(
    ps_object,
    exp_group, # must be as factor
    time_group, # must be as factor
    sample_name,
    main_level = 'Phylum',
    sub_level = 'Family',
    threshold = 1,
    smoothing = FALSE,
    average_relab_per_group = TRUE,
    sample_name_unique = NULL,
    n_phy,
    hues = c("Oranges", "Greens", "Blues", "Purples"),
    color_bias = 2,
    custom_theme = NULL,
    additionnalAes = NULL,
    stats = FALSE
){
  
  #### 1. Data Validation & Preparation ####
  if (!("phyloseq" %in% class(ps_object)))
    stop("ps_object must be a phyloseq-class object.")
  if (!main_level %in% colnames(tax_table(ps_object)))
    stop("main_level column not found in the tax_table.")
  if (!sub_level %in% colnames(tax_table(ps_object)))
    stop("sub_level column not found in the tax_table.")
  if (!exp_group %in% names(sample_data(ps_object)))
    stop("exp_group column not found in the sample_data.")
  if(!time_group %in% names(sample_data(ps_object)))
    stop("time_variable column not found in sample_data (needed for timepoints design).")
  
  # Ensure taxa are rows
  if (!taxa_are_rows(ps_object)) {
    ps_object <- t(ps_object)
  }
  
  # Creates levels possible with time variable and group variable
  sample_data(ps_object)$combined_group <- as.factor(paste0(sample_data(ps_object)[[exp_group]],
                                                            ":", sample_data(ps_object)[[time_group]]))
  
  # Transform counts to relative abundance (%)
  ps_prop <- transform_sample_counts(ps_object, function(OTU) ((OTU / sum(OTU)) * 100))
  
  # Extract OTU, taxonomy, and metadata tables
  otu <- as.data.frame(otu_table(ps_prop))
  tax <- as.data.frame(tax_table(ps_prop))
  meta <- data.frame(sample_data(ps_prop))
  
  # Agglomerate at the desired sub-level
  ps_f <- tax_glom(ps_prop, taxrank = sub_level, NArm = FALSE)
  otu_f <- as.data.frame(otu_table(ps_f))
  tax_f <- as.data.frame(tax_table(ps_f))
  
  # Bind taxonomy and OTU tables
  otu_tax_f <- cbind(tax_f, otu_f)
  
  # Identify columns corresponding to sample abundances
  position <- ncol(tax) + 1:(ncol(otu_tax_f) - ncol(tax))
  message(paste0('\n', length(position), ' samples are analyzed \n'))
  
  # Compute mean abundance and flag features above threshold
  otu_tax_f <- otu_tax_f %>%
    rowwise() %>%
    mutate(Mean = mean(c_across(all_of(position)))) %>%
    mutate(high_abundance = case_when(Mean > threshold ~ TRUE,
                                      TRUE ~ FALSE))
  
  # Identify the top n_phy groups (by main_level)
  topx <- otu_tax_f %>%
    dplyr::group_by(!!as.name(main_level)) %>%
    dplyr::summarise(sum_top = sum(c_across(all_of(position)))) %>%
    dplyr::arrange(desc(sum_top)) %>%
    dplyr::slice_head(n = n_phy)
  
  # Mark taxa that belong to the top groups
  otu_tax_f <- mutate(otu_tax_f,
                      selected_top = !!as.name(main_level) %in% pull(topx[, main_level]))
  
  # Replace NAs with "Unknown"
  otu_tax_f[is.na(otu_tax_f)] <- "Unknown"
  
  # Create a column for legend labels (plot_taxa)
  otu_tax_f <- otu_tax_f %>%
    mutate(
      plot_taxa = case_when(
        high_abundance == TRUE & selected_top == TRUE & !!as.name(sub_level) == "Unknown" ~ paste0("Unknown ", !!as.name(main_level)),
        high_abundance == TRUE & selected_top == TRUE & !!as.name(sub_level) != "Unknown" ~ !!as.name(sub_level),
        high_abundance == FALSE & selected_top == TRUE ~ paste0("Others ", !!as.name(main_level)),
        selected_top == FALSE ~ paste0("Others ")
      )
    )
  
  # Group unknowns and "others" relative abundance values
  # 1. Identify numeric/sample columns (those you want to sum)
  sample_cols <- c(colnames(otu_f), "Mean")  # Add all sample columns here
  
  # 2. Identify metadata columns (those you want to preserve, without summing)
  metadata_cols <- setdiff(colnames(otu_tax_f), c(sample_cols, "plot_taxa"))
  
  # 3. Group by plot_taxa and keep one representative row for metadata columns
  otu_tax_f <- otu_tax_f %>%
    group_by(plot_taxa) %>%
    summarise(
      across(all_of(sample_cols), sum, na.rm = TRUE),
      across(all_of(metadata_cols), ~ first(.x)),  # keep first value for each metadata column
      .groups = "drop"
    )
  
  if (nrow(topx) != length(hues)) {
    message('the number of colors chosen (', length(hues),
            ') is different from the defined number of features to plot (', nrow(topx), ')')
  }
  
  # Add rowmeans for groups of interest => provides mean relative abundance per group within each timepoint
  for(group in levels(sample_data(ps_object)$combined_group)){
    otu_tax_f[[group]]<- rowMeans(otu_tax_f[,sample_data(ps_object)[[sample_name]][sample_data(ps_object)$combined_group == group]])
  }

  #### 2. Color Assignment ####
  df <- as_tibble(matrix(nrow = 0, ncol = length(colnames(otu_tax_f))),
                  .name_repair = ~ colnames(otu_tax_f))
  main_level_col <- c()
  
  for (i in 1:nrow(topx)) {
    top_loop <- pull(topx[, main_level])[i]
    temp <- otu_tax_f %>%
      dplyr::filter(!!as.name(main_level) == top_loop) %>%
      dplyr::arrange(desc(Mean))
    message('Among the ', main_level, ' ', top_loop, ' ',
            length(unique(temp$plot_taxa)),
            ' features at ', sub_level, ' level will be plotted')
    getPalette <- colorRampPalette(brewer.pal(length(unique(temp$plot_taxa)), hues[i]), bias = color_bias)
    col <- getPalette(length(unique(temp$plot_taxa)) + 1)
    if (length(col) >= 3) {
      main_level_col[i] <- col[length(col) - 1]
    } else {
      main_level_col[i] <- col[length(col)]
    }
    for (u in 1:length(unique(temp$plot_taxa))) {
      temp[which(temp$plot_taxa == unique(temp$plot_taxa)[u]), 'MyColors'] <- col[u + 1]
    }
    df <- rbind(df, temp)
  }
  
  # Assign black to features not in the top groups
  unselescted <- otu_tax_f %>% filter(!(!!as.name(main_level) %in% pull(topx[, main_level])))
  unselescted$MyColors <- '#000000'
  df <- rbind(df, unselescted)
  
  # Order the data and set factor levels for plotting
  df <- df %>% ungroup() %>% dplyr::arrange(desc(selected_top), !!as.name(main_level), desc(Mean))
  df$plot_taxa <- factor(df$plot_taxa, levels = unique(df$plot_taxa))
  
  # Melt the data frame (long format)
  if(average_relab_per_group){
    
    group_cols <- intersect(colnames(df), unique(meta$combined_group)) # group-level average columns only
    
    df_long <- melt(df,
                    id = c("plot_taxa", "MyColors", main_level),
                    measure.vars = group_cols,
                    variable.name = "combined_group")
    
    df_long[, main_level] <- ifelse(df_long[, main_level] %in% pull(topx[, main_level]),
                                    df_long[, main_level],
                                    "Others ")
    
  }else{
    
    df_long <- melt(df,
                    id = c(sub_level, "MyColors", main_level),
                    measure.vars = meta[, sample_name],
                    variable.name = sample_name)

    df_long <- left_join(df_long, meta, by = sample_name)

    df_long[, main_level] <- ifelse(df_long[, main_level] %in% pull(topx[, main_level]),
                                    df_long[, main_level],
                                    "Others ")
  }

  #### 3. Prepare Colors for Plotting ####
  MyColors <- df_long$MyColors
  names(MyColors) <- df_long$plot_taxa
  MyColors2 <- unique(df_long$MyColors)
  names(MyColors2) <- unique(df_long$plot_taxa)
  main_level_col[length(main_level_col) + 1] <- '#000000'
  df_long[, main_level] <- factor(df_long[, main_level], levels = unique(df_long[, main_level]))
  vec1 <- unique(df_long[, main_level])
  vec2 <- c(pull(topx[, main_level]), paste0("Others"))
  core_text_vec1 <- gsub("(<[^>]*>|\\*|\\s+$)", "", vec1)
  core_text_vec1 <- trimws(core_text_vec1)
  order_index <- match(core_text_vec1, vec2)
  main_level_col <- main_level_col[order_index]
  names(main_level_col) <- as.character(vec1)
  
  #### 4. Plotting ####
  if (average_relab_per_group) {
    df_long <- merge(df_long, meta, by = "combined_group") # Merge metadata to retrieve the timepoint and group column
    
    # Ensure your time variable is numeric (so geom_area draws it as continuous)
    df_long[[time_group]] <- as.numeric(as.character(df_long[[time_group]]))
    
    # Make sure plot_taxa is a factor in the stacking order you want
    df_long$plot_taxa <- factor(df_long$plot_taxa,
                                levels = unique(df_long$plot_taxa))
    
    if(smoothing){
      
      # Average rel_ab per group and timepoint
      df_summary <- df_long %>%
        group_by(.data[[time_group]], plot_taxa, .data[[exp_group]]) %>%
        summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
      
      # Add back other columns that are constant per group (like sub_level, main_level, etc.)
      metadata_cols <- df_long %>%
        select(plot_taxa,.data[[main_level]],MyColors) %>%
        distinct()
      
      df_summary <- df_summary %>%
        left_join(metadata_cols, by = "plot_taxa")
      colnames(df_summary)[1] <- "time_group"
      
      # Keep only what's needed for interpolation
      interp_core <- df_summary %>%
        select(time_group, plot_taxa, !!sym(exp_group), mean_value)
      
      # Do the interpolation on core data
      full_time <- seq(min(interp_core$time_group),
                       max(interp_core$time_group),
                       length.out = 10)
      
      df_interp <- interp_core %>%
        group_by(plot_taxa, !!sym(exp_group)) %>%
        complete(time_group = full_time) %>%
        arrange(plot_taxa, !!sym(exp_group), time_group) %>%
        mutate(mean_value = approx(
          x = time_group[!is.na(mean_value)],
          y = mean_value[!is.na(mean_value)],
          xout = time_group,
          rule = 2
        )$y) %>%
        ungroup()
      
      # Join metadata back in
      df_interp <- df_interp %>%
        left_join(distinct(df_summary, plot_taxa, !!sym(main_level), MyColors), by = "plot_taxa")
        
        # Build the ggplot
        p <- ggplot(df_interp,
                    aes(x = time_group,
                        y = mean_value,
                        fill = plot_taxa,
                        group = plot_taxa)) +  # Alpha goes here
          facet_wrap(as.formula(paste("~", exp_group)), ncol = 1)+
          geom_stream(type = "proportional", bw = 0.9)+
          scale_fill_manual(name = sub_level, values = MyColors2)+
          labs(x = time_group,
               y = "Relative abundance (%)",
               fill = sub_level)+
          guides(fill = guide_legend(reverse = FALSE, title = sub_level))+
          theme_minimal()
        
        if(isFALSE(is.null(additionnalAes))){
          p <- Reduce("+", c(list(p), additionnalAes))}
          
          final_plot <- p+custom_theme
          # ggsave(plot = plot, filename = paste0(path, "/", main_level, "/", main_taxon, "_", sub_level, "_abundance.png"), width = dim[1], height = dim[2], dpi = 300, bg = "white")
      
    }else{
      
      # Build the ggplot
      final_plot <- ggplot(df_long,
                  aes(x = .data[[time_group]],
                      y = value,
                      fill = plot_taxa,
                      group = plot_taxa,
                      alpha = .data[[main_level]])) +  # Alpha goes here
        facet_wrap(as.formula(paste("~", exp_group)), ncol = 1) +
        geom_area(show.legend = TRUE) +
        scale_alpha_manual(
          values = rep(1, length(unique(df_long[[main_level]]))),  # Set to 1 to show legend, or actual transparency values if desired
          guide = guide_legend(order = 1, override.aes = list(fill = main_level_col)),
          name = main_level
        ) +
        scale_fill_manual(name = sub_level, values = MyColors2) +
        labs(x = time_group,
             y = "Relative abundance (%)",
             fill = sub_level,
             alpha = main_level) +
        guides(fill = guide_legend(reverse = FALSE, title = sub_level, order = 2)) +
        theme_minimal()+
        custom_theme
      
      if(isFALSE(is.null(additionnalAes))){
        p <- Reduce("+", c(list(p), additionnalAes))
      }
      
    }
    
  }
  else{
    
    # Ensure your time variable is numeric (so geom_area draws it as continuous)
    df_long[[time_group]] <- as.numeric(as.character(df_long[[time_group]]))
    

    
    # Make sure plot_taxa is a factor in the stacking order you want
    df_long[[sub_level]] <- factor(df_long[[sub_level]],
                                levels = unique(df_long[[sub_level]]))
    
    colnames(df_long)[13] <- "group"

    plots <- df_long %>%
      split(.$group) %>%
      lapply(function(sub_df) {
        final_plot <- ggplot(sub_df,
                            aes(x = .data[[time_group]],
                                y = value,
                                fill = .data[[sub_level]],
                                group = .data[[sub_level]],
                                alpha = .data[[main_level]])) +
          facet_wrap(~ .data[[sample_name_unique]], ncol = 1, scales = "free_y") +
          geom_area(position = "stack") +
          ggtitle(unique(sub_df$group)) +
          scale_alpha_manual(
            values = rep(1, length(unique(df_long[[main_level]]))),
            guide = guide_legend(order = 1, override.aes = list(fill = main_level_col)),
            name = main_level
          ) +
          scale_fill_manual(name = sub_level, values = MyColors2) +
          labs(x = time_group,
               y = "Relative abundance (%)",
               fill = sub_level,
               alpha = main_level) +
          guides(fill = guide_legend(reverse = FALSE, title = sub_level, order = 2)) +
          theme_minimal()
      })
    
    # # Remove legends from all but the first plot
    # plots_shared_legend <- plots
    # plots_shared_legend[-1] <- lapply(plots_shared_legend[-1], function(p) p + theme(legend.position = "none"))
    # wrap_plots(plots_shared_legend, ncol = 4, guides = "collect") &
    #   theme(legend.position = "bottom")

    
    
    # # Prepare variables for plotting
    # sample_var <- ensym(sample_name)
    # group_var <- ensym(exp_group)
    # 
    # df_long$facet_var <- interaction(df_long[[as.character(group_var)]], df_long[[as.character(sample_var)]])
    
    # Build the ggplot
    # p <- ggplot(df_long,
    #             aes(x = .data[[time_group]],
    #                 y = value,
    #                 fill = .data[[sub_level]],
    #                 group = .data[[sub_level]],
    #                 alpha = .data[[main_level]])) +  # Alpha goes here
    #   # facet_wrap(as.formula(paste("~", exp_group)), ncol = 1) +
    #   facet_grid2(
    #     rows = vars(.data[[sample_name]]), cols = vars(.data[[exp_group]]),
    #     scales = "free", space = "free",
    #     switch = "y")+                       
    #     # strip = strip_themed(background_x = elem_list_rect(fill = c(custom_colors_group[1],custom_colors_group[2])),
    #     #                      background_y = elem_list_rect(fill = custom_col_cat),
    #     #                      text_y = elem_list_text(angle = c(0, 0)),
    #     #                      by_layer_y   = FALSE))
    #   geom_area(show.legend = TRUE)
      # facet_wrap(~ facet_var, ncol = 4)
      
    
    #   scale_alpha_manual(
    #     values = rep(1, length(unique(df_long[[main_level]]))),  # Set to 1 to show legend, or actual transparency values if desired
    #     guide = guide_legend(order = 1, override.aes = list(fill = main_level_col)),
    #     name = main_level
    #   ) +
    #   scale_fill_manual(name = sub_level, values = MyColors2) +
    #   labs(x = time_group,
    #        y = "Relative abundance (%)",
    #        fill = sub_level,
    #        alpha = main_level) +
    #   guides(fill = guide_legend(reverse = FALSE, title = sub_level, order = 2)) +
    #   theme_minimal()
    # 
    # if(isFALSE(is.null(additionnalAes))){
    #   p <- Reduce("+", c(list(p), additionnalAes))
    # }
    
    # return(p+custom_theme)


    
  }
  #### 5. Stats ####
  if (stats) {
   subTaxaList <- unique(df_long$plot_taxa)
   mainTaxaList <- unique(df_long[[main_level]])
   
   sub_stats <- timelineStatTable(ps_object, exp_group, time_group, sub_level,
                                 top_n = 30, p_col = "padj", p_sig = 0.05,
                                 title = NULL, taxaList = subTaxaList)
   
   main_stats <- timelineStatTable(ps_object, exp_group, time_group, main_level,
                                  top_n = 30, p_col = "padj", p_sig = 0.05,
                                  title = NULL, taxaList = mainTaxaList)
    
  }
  
  # Return final objects
  chronobiome <- NULL
  chronobiome[["main_fig"]] <- final_plot
  chronobiome[["main_names"]] <- unique(df_long[[main_level]])
  chronobiome[["sub_names"]] <- unique(df_long[["plot_taxa"]])
  if(stats){
    chronobiome[["sub_stat_df"]] <- sub_stats$res_all
    chronobiome[["sub_stat_plot"]] <- sub_stats$p
    chronobiome[["main_stat_df"]] <- main_stats$res_all
    chronobiome[["main_stat_plot"]] <- main_stats$p
  }
  return(chronobiome)
}


# Function that creates a timeline for gut microbiota relative abundance composition - specific to a main_level taxa and its assoicated sub_level taxons, multiple graphs
plot_timeline_taxa <- function(
    ps_object,
    exp_group, # must be as factor
    time_group, # must be as factor
    sample_name,
    main_level = 'Phylum',
    sub_level = 'Family',
    threshold = 1,
    average_relab_per_group = TRUE,
    smoothing = FALSE,
    n_phy,
    hues = c("Oranges", "Greens", "Blues", "Purples"),
    color_bias = 2,
    custom_theme = NULL,
    additionnalAes = NULL,
    dim = c(6,6),
    path = NULL
){
  
  #### 1. Data Validation & Preparation ####
  if (!("phyloseq" %in% class(ps_object)))
    stop("ps_object must be a phyloseq-class object.")
  if (!main_level %in% colnames(tax_table(ps_object)))
    stop("main_level column not found in the tax_table.")
  if (!sub_level %in% colnames(tax_table(ps_object)))
    stop("sub_level column not found in the tax_table.")
  if (!exp_group %in% names(sample_data(ps_object)))
    stop("exp_group column not found in the sample_data.")
  if(!time_group %in% names(sample_data(ps_object)))
    stop("time_variable column not found in sample_data (needed for timepoints design).")
  if(is.null(path))
    stop("No defined path to save the graphs")
  
  # Create folder where to save graphs
  existingDirCheck(paste0(path, "/", main_level))
  
  # Ensure taxa are rows
  if (!taxa_are_rows(ps_object)) {
    ps_object <- t(ps_object)
  }
  
  # Creates levels possible with time variable and group variable
  sample_data(ps_object)$combined_group <- as.factor(paste0(sample_data(ps_object)[[exp_group]],
                                                            ":", sample_data(ps_object)[[time_group]]))
  
  # Transform counts to relative abundance (%)
  ps_prop <- transform_sample_counts(ps_object, function(OTU) ((OTU / sum(OTU)) * 100))
  
  # Extract OTU, taxonomy, and metadata tables
  otu <- as.data.frame(otu_table(ps_prop))
  tax <- as.data.frame(tax_table(ps_prop))
  meta <- data.frame(sample_data(ps_prop))
  
  # Agglomerate at the desired sub-level
  ps_f <- tax_glom(ps_prop, taxrank = sub_level, NArm = FALSE)
  otu_f <- as.data.frame(otu_table(ps_f))
  tax_f <- as.data.frame(tax_table(ps_f))
  
  # Bind taxonomy and OTU tables
  otu_tax_f <- cbind(tax_f, otu_f)
  
  # Identify columns corresponding to sample abundances
  position <- ncol(tax) + 1:(ncol(otu_tax_f) - ncol(tax))
  message(paste0('\n', length(position), ' samples are analyzed \n'))
  
  # Compute mean abundance and flag features above threshold
  otu_tax_f <- otu_tax_f %>%
    rowwise() %>%
    mutate(Mean = mean(c_across(all_of(position)))) %>%
    mutate(high_abundance = case_when(Mean > threshold ~ TRUE,
                                      TRUE ~ FALSE))
  
  # Identify the top n_phy groups (by main_level)
  topx <- otu_tax_f %>%
    dplyr::group_by(!!as.name(main_level)) %>%
    dplyr::summarise(sum_top = sum(c_across(all_of(position)))) %>%
    dplyr::arrange(desc(sum_top)) %>%
    dplyr::slice_head(n = n_phy)
  
  # Mark taxa that belong to the top groups
  otu_tax_f <- mutate(otu_tax_f,
                      selected_top = !!as.name(main_level) %in% pull(topx[, main_level]))
  
  # Replace NAs with "Unknown"
  otu_tax_f[is.na(otu_tax_f)] <- "Unknown"
  
  # Create a column for legend labels (plot_taxa)
  otu_tax_f <- otu_tax_f %>%
    mutate(
      plot_taxa = case_when(
        high_abundance == TRUE & selected_top == TRUE & !!as.name(sub_level) == "Unknown" ~ paste0("Unknown ", !!as.name(main_level)),
        high_abundance == TRUE & selected_top == TRUE & !!as.name(sub_level) != "Unknown" ~ !!as.name(sub_level),
        high_abundance == FALSE & selected_top == TRUE ~ paste0("Others ", !!as.name(main_level)),
        selected_top == FALSE ~ paste0("Others ")
      )
    )
  
  # Group unknowns and "others" relative abundance values
  # 1. Identify numeric/sample columns (those you want to sum)
  sample_cols <- c(colnames(otu_f), "Mean")  # Add all sample columns here
  
  # 2. Identify metadata columns (those you want to preserve, without summing)
  metadata_cols <- setdiff(colnames(otu_tax_f), c(sample_cols, "plot_taxa"))
  
  # 3. Group by plot_taxa and keep one representative row for metadata columns
  otu_tax_f <- otu_tax_f %>%
    group_by(plot_taxa) %>%
    summarise(
      across(all_of(sample_cols), sum, na.rm = TRUE),
      across(all_of(metadata_cols), ~ first(.x)),  # keep first value for each metadata column
      .groups = "drop"
    )
  
  if (nrow(topx) != length(hues)) {
    message('the number of colors chosen (', length(hues),
            ') is different from the defined number of features to plot (', nrow(topx), ')')
  }
  
  # Add rowmeans for groups of interest => provides mean relative abundance per group within each timepoint
  for(group in levels(sample_data(ps_object)$combined_group)){
    otu_tax_f[[group]]<- rowMeans(otu_tax_f[,sample_data(ps_object)[[sample_name]][sample_data(ps_object)$combined_group == group]])
  }
  
  #### 2. Color Assignment ####
  df <- as_tibble(matrix(nrow = 0, ncol = length(colnames(otu_tax_f))),
                  .name_repair = ~ colnames(otu_tax_f))
  main_level_col <- c()
  
  for (i in 1:nrow(topx)) {
    top_loop <- pull(topx[, main_level])[i]
    temp <- otu_tax_f %>%
      dplyr::filter(!!as.name(main_level) == top_loop) %>%
      dplyr::arrange(desc(Mean))
    message('Among the ', main_level, ' ', top_loop, ' ',
            length(unique(temp$plot_taxa)),
            ' features at ', sub_level, ' level will be plotted')
    getPalette <- colorRampPalette(brewer.pal(length(unique(temp$plot_taxa)), hues[i]), bias = color_bias)
    col <- getPalette(length(unique(temp$plot_taxa)) + 1)
    if (length(col) >= 3) {
      main_level_col[i] <- col[length(col) - 1]
    } else {
      main_level_col[i] <- col[length(col)]
    }
    for (u in 1:length(unique(temp$plot_taxa))) {
      temp[which(temp$plot_taxa == unique(temp$plot_taxa)[u]), 'MyColors'] <- col[u + 1]
    }
    df <- rbind(df, temp)
  }
  
  # Assign black to features not in the top groups
  unselescted <- otu_tax_f %>% filter(!(!!as.name(main_level) %in% pull(topx[, main_level])))
  unselescted$MyColors <- '#000000'
  df <- rbind(df, unselescted)
  
  # Order the data and set factor levels for plotting
  df <- df %>% ungroup() %>% dplyr::arrange(desc(selected_top), !!as.name(main_level), desc(Mean))
  df$plot_taxa <- factor(df$plot_taxa, levels = unique(df$plot_taxa))
  
  # Melt the data frame (long format)
  if(average_relab_per_group){
    
    group_cols <- intersect(colnames(df), unique(meta$combined_group)) # group-level average columns only
    
    df_long <- melt(df,
                    id = c("plot_taxa", "MyColors", main_level),
                    measure.vars = group_cols,
                    variable.name = "combined_group")
    
    df_long[, main_level] <- ifelse(df_long[, main_level] %in% pull(topx[, main_level]),
                                    df_long[, main_level],
                                    "Others ")
    
  }else{
    df_long <- melt(df,
                    id = c(sub_level, "MyColors", main_level),
                    measure.vars = meta[, sample_name],
                    variable.name = sample_name)
    df_long <- left_join(df_long, meta, by = sample_name)
    df_long[, main_level] <- ifelse(df_long[, main_level] %in% pull(topx[, main_level]),
                                    df_long[, main_level],
                                    "Others ")
  }
  
  
  #### 3. Prepare Colors for Plotting ####
  MyColors <- df_long$MyColors
  names(MyColors) <- df_long$plot_taxa
  MyColors2 <- unique(df_long$MyColors)
  names(MyColors2) <- unique(df_long$plot_taxa)
  main_level_col[length(main_level_col) + 1] <- '#000000'
  df_long[, main_level] <- factor(df_long[, main_level], levels = unique(df_long[, main_level]))
  vec1 <- unique(df_long[, main_level])
  vec2 <- c(pull(topx[, main_level]), paste0("Others"))
  core_text_vec1 <- gsub("(<[^>]*>|\\*|\\s+$)", "", vec1)
  core_text_vec1 <- trimws(core_text_vec1)
  order_index <- match(core_text_vec1, vec2)
  main_level_col <- main_level_col[order_index]
  names(main_level_col) <- as.character(vec1)
  
  #### 4. Plotting ####
  if (average_relab_per_group) {
    df_long <- merge(df_long, meta, by = "combined_group") # Merge metadata to retrieve the timepoint and group column
    
    # Ensure your time variable is numeric (so geom_area draws it as continuous)
    df_long[[time_group]] <- as.numeric(as.character(df_long[[time_group]]))
    
    # Make sure plot_taxa is a factor in the stacking order you want
    df_long$plot_taxa <- factor(df_long$plot_taxa,
                                levels = unique(df_long$plot_taxa))
    
    # Average rel_ab per group and timepoint
    df_summary <- df_long %>%
    group_by(.data[[time_group]], plot_taxa, .data[[exp_group]]) %>%
      summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
    
    # Add back other columns that are constant per group (like sub_level, main_level, etc.)
    metadata_cols <- df_long %>%
      select(plot_taxa,.data[[main_level]],MyColors) %>%
      distinct()
    
    df_summary <- df_summary %>%
      left_join(metadata_cols, by = "plot_taxa")
    colnames(df_summary)[1] <- "time_group"
    
    # Keep only what's needed for interpolation
    interp_core <- df_summary %>%
      select(time_group, plot_taxa, !!sym(exp_group), mean_value)
    
    # Do the interpolation on core data
    full_time <- seq(min(interp_core$time_group),
                     max(interp_core$time_group),
                     length.out = 20)
    
    df_interp <- interp_core %>%
      group_by(plot_taxa, !!sym(exp_group)) %>%
      complete(time_group = full_time) %>%
      arrange(plot_taxa, !!sym(exp_group), time_group) %>%
      mutate(mean_value = approx(
        x = time_group[!is.na(mean_value)],
        y = mean_value[!is.na(mean_value)],
        xout = time_group,
        rule = 2
      )$y) %>%
      ungroup()
    
    # Join metadata back in
    df_interp <- df_interp %>%
      left_join(distinct(df_summary, plot_taxa, !!sym(main_level), MyColors), by = "plot_taxa")
    
    # Create plot for each main taxa lvl
    for(main_taxon in levels(df_interp[[main_level]])){
      
      # Subset df for taxon of interest
      df_sub <- df_interp[df_interp[[main_level]] == main_taxon,]
      
      # Build the ggplot
      p <- ggplot(df_sub,
                  aes(x = time_group,
                      y = mean_value,
                      fill = plot_taxa,
                      group = plot_taxa)) +  # Alpha goes here
        facet_wrap(as.formula(paste("~", exp_group)), ncol = 1)
      
      if(smoothing){
        p <- p+
          geom_stream(type = "ridge", bw = 0.9)
      }else{
        p <- p+
          geom_area(show.legend = TRUE)
      }
      
      p <- p+
        scale_fill_manual(name = sub_level, values = MyColors2) +
        labs(x = time_group,
             y = "Relative abundance (%)",
             fill = sub_level,
             title = paste(main_level, main_taxon)) +
        guides(fill = guide_legend(reverse = FALSE, title = sub_level)) +
        theme_minimal()
        
      
      if(isFALSE(is.null(additionnalAes))){
        p <- Reduce("+", c(list(p), additionnalAes))
      
      plot <- p+custom_theme
      ggsave(plot = plot, filename = paste0(path, "/", main_level, "/", main_taxon, "_", sub_level, "_abundance.png"), width = dim[1], height = dim[2], dpi = 300, bg = "white")
      }
    }
  }
}

plot_timeline_selected_taxa <- function(
    ps_object,
    exp_group, # must be as factor
    time_group, # must be as factor
    sample_name,
    main_level = 'Phylum',
    sub_level = 'Family',
    selected_taxa = 'ASV1',
    threshold = 1,
    smoothing = FALSE,
    average_relab_per_group = TRUE,
    sample_name_unique = NULL,
    n_phy,
    hues = c("Oranges", "Greens", "Blues", "Purples"),
    color_bias = 2,
    custom_theme = NULL,
    additionnalAes = NULL,
    stats = FALSE
){
  
  #### 1. Data Validation & Preparation ####
  if (!("phyloseq" %in% class(ps_object)))
    stop("ps_object must be a phyloseq-class object.")
  if (!main_level %in% colnames(tax_table(ps_object)))
    stop("main_level column not found in the tax_table.")
  if (!sub_level %in% colnames(tax_table(ps_object)))
    stop("sub_level column not found in the tax_table.")
  if (!exp_group %in% names(sample_data(ps_object)))
    stop("exp_group column not found in the sample_data.")
  if(!time_group %in% names(sample_data(ps_object)))
    stop("time_variable column not found in sample_data (needed for timepoints design).")
  
  # Ensure taxa are rows
  if (!taxa_are_rows(ps_object)) {
    ps_object <- t(ps_object)
  }
  
  # Creates levels possible with time variable and group variable
  sample_data(ps_object)$combined_group <- as.factor(paste0(sample_data(ps_object)[[exp_group]],
                                                            ":", sample_data(ps_object)[[time_group]]))
  
  # Transform counts to relative abundance (%)
  ps_prop <- transform_sample_counts(ps_object, function(OTU) ((OTU / sum(OTU)) * 100))
  
  # Extract OTU, taxonomy, and metadata tables
  otu <- as.data.frame(otu_table(ps_prop))
  tax <- as.data.frame(tax_table(ps_prop))
  meta <- data.frame(sample_data(ps_prop))
  
  # Agglomerate at the desired sub-level
  ps_f <- tax_glom(ps_prop, taxrank = sub_level, NArm = FALSE)
  otu_f <- as.data.frame(otu_table(ps_f))
  tax_f <- as.data.frame(tax_table(ps_f))
  
  # Bind taxonomy and OTU tables
  otu_tax_f <- cbind(tax_f, otu_f)
  
  # Keep only taxa of interest
  otu_tax_f <- otu_tax_f[otu_tax_f[[main_level]] %in% selected_taxa,]
  
  # Identify columns corresponding to sample abundances
  position <- ncol(tax) + 1:(ncol(otu_tax_f) - ncol(tax))
  message(paste0('\n', length(position), ' samples are analyzed \n'))
  
  # Compute mean abundance and flag features above threshold
  otu_tax_f <- otu_tax_f %>%
    rowwise() %>%
    mutate(Mean = mean(c_across(all_of(position)))) %>%
    mutate(high_abundance = case_when(Mean > threshold ~ TRUE,
                                      TRUE ~ FALSE))
  
  # Identify the top n_phy groups (by main_level)
  topx <- otu_tax_f %>%
    dplyr::group_by(!!as.name(main_level)) %>%
    dplyr::summarise(sum_top = sum(c_across(all_of(position)))) %>%
    dplyr::arrange(desc(sum_top)) %>%
    dplyr::slice_head(n = n_phy)
  
  # Mark taxa that belong to the top groups
  otu_tax_f <- mutate(otu_tax_f,
                      selected_top = !!as.name(main_level) %in% pull(topx[, main_level]))
  
  # Replace NAs with "Unknown"
  otu_tax_f[is.na(otu_tax_f)] <- "Unknown"
  
  # Create a column for legend labels (plot_taxa)
  otu_tax_f <- otu_tax_f %>%
    mutate(
      plot_taxa = case_when(
        high_abundance == TRUE & selected_top == TRUE & !!as.name(sub_level) == "Unknown" ~ paste0("Unknown ", !!as.name(main_level)),
        high_abundance == TRUE & selected_top == TRUE & !!as.name(sub_level) != "Unknown" ~ !!as.name(sub_level),
        high_abundance == FALSE & selected_top == TRUE ~ paste0("Others ", !!as.name(main_level)),
        selected_top == FALSE ~ paste0("Others ")
      )
    )
  
  # Group unknowns and "others" relative abundance values
  # 1. Identify numeric/sample columns (those you want to sum)
  sample_cols <- c(colnames(otu_f), "Mean")  # Add all sample columns here
  
  # 2. Identify metadata columns (those you want to preserve, without summing)
  metadata_cols <- setdiff(colnames(otu_tax_f), c(sample_cols, "plot_taxa"))
  
  # 3. Group by plot_taxa and keep one representative row for metadata columns
  otu_tax_f <- otu_tax_f %>%
    group_by(plot_taxa) %>%
    summarise(
      across(all_of(sample_cols), sum, na.rm = TRUE),
      across(all_of(metadata_cols), ~ first(.x)),  # keep first value for each metadata column
      .groups = "drop"
    )
  
  if (nrow(topx) != length(hues)) {
    message('the number of colors chosen (', length(hues),
            ') is different from the defined number of features to plot (', nrow(topx), ')')
  }
  
  # Add rowmeans for groups of interest => provides mean relative abundance per group within each timepoint
  for(group in levels(sample_data(ps_object)$combined_group)){
    otu_tax_f[[group]]<- rowMeans(otu_tax_f[,sample_data(ps_object)[[sample_name]][sample_data(ps_object)$combined_group == group]])
  }
  
  #### 2. Color Assignment ####
  df <- as_tibble(matrix(nrow = 0, ncol = length(colnames(otu_tax_f))),
                  .name_repair = ~ colnames(otu_tax_f))
  main_level_col <- c()
  
  for (i in 1:nrow(topx)) {
    top_loop <- pull(topx[, main_level])[i]
    temp <- otu_tax_f %>%
      dplyr::filter(!!as.name(main_level) == top_loop) %>%
      dplyr::arrange(desc(Mean))
    message('Among the ', main_level, ' ', top_loop, ' ',
            length(unique(temp$plot_taxa)),
            ' features at ', sub_level, ' level will be plotted')
    getPalette <- colorRampPalette(brewer.pal(length(unique(temp$plot_taxa)), hues[i]), bias = color_bias)
    col <- getPalette(length(unique(temp$plot_taxa)) + 1)
    if (length(col) >= 3) {
      main_level_col[i] <- col[length(col) - 1]
    } else {
      main_level_col[i] <- col[length(col)]
    }
    for (u in 1:length(unique(temp$plot_taxa))) {
      temp[which(temp$plot_taxa == unique(temp$plot_taxa)[u]), 'MyColors'] <- col[u + 1]
    }
    df <- rbind(df, temp)
  }
  
  # Assign black to features not in the top groups
  unselescted <- otu_tax_f %>% filter(!(!!as.name(main_level) %in% pull(topx[, main_level])))
  unselescted$MyColors <- '#000000'
  df <- rbind(df, unselescted)
  
  # Order the data and set factor levels for plotting
  df <- df %>% ungroup() %>% dplyr::arrange(desc(selected_top), !!as.name(main_level), desc(Mean))
  df$plot_taxa <- factor(df$plot_taxa, levels = unique(df$plot_taxa))
  
  # Melt the data frame (long format)
  if(average_relab_per_group){
    
    group_cols <- intersect(colnames(df), unique(meta$combined_group)) # group-level average columns only
    
    df_long <- melt(df,
                    id = c("plot_taxa", "MyColors", main_level),
                    measure.vars = group_cols,
                    variable.name = "combined_group")
    
    df_long[, main_level] <- ifelse(df_long[, main_level] %in% pull(topx[, main_level]),
                                    df_long[, main_level],
                                    "Others ")
    
  }else{
    
    df_long <- melt(df,
                    id = c(sub_level, "MyColors", main_level),
                    measure.vars = meta[, sample_name],
                    variable.name = sample_name)
    
    df_long <- left_join(df_long, meta, by = sample_name)
    
    df_long[, main_level] <- ifelse(df_long[, main_level] %in% pull(topx[, main_level]),
                                    df_long[, main_level],
                                    "Others ")
  }
  
  #### 3. Prepare Colors for Plotting ####
  MyColors <- df_long$MyColors
  names(MyColors) <- df_long$plot_taxa
  MyColors2 <- unique(df_long$MyColors)
  names(MyColors2) <- unique(df_long$plot_taxa)
  main_level_col[length(main_level_col) + 1] <- '#000000'
  df_long[, main_level] <- factor(df_long[, main_level], levels = unique(df_long[, main_level]))
  vec1 <- unique(df_long[, main_level])
  vec2 <- c(pull(topx[, main_level]), paste0("Others"))
  core_text_vec1 <- gsub("(<[^>]*>|\\*|\\s+$)", "", vec1)
  core_text_vec1 <- trimws(core_text_vec1)
  order_index <- match(core_text_vec1, vec2)
  main_level_col <- main_level_col[order_index]
  names(main_level_col) <- as.character(vec1)
  
  #### 4. Plotting ####
  if (average_relab_per_group) {
    df_long <- merge(df_long, meta, by = "combined_group") # Merge metadata to retrieve the timepoint and group column
    
    # Ensure your time variable is numeric (so geom_area draws it as continuous)
    df_long[[time_group]] <- as.numeric(as.character(df_long[[time_group]]))
    
    # Make sure plot_taxa is a factor in the stacking order you want
    df_long$plot_taxa <- factor(df_long$plot_taxa,
                                levels = unique(df_long$plot_taxa))
    
    
    if(smoothing){
      
      # Average rel_ab per group and timepoint
      df_summary <- df_long %>%
        group_by(.data[[time_group]], plot_taxa, .data[[exp_group]]) %>%
        summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
      
      # Add back other columns that are constant per group (like sub_level, main_level, etc.)
      metadata_cols <- df_long %>%
        select(plot_taxa,.data[[main_level]],MyColors) %>%
        distinct()
      
      df_summary <- df_summary %>%
        left_join(metadata_cols, by = "plot_taxa")
      colnames(df_summary)[1] <- "time_group"
      
      # Keep only what's needed for interpolation
      interp_core <- df_summary %>%
        select(time_group, plot_taxa, !!sym(exp_group), mean_value)
      
      # Do the interpolation on core data
      full_time <- seq(min(interp_core$time_group),
                       max(interp_core$time_group),
                       length.out = 10)
      
      df_interp <- interp_core %>%
        group_by(plot_taxa, !!sym(exp_group)) %>%
        complete(time_group = full_time) %>%
        arrange(plot_taxa, !!sym(exp_group), time_group) %>%
        mutate(mean_value = approx(
          x = time_group[!is.na(mean_value)],
          y = mean_value[!is.na(mean_value)],
          xout = time_group,
          rule = 2
        )$y) %>%
        ungroup()
      
      # Join metadata back in
      df_interp <- df_interp %>%
        left_join(distinct(df_summary, plot_taxa, !!sym(main_level), MyColors), by = "plot_taxa")
      
      # Build the ggplot
      p <- ggplot(df_interp,
                  aes(x = time_group,
                      y = mean_value,
                      fill = plot_taxa,
                      group = plot_taxa)) +  # Alpha goes here
        facet_wrap(as.formula(paste("~", exp_group)), ncol = 1)+
        geom_stream(type = "proportional", bw = 0.9)+
        scale_fill_manual(name = sub_level, values = MyColors2)+
        labs(x = time_group,
             y = "Relative abundance (%)",
             fill = sub_level)+
        guides(fill = guide_legend(reverse = FALSE, title = sub_level))+
        theme_minimal()
      
      if(isFALSE(is.null(additionnalAes))){
        p <- Reduce("+", c(list(p), additionnalAes))}
      
      final_plot <- p+custom_theme
      # ggsave(plot = plot, filename = paste0(path, "/", main_level, "/", main_taxon, "_", sub_level, "_abundance.png"), width = dim[1], height = dim[2], dpi = 300, bg = "white")
      
    }else{
      
      # Build the ggplot
      final_plot <- ggplot(df_long,
                           aes(x = .data[[time_group]],
                               y = value,
                               fill = plot_taxa,
                               group = plot_taxa,
                               alpha = .data[[main_level]])) +  # Alpha goes here
        facet_wrap(as.formula(paste("~", exp_group)), ncol = 1) +
        geom_area(show.legend = TRUE) +
        scale_alpha_manual(
          values = rep(1, length(unique(df_long[[main_level]]))),  # Set to 1 to show legend, or actual transparency values if desired
          guide = guide_legend(order = 1, override.aes = list(fill = main_level_col)),
          name = main_level
        ) +
        scale_fill_manual(name = sub_level, values = MyColors2) +
        labs(x = time_group,
             y = "Relative abundance (%)",
             fill = sub_level,
             alpha = main_level) +
        guides(fill = guide_legend(reverse = FALSE, title = sub_level, order = 2)) +
        theme_minimal()+
        custom_theme
      
      if(isFALSE(is.null(additionnalAes))){
        p <- Reduce("+", c(list(p), additionnalAes))
      }
      
    }
    
  }
  else{
    
    # Ensure your time variable is numeric (so geom_area draws it as continuous)
    df_long[[time_group]] <- as.numeric(as.character(df_long[[time_group]]))
    
    
    
    # Make sure plot_taxa is a factor in the stacking order you want
    df_long[[sub_level]] <- factor(df_long[[sub_level]],
                                   levels = unique(df_long[[sub_level]]))
    
    colnames(df_long)[13] <- "group"
    
    plots <- df_long %>%
      split(.$group) %>%
      lapply(function(sub_df) {
        final_plot <- ggplot(sub_df,
                             aes(x = .data[[time_group]],
                                 y = value,
                                 fill = .data[[sub_level]],
                                 group = .data[[sub_level]],
                                 alpha = .data[[main_level]])) +
          facet_wrap(~ .data[[sample_name_unique]], ncol = 1, scales = "free_y") +
          geom_area(position = "stack") +
          ggtitle(unique(sub_df$group)) +
          scale_alpha_manual(
            values = rep(1, length(unique(df_long[[main_level]]))),
            guide = guide_legend(order = 1, override.aes = list(fill = main_level_col)),
            name = main_level
          ) +
          scale_fill_manual(name = sub_level, values = MyColors2) +
          labs(x = time_group,
               y = "Relative abundance (%)",
               fill = sub_level,
               alpha = main_level) +
          guides(fill = guide_legend(reverse = FALSE, title = sub_level, order = 2)) +
          theme_minimal()
      })
    
    # # Remove legends from all but the first plot
    # plots_shared_legend <- plots
    # plots_shared_legend[-1] <- lapply(plots_shared_legend[-1], function(p) p + theme(legend.position = "none"))
    # wrap_plots(plots_shared_legend, ncol = 4, guides = "collect") &
    #   theme(legend.position = "bottom")
    
    
    
    # # Prepare variables for plotting
    # sample_var <- ensym(sample_name)
    # group_var <- ensym(exp_group)
    # 
    # df_long$facet_var <- interaction(df_long[[as.character(group_var)]], df_long[[as.character(sample_var)]])
    
    # Build the ggplot
    # p <- ggplot(df_long,
    #             aes(x = .data[[time_group]],
    #                 y = value,
    #                 fill = .data[[sub_level]],
    #                 group = .data[[sub_level]],
    #                 alpha = .data[[main_level]])) +  # Alpha goes here
    #   # facet_wrap(as.formula(paste("~", exp_group)), ncol = 1) +
    #   facet_grid2(
    #     rows = vars(.data[[sample_name]]), cols = vars(.data[[exp_group]]),
    #     scales = "free", space = "free",
    #     switch = "y")+                       
    #     # strip = strip_themed(background_x = elem_list_rect(fill = c(custom_colors_group[1],custom_colors_group[2])),
    #     #                      background_y = elem_list_rect(fill = custom_col_cat),
    #     #                      text_y = elem_list_text(angle = c(0, 0)),
    #     #                      by_layer_y   = FALSE))
    #   geom_area(show.legend = TRUE)
    # facet_wrap(~ facet_var, ncol = 4)
    
    
    #   scale_alpha_manual(
    #     values = rep(1, length(unique(df_long[[main_level]]))),  # Set to 1 to show legend, or actual transparency values if desired
    #     guide = guide_legend(order = 1, override.aes = list(fill = main_level_col)),
    #     name = main_level
    #   ) +
    #   scale_fill_manual(name = sub_level, values = MyColors2) +
    #   labs(x = time_group,
    #        y = "Relative abundance (%)",
    #        fill = sub_level,
    #        alpha = main_level) +
    #   guides(fill = guide_legend(reverse = FALSE, title = sub_level, order = 2)) +
    #   theme_minimal()
    # 
    # if(isFALSE(is.null(additionnalAes))){
    #   p <- Reduce("+", c(list(p), additionnalAes))
    # }
    
    # return(p+custom_theme)
    
    
    
  }
  #### 5. Stats ####
  if (stats) {
    subTaxaList <- unique(df_long$plot_taxa)
    mainTaxaList <- unique(df_long[[main_level]])
    
    sub_stats <- timelineStatTable(ps_object, exp_group, time_group, sub_level,
                                   top_n = 30, p_col = "padj", p_sig = 0.05,
                                   title = NULL, taxaList = subTaxaList)
    
    main_stats <- timelineStatTable(ps_object, exp_group, time_group, main_level,
                                    top_n = 30, p_col = "padj", p_sig = 0.05,
                                    title = NULL, taxaList = mainTaxaList)
    
  }
  
  # Return final objects
  chronobiome <- NULL
  chronobiome[["main_fig"]] <- final_plot
  chronobiome[["main_names"]] <- unique(df_long[[main_level]])
  chronobiome[["sub_names"]] <- unique(df_long[["plot_taxa"]])
  if(stats){
    chronobiome[["sub_stat_df"]] <- sub_stats$res_all
    chronobiome[["sub_stat_plot"]] <- sub_stats$p
    chronobiome[["main_stat_df"]] <- main_stats$res_all
    chronobiome[["main_stat_plot"]] <- main_stats$p
  }
  return(chronobiome)
}


# Function that produces a stat table for comparison between 2 groups with DESeq2 for different timepoints
# Timepoints and group variables must be ordered as factors
timelineStatTable <- function(ps, groupVar, timeVar, taxonLevel = "Family",
                              top_n = 30, p_col = "padj", p_sig = 0.05,
                              title = NULL, taxaList = NULL, cap_quantile = 1.0, customColors = c("blue","red")) {
  
  if(taxonLevel != "Species"){
    ps <- tax_glom(ps, taxrank = taxonLevel, NArm = FALSE)
  }
  
  res_all <- data.frame() # to store all results
  
  for(timePoint in levels(sample_data(ps)[[timeVar]])){
    
    # Subset to this timepoint
    ps_subset <- prune_samples(sample_data(ps)[[timeVar]] == timePoint, ps)
    ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)
    
    # DESeq2 analysis
    deseq <- phyloseq_to_deseq2(ps_subset, as.formula(paste("~", groupVar)))
    deseq <- DESeq(deseq, test = "Wald", fitType = "parametric", quiet = TRUE)
    res <- results(deseq)
    
    # Convert to data frame and add taxonomy
    sigtab <- cbind(as.data.frame(res),
                    as(tax_table(ps_subset)[rownames(res), ], "matrix"))
    sigtab$padj[is.na(sigtab$padj)] <- 1
    sigtab[[timeVar]] <- timePoint
    
    # rbind to the main table
    res_all <- rbind(res_all, sigtab)
  }
  
  dat <- within(res_all, {
    taxon <- ifelse(is.na(get(taxonLevel)) | get(taxonLevel) == "",
                    "Unclassified", get(taxonLevel))
    tp    <- get(timeVar)
    coef  <- as.numeric(log2FoldChange)
    pval  <- as.numeric(get(p_col)); pval[is.na(pval)] <- 1
    sig   <- ifelse(pval < p_sig, "*", "")
    mag   <- abs(coef)  # size mapped to magnitude of effect
  })
  
  # Select taxons to be displayed
  if(!is.null(taxaList)){
    dat <- dat[dat$taxon %in% taxaList,]
  }
  
  # keep most varying taxa for readability
  keep <- aggregate(mag ~ taxon, dat, mean, na.rm = TRUE)
  keep <- head(keep[order(-keep$mag), "taxon"], top_n)
  dat  <- dat[dat$taxon %in% keep, , drop = FALSE]
  
  # order axes
  ord  <- aggregate(mag ~ taxon, dat, mean, na.rm = TRUE)
  dat$taxon <- factor(dat$taxon, levels = rev(ord[order(ord$mag), "taxon"]))
  if(!is.null(taxaList)){
    dat$taxon <- factor(dat$taxon, levels = rev(taxaList))
  }
  if (is.numeric(dat$tp)) dat$tp <- factor(dat$tp, levels = sort(unique(dat$tp)))
  
  # cap color limits to avoid one extreme dominating; default = max (no cap)
  L <- stats::quantile(abs(dat$coef), probs = cap_quantile, na.rm = TRUE)
  if (!is.finite(L) || L <= 0) L <- 1
  
  p <- ggplot(dat, aes(x = tp, y = taxon)) +
    # size = |log2FC|, fill = signed log2FC
    geom_point(aes(size = mag, fill = coef),
               shape = 21, stroke = 0.2, color = "grey20", alpha = 0.9) +
    # significance star (no p-value color scale)
    geom_text(aes(label = sig), size = 5, color = "white") +
    scale_size_continuous(name = "|log2FC|", range = c(1.5, 12)) +
    scale_fill_gradient2(
      name = "log2FC",
      low = customColors[1], mid = "white", high = customColors[2],
      midpoint = 0, limits = c(-L, L), oob = scales::squish
    ) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  return(list(table = res_all, plot = p))
}




