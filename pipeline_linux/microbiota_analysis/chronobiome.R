# Function that creates a timeline for gut microbiota relative abundance composition

plot_timeline_2_groups <- function(
    ps_object,
    exp_group, # must be as factor
    time_group, # must be as factor
    sample_name,
    main_level = 'Phylum',
    sub_level = 'Family',
    threshold = 1,
    average_relab_per_group = TRUE,
    n_phy,
    hues = c("Oranges", "Greens", "Blues", "Purples"),
    color_bias = 2,
    differential_analysis = FALSE,
    test = c("Wald", "LRT")[1],
    fdr_threshold = 0.05,
    sig_lab = FALSE,
    fitType = c("parametric", "local", "mean", "glmGamPoi")[1],
    sfType = c("ratio", "poscounts", "iterate")[1],
    betaPrior = FALSE,
    reduced = FALSE,
    quiet = TRUE,
    minReplicatesForReplace = 7,
    modelMatrixType = c("standard", "expanded")[1],
    useT = FALSE,
    minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5,
    parallel = FALSE
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

  #### 3. Differential Analysis (Optional) ####
  significant_features_sub <- NULL
  significant_features_main <- NULL
  
  if (differential_analysis) {
    # -- Two-group (if only 2 groups and not using multiple comparisons) --
    if (!mult_comp && length(unique(meta[[exp_group]])) == 2) {
      # Sub-level analysis
      fam_glom <- tax_glom(ps_object, taxrank = sub_level)
      if (sum(rowSums(otu_table(fam_glom)) == 0) > 0) {
        fam_glom <- prune_taxa(rowSums(otu_table(fam_glom)) > 0, fam_glom)
      }
      diagdds <- phyloseq_to_deseq2(fam_glom, formula(paste("~", exp_group)))
      if (test == "Wald") {
        diag <- DESeq(diagdds, test = test, fitType = fitType, sfType = sfType,
                      betaPrior = betaPrior, quiet = quiet,
                      minReplicatesForReplace = minReplicatesForReplace,
                      modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                      parallel = parallel)
      } else {
        diag <- DESeq(diagdds, test = test, fitType = fitType, sfType = sfType,
                      betaPrior = betaPrior, reduced = formula(reduced), quiet = quiet,
                      minReplicatesForReplace = minReplicatesForReplace,
                      modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                      parallel = parallel)
      }
      res <- results(diag)
      significant_features <- subset(res, padj < fdr_threshold)
      if (nrow(significant_features) == 0) {
        message("No significant feature detected at ", sub_level, " level.")
      } else {
        significant_features <- as.data.frame(cbind(significant_features,
                                                    tax_table(fam_glom)[rownames(tax_table(fam_glom)) %in%
                                                                          rownames(significant_features), ]))
        significant_features_sub <- significant_features
        df_long$differential_abundance <- FALSE
        df_long$differential_abundance[df_long$plot_taxa %in% significant_features[, sub_level]] <- TRUE
        significant_features$stars <- ""
        if (sig_lab == TRUE) {
          significant_features$stars <- symnum(significant_features$padj,
                                               symbols = c("***", "**", "*", ""),
                                               cutpoints = c(0, .001, .01, .05, 1),
                                               corr = FALSE)
          star_vec <- significant_features$stars[match(df_long$plot_taxa, significant_features[, sub_level])]
          star_vec[is.na(star_vec)] <- ""
          df_long$plot_taxa <- paste0(df_long$plot_taxa, " ", star_vec)
        }
        df_long$legend_label <- ifelse(df_long$differential_abundance,
                                       paste0("<b>", df_long$plot_taxa, "</b>"),
                                       as.character(df_long$plot_taxa))
        df_long$plot_taxa <- df_long$legend_label
        df_long$plot_taxa <- factor(df_long$plot_taxa, levels = unique(df_long$plot_taxa))
      }
      
      # Main-level analysis
      main_glom <- tax_glom(ps_object, taxrank = main_level)
      if (sum(rowSums(otu_table(main_glom)) == 0) > 0) {
        main_glom <- prune_taxa(rowSums(otu_table(main_glom)) > 0, main_glom)
      }
      diagdds_main <- phyloseq_to_deseq2(main_glom, formula(paste("~", exp_group)))
      if (test == "Wald") {
        diag_main <- DESeq(diagdds_main, test = test, fitType = fitType, sfType = sfType,
                           betaPrior = betaPrior, quiet = quiet,
                           minReplicatesForReplace = minReplicatesForReplace,
                           modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                           parallel = parallel)
      } else {
        diag_main <- DESeq(diagdds_main, test = test, fitType = fitType, sfType = sfType,
                           betaPrior = betaPrior, reduced = formula(reduced), quiet = quiet,
                           minReplicatesForReplace = minReplicatesForReplace,
                           modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                           parallel = parallel)
      }
      res_main <- results(diag_main)
      significant_features_main <- subset(res_main, padj < fdr_threshold)
      if (nrow(significant_features_main) == 0) {
        message("No significant feature detected at ", main_level, " level.")
      } else {
        significant_features_main <- as.data.frame(cbind(significant_features_main,
                                                         tax_table(main_glom)[rownames(tax_table(main_glom)) %in%
                                                                                rownames(significant_features_main), ]))
        df_long$differential_abundance_main <- FALSE
        df_long$differential_abundance_main[df_long[, main_level] %in% significant_features_main[, main_level]] <- TRUE
        df_long$legend_label_main <- df_long[, main_level]
        significant_features_main$stars <- ""
        if (sig_lab == TRUE) {
          significant_features_main$stars <- symnum(significant_features_main$padj,
                                                    symbols = c("***", "**", "*", ""),
                                                    cutpoints = c(0, .001, .01, .05, 1),
                                                    corr = FALSE)
          star_vec_main <- significant_features_main$stars[match(df_long[, main_level],
                                                                 significant_features_main[, main_level])]
          star_vec_main[is.na(star_vec_main)] <- ""
          df_long[, main_level] <- paste0(df_long[, main_level], " ", star_vec_main)
          significant_features_main <- significant_features_main[, -ncol(significant_features_main)]
        }
        df_long[, main_level] <- ifelse(df_long$differential_abundance_main,
                                        paste0("<b>", df_long[, main_level], "</b>"),
                                        as.character(df_long[, main_level]))
        df_long[, main_level] <- factor(df_long[, main_level], levels = unique(df_long[, main_level]))
      }
    }
    
    # -- Multiple Comparisons Differential Analysis (if mult_comp is TRUE) --
    if (mult_comp) {
      # Sub-level multiple comparisons
      fam_glom <- tax_glom(ps_object, taxrank = sub_level)
      if (sum(rowSums(otu_table(fam_glom)) == 0) > 0) {
        fam_glom <- prune_taxa(rowSums(otu_table(fam_glom)) > 0, fam_glom)
      }
      if (twoFactor) {
        diagdds <- phyloseq_to_deseq2(fam_glom, formula(paste("~", fac1, "+", fac2, "+", fac1, ":", fac2)))
        colData(diagdds)[[fac1]] <- relevel(colData(diagdds)[[fac1]], ref = refFac1)
        colData(diagdds)[[fac2]] <- relevel(colData(diagdds)[[fac2]], ref = refFac2)
      } else {
        diagdds <- phyloseq_to_deseq2(fam_glom, formula(paste("~", exp_group)))
      }
      if (test == "Wald") {
        diag <- DESeq(diagdds, test = test, fitType = fitType, sfType = sfType,
                      betaPrior = betaPrior, quiet = quiet,
                      minReplicatesForReplace = minReplicatesForReplace,
                      modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                      parallel = parallel)
      } else {
        diag <- DESeq(diagdds, test = test, fitType = fitType, sfType = sfType,
                      betaPrior = betaPrior, reduced = formula(reduced), quiet = quiet,
                      minReplicatesForReplace = minReplicatesForReplace,
                      modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                      parallel = parallel)
      }
      comparisons <- combn(levels(meta[[exp_group]]), 2, simplify = FALSE)
      if (!is.null(selected_comparisons)) {
        comparisons <- Filter(function(cmp) {
          any(sapply(selected_comparisons, function(sel) all(sel == cmp)))
        }, comparisons)
      }
      comparison_names <- sapply(comparisons, function(cmp) paste(cmp, collapse = "_vs_"))
      results_list <- list()
      if (twoFactor) {
        res_subset1 <- results(diag, contrast = list(c(resultsNames(diag)[3])))
        res_subset2 <- results(diag, contrast = list(c(resultsNames(diag)[3], resultsNames(diag)[4])))
        res_subset3 <- results(diag, contrast = list(c(resultsNames(diag)[2])))
        res_subset4 <- results(diag, contrast = list(c(resultsNames(diag)[2], resultsNames(diag)[4])))
        results_list <- append(results_list, list(res_subset1, res_subset2, res_subset3, res_subset4))
      } else {
        results_list <- lapply(comparisons, function(cmp) {
          contrast_vec <- c(exp_group, cmp[1], cmp[2])
          tryCatch({
            res <- results(diag, contrast = contrast_vec)
            return(res)
          }, error = function(e) {
            message("Failed to compute results for comparison: ", paste(cmp, collapse = " vs "),
                    "\nError: ", e$message)
            return(NULL)
          })
        })
      }
      results_list <- setNames(results_list, comparison_names)
      significant_features_sub <- list()
      for (i in seq_along(results_list)) {
        res <- results_list[[i]]
        comparison_name <- names(results_list)[i]
        significant_features <- subset(res, padj < fdr_threshold)
        if (nrow(significant_features) == 0) {
          message("No significant feature detected for comparison ", comparison_name,
                  " at the ", sub_level, " level.")
        } else {
          significant_features <- cbind(significant_features,
                                        tax_table(fam_glom)[rownames(tax_table(fam_glom)) %in%
                                                              rownames(significant_features), , drop = FALSE])
          significant_features <- as.data.frame(significant_features)
          significant_features_sub[[comparison_name]] <- significant_features
        }
      }
      
      # Main-level multiple comparisons
      main_glom <- tax_glom(ps_object, taxrank = main_level)
      if (sum(rowSums(otu_table(main_glom)) == 0) > 0) {
        main_glom <- prune_taxa(rowSums(otu_table(main_glom)) > 0, main_glom)
      }
      if (twoFactor) {
        diagdds <- phyloseq_to_deseq2(main_glom, formula(paste("~", fac1, "+", fac2, "+", fac1, ":", fac2)))
        colData(diagdds)[[fac1]] <- relevel(colData(diagdds)[[fac1]], ref = refFac1)
        colData(diagdds)[[fac2]] <- relevel(colData(diagdds)[[fac2]], ref = refFac2)
      } else {
        diagdds <- phyloseq_to_deseq2(main_glom, formula(paste("~", exp_group)))
      }
      if (test == "Wald") {
        diag <- DESeq(diagdds, test = test, fitType = fitType, sfType = sfType,
                      betaPrior = betaPrior, quiet = quiet,
                      minReplicatesForReplace = minReplicatesForReplace,
                      modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                      parallel = parallel)
      } else {
        diag <- DESeq(diagdds, test = test, fitType = fitType, sfType = sfType,
                      betaPrior = betaPrior, reduced = formula(reduced), quiet = quiet,
                      minReplicatesForReplace = minReplicatesForReplace,
                      modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                      parallel = parallel)
      }
      comparisons <- combn(levels(meta[[exp_group]]), 2, simplify = FALSE)
      if (!is.null(selected_comparisons)) {
        comparisons <- Filter(function(cmp) {
          any(sapply(selected_comparisons, function(sel) all(sel == cmp)))
        }, comparisons)
      }
      comparison_names <- sapply(comparisons, function(cmp) paste(cmp, collapse = "_vs_"))
      results_list <- list()
      if (twoFactor) {
        res_subset1 <- results(diag, contrast = list(c(resultsNames(diag)[3])))
        res_subset2 <- results(diag, contrast = list(c(resultsNames(diag)[3], resultsNames(diag)[4])))
        res_subset3 <- results(diag, contrast = list(c(resultsNames(diag)[2])))
        res_subset4 <- results(diag, contrast = list(c(resultsNames(diag)[2], resultsNames(diag)[4])))
        results_list <- append(results_list, list(res_subset1, res_subset2, res_subset3, res_subset4))
      } else {
        results_list <- lapply(comparisons, function(cmp) {
          contrast_vec <- c(exp_group, cmp[1], cmp[2])
          tryCatch({
            res <- results(diag, contrast = contrast_vec)
            return(res)
          }, error = function(e) {
            message("Failed to compute results for comparison: ", paste(cmp, collapse = " vs "),
                    "\nError: ", e$message)
            return(NULL)
          })
        })
      }
      results_list <- setNames(results_list, comparison_names)
      significant_features_main <- list()
      for (i in seq_along(results_list)) {
        res <- results_list[[i]]
        comparison_name <- names(results_list)[i]
        significant_features <- subset(res, padj < fdr_threshold)
        if (nrow(significant_features) == 0) {
          message("No significant feature detected for comparison ", comparison_name,
                  " at the ", main_level, " level.")
        } else {
          significant_features <- cbind(significant_features,
                                        tax_table(main_glom)[rownames(tax_table(main_glom)) %in%
                                                               rownames(significant_features), , drop = FALSE])
          significant_features <- as.data.frame(significant_features)
          significant_features_main[[comparison_name]] <- significant_features
        }
      }
    }
    
    # -- Timepoints-Specific Differential Analysis --
    if (timePoints) {
      ## Sub-level analysis per timepoint
      fam_glom <- tax_glom(ps_object, taxrank = sub_level)
      if (sum(rowSums(otu_table(fam_glom)) == 0) > 0) {
        fam_glom <- prune_taxa(rowSums(otu_table(fam_glom)) > 0, fam_glom)
      }
      comparisons <- combn(levels(meta[[combined_group]]), 2, simplify = FALSE)
      if (!is.null(selected_comparisons)) {
        comparisons <- Filter(function(cmp) {
          any(sapply(selected_comparisons, function(sel) all(sel == cmp)))
        }, comparisons)
      }
      comparison_names <- sapply(comparisons, function(cmp) paste(cmp, collapse = "_vs_"))
      results_list <- list()
      for(i in seq_along(levels(sample_data(ps_object)[[time_variable]]))) {
        fam_subset <- prune_samples(sample_data(ps_object)[[time_variable]] == 
                                      levels(sample_data(ps_object)[[time_variable]])[i], fam_glom)
        diagdds <- phyloseq_to_deseq2(fam_subset, formula(paste("~", exp_group)))
        if (test == "Wald") {
          diag <- DESeq(diagdds, test = test, fitType = fitType, sfType = sfType,
                        betaPrior = betaPrior, quiet = quiet,
                        minReplicatesForReplace = minReplicatesForReplace,
                        modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                        parallel = parallel)
        } else {
          diag <- DESeq(diagdds, test = test, fitType = fitType, sfType = sfType,
                        betaPrior = betaPrior, reduced = formula(reduced), quiet = quiet,
                        minReplicatesForReplace = minReplicatesForReplace,
                        modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                        parallel = parallel)
        }
        results_list[[comparison_names[i]]] <- results(diag, contrast = list(c(resultsNames(diag)[2])))
      }
      significant_features_sub <- list()
      for (i in seq_along(results_list)) {
        res <- results_list[[i]]
        comparison_name <- names(results_list)[i]
        significant_features <- subset(res, padj < fdr_threshold)
        if (nrow(significant_features) == 0) {
          message("No significant feature detected for comparison ", comparison_name,
                  " at the ", sub_level, " level (timepoint).")
        } else {
          significant_features <- cbind(significant_features,
                                        tax_table(fam_subset)[rownames(tax_table(fam_subset)) %in%
                                                                rownames(significant_features), , drop = FALSE])
          significant_features <- as.data.frame(significant_features)
          significant_features_sub[[comparison_name]] <- significant_features
        }
      }
      
      ## Main-level analysis per timepoint
      main_glom <- tax_glom(ps_object, taxrank = main_level)
      if (sum(rowSums(otu_table(main_glom)) == 0) > 0) {
        main_glom <- prune_taxa(rowSums(otu_table(main_glom)) > 0, main_glom)
      }
      comparisons <- combn(levels(meta[[combined_group]]), 2, simplify = FALSE)
      if (!is.null(selected_comparisons)) {
        comparisons <- Filter(function(cmp) {
          any(sapply(selected_comparisons, function(sel) all(sel == cmp)))
        }, comparisons)
      }
      comparison_names <- sapply(comparisons, function(cmp) paste(cmp, collapse = "_vs_"))
      results_list <- list()
      for(i in seq_along(levels(sample_data(ps_object)[[time_variable]]))) {
        main_subset <- prune_samples(sample_data(ps_object)[[time_variable]] == 
                                       levels(sample_data(ps_object)[[time_variable]])[i], main_glom)
        diagdds <- phyloseq_to_deseq2(main_subset, formula(paste("~", exp_group)))
        if (test == "Wald") {
          diag <- DESeq(diagdds, test = test, fitType = fitType, sfType = sfType,
                        betaPrior = betaPrior, quiet = quiet,
                        minReplicatesForReplace = minReplicatesForReplace,
                        modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                        parallel = parallel)
        } else {
          diag <- DESeq(diagdds, test = test, fitType = fitType, sfType = sfType,
                        betaPrior = betaPrior, reduced = formula(reduced), quiet = quiet,
                        minReplicatesForReplace = minReplicatesForReplace,
                        modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                        parallel = parallel)
        }
        results_list[[comparison_names[i]]] <- results(diag, contrast = list(c(resultsNames(diag)[2])))
      }
      significant_features_main <- list()
      for (i in seq_along(results_list)) {
        res <- results_list[[i]]
        comparison_name <- names(results_list)[i]
        significant_features <- subset(res, padj < fdr_threshold)
        if (nrow(significant_features) == 0) {
          message("No significant feature detected for comparison ", comparison_name,
                  " at the ", main_level, " level (timepoint).")
        } else {
          significant_features <- cbind(significant_features,
                                        tax_table(main_subset)[rownames(tax_table(main_subset)) %in%
                                                                 rownames(significant_features), , drop = FALSE])
          significant_features <- as.data.frame(significant_features)
          significant_features_main[[comparison_name]] <- significant_features
        }
      }
    }
  }
  
  
  #### 4. Prepare Colors for Plotting ####
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
  
  #### 5. Plotting ####
  if (average_relab_per_group) {
    df_long <- merge(df_long, meta, by = "combined_group") # Merge metadata to retrieve the timepoint and group column
    View(df_long)
    
    # Ensure your time variable is numeric (so geom_area draws it as continuous)
    df_long[[time_group]] <- as.numeric(as.character(df_long[[time_group]]))
    
    # Make sure plot_taxa is a factor in the stacking order you want
    df_long$plot_taxa <- factor(df_long$plot_taxa,
                                levels = unique(df_long$plot_taxa))
    
    # Build the ggplot
    p <- ggplot(df_long,
                aes(x = .data[[time_group]],
                    y = value,
                    fill = plot_taxa,
                    group = plot_taxa)) +
      geom_area() +
      # use your named‐vector of colours (plot_taxa → hex) that you built earlier
      scale_fill_manual(name = sub_level, values = MyColors2) +
      # one facet per experimental group (e.g. diet)
      facet_wrap(as.formula(paste("~", exp_group)), ncol = 1) +
      labs(x    = time_group,
           y    = "Relative abundance (%)",
           fill = sub_level) +
      theme_minimal()
    
    print(p)
    
    
  }
    
  
  
}
