library(phyloseq)
library(ggplot2)
library(ggtext)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(DESeq2)

# Modified plot_microbiota function from StackbarExtended package, adding 
plot_microbiota_ext <- function (ps_object = ps, exp_group = "group", subset_group = NULL, 
          sample_name = "SampleID", main_level = "Phylum", sub_level = "Family", 
          threshold = 1, n_phy = 4, mean_group = FALSE, hues = c("Oranges", 
                                                                 "Greens", "Blues", "Purples"), color_bias = 2, n_row = 1, 
          n_col = NULL, text_size = 9, legend_size = 7, x_axis_size = 8, 
          differential_analysis = FALSE, mult_comp = FALSE, selected_comparisons = NULL, 
          test = c("Wald", "LRT")[1], fdr_threshold = 0.05, sig_lab = FALSE, 
          fitType = c("parametric", "local", "mean", "glmGamPoi")[1], 
          sfType = c("ratio", "poscounts", "iterate")[1], betaPrior = FALSE, 
          reduced = FALSE, quiet = TRUE, minReplicatesForReplace = 7, 
          modelMatrixType = c("standard", "expanded")[1], useT = FALSE, 
          minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5, parallel = FALSE) 
{
  if (!("phyloseq" %in% class(ps_object))) 
    stop("ps_object must be a phyloseq-class object.")
  if (!main_level %in% colnames(tax_table(ps_object))) 
    stop("main_level column not found in the tax_table.")
  if (!sub_level %in% colnames(tax_table(ps_object))) 
    stop("sub_level column not found in the tax_table.")
  if (!exp_group %in% names(sample_data(ps_object))) 
    stop("exp_group column not found in the sample_data.")
  if (!taxa_are_rows(ps_object)) {
    ps_object <- t(ps_object)
  }
  if (is.null(subset_group) == F) {
    keep_samples = as.character(get_variable(ps_object, 
                                             exp_group)) %in% subset_group
    ps_object = prune_samples(keep_samples, ps_object)
  }
  ps_prop <- transform_sample_counts(ps_object, function(OTU) ((OTU/sum(OTU)) * 
                                                                 100))
  otu <- as.data.frame(otu_table(ps_prop))
  tax <- as.data.frame(tax_table(ps_prop))
  meta <- data.frame(sample_data(ps_prop))
  otu <- otu %>% filter_all(any_vars(. != 0))
  tax <- subset(tax, rownames(tax) %in% rownames(otu))
  ps_f <- tax_glom(ps_prop, taxrank = sub_level, NArm = FALSE)
  otu_f <- as.data.frame(otu_table(ps_f))
  tax_f <- as.data.frame(tax_table(ps_f))
  otu_tax_f <- cbind(tax_f, otu_f)
  position <- ncol(tax) + 1:(ncol(otu_tax_f) - ncol(tax))
  message(paste0("\n", length(position), " samples are analyzed \n"))
  otu_tax_f <- otu_tax_f %>% rowwise() %>% mutate(Mean = mean(c_across(all_of(position)))) %>% 
    mutate(high_abundance = case_when(Mean > threshold ~ 
                                        TRUE, TRUE ~ FALSE))
  topx <- otu_tax_f %>% dplyr::group_by(!!as.name(main_level)) %>% 
    dplyr::summarise(sum_top = sum(c_across(all_of(position)))) %>% 
    dplyr::arrange(desc(sum_top)) %>% dplyr::slice_head(n = n_phy)
  otu_tax_f <- mutate(otu_tax_f, selected_top = !!as.name(main_level) %in% 
                        pull(topx[, main_level]))
  otu_tax_f[is.na(otu_tax_f)] <- "Unknown"
  otu_tax_f <- otu_tax_f %>% mutate(plot_taxa = case_when(high_abundance == 
                                                            TRUE & selected_top == TRUE & !!as.name(sub_level) == 
                                                            "Unknown" ~ paste0("Unknown ", !!as.name(main_level)), 
                                                          high_abundance == TRUE & selected_top == TRUE & !!as.name(sub_level) != 
                                                            "Unknown" ~ !!as.name(sub_level), high_abundance == 
                                                            FALSE & selected_top == TRUE ~ paste0("Others ", 
                                                                                                  !!as.name(main_level)), selected_top == FALSE ~ 
                                                            paste0("Others ")))
  if (nrow(topx) != length(hues)) {
    message("the number of colors choosen (", length(hues), 
            ") is different from the defined number of features to plot (", 
            nrow(topx), ")")
  }
  df <- as_tibble(matrix(nrow = 0, ncol = length(colnames(otu_tax_f))), 
                  .name_repair = ~colnames(otu_tax_f))
  main_level_col <- c()
  i <- 1
  for (i in 1:nrow(topx)) {
    top_loop <- pull(topx[, main_level])[i]
    temp <- otu_tax_f %>% dplyr::filter(!!as.name(main_level) == 
                                          top_loop) %>% dplyr::arrange(desc(Mean))
    message("Among the ", main_level, " ", top_loop, " ", 
            length(unique(temp$plot_taxa)), " features at ", 
            sub_level, " level will be plotted")
    getPalette = colorRampPalette(brewer.pal(length(unique(temp$plot_taxa)), 
                                             hues[i]), bias = color_bias)
    col <- getPalette(length(unique(temp$plot_taxa)) + 1)
    if (length(col) >= 3) {
      main_level_col[i] <- col[length(col) - 1]
    }
    else {
      main_level_col[i] <- col[length(col)]
    }
    u <- 1
    for (u in 1:length(unique(temp$plot_taxa))) {
      print(unique(temp$plot_taxa)[u])
      temp[which(temp$plot_taxa == unique(temp$plot_taxa)[u]), 
           "MyColors"] <- col[u + 1]
    }
    df <- rbind(df, temp)
  }
  unselescted <- otu_tax_f %>% filter(!(!!as.name(main_level) %in% 
                                          pull(topx[, main_level])))
  unselescted$MyColors <- "#000000"
  df <- rbind(df, unselescted)
  df <- df %>% ungroup %>% dplyr::arrange(desc(selected_top), 
                                          !!as.name(main_level), desc(Mean))
  nrow(df) == nrow(otu_tax_f)
  df$plot_taxa <- factor(df$plot_taxa, levels = unique(df$plot_taxa))
  df_long <- melt(df, id = c("plot_taxa", "MyColors", main_level), 
                  measure.vars = meta[, sample_name], variable.name = sample_name)
  df_long <- left_join(df_long, meta, by = sample_name)
  df_long[, main_level] <- ifelse(df_long[, main_level] %in% 
                                    pull(topx[, main_level]), df_long[, main_level], "Others ")
  if (differential_analysis && length(unique(meta[[exp_group]])) != 
      2) {
    print("The number of experimental group to test is not equal to 2, no statistical significance will appear in the legend")
    mult_comp <- T
  }
  if (differential_analysis && length(unique(meta[[exp_group]])) == 
      2) {
    fam_glom <- tax_glom(ps_object, taxrank = sub_level)
    if (sum(rowSums(otu_table(fam_glom)) == 0) > 0) {
      fam_glom <- prune_taxa(rowSums(otu_table(fam_glom)) > 
                               0, fam_glom)
    }
    diagdds = phyloseq_to_deseq2(fam_glom, formula(paste("~", 
                                                         exp_group)))
    if (test == "Wald") {
      diag = DESeq(diagdds, test = test, fitType = fitType, 
                   sfType = sfType, betaPrior = betaPrior, quiet = quiet, 
                   minReplicatesForReplace = minReplicatesForReplace, 
                   modelMatrixType = modelMatrixType, useT = useT, 
                   minmu = minmu, parallel = parallel)
    }
    else {
      diag = DESeq(diagdds, test = test, fitType = fitType, 
                   sfType = sfType, betaPrior = betaPrior, reduced = formula(reduced), 
                   quiet = quiet, minReplicatesForReplace = minReplicatesForReplace, 
                   modelMatrixType = modelMatrixType, useT = useT, 
                   minmu = minmu, parallel = parallel)
    }
    results <- results(diag)
    significant_features <- subset(results, padj < fdr_threshold)
    if (nrow(significant_features) == 0) {
      message("no significant feature detected at ", sub_level, 
              " level.")
    }
    else {
      significant_features <- as.data.frame(cbind(significant_features, 
                                                  tax_table(fam_glom)[rownames(tax_table(fam_glom)) %in% 
                                                                        rownames(significant_features)]))
      significant_features_sub <- significant_features
      df_long$differential_abundance <- FALSE
      df_long$differential_abundance[df_long$plot_taxa %in% 
                                       significant_features[, sub_level]] <- TRUE
      significant_features$stars <- ""
      if (sig_lab == T) {
        significant_features$stars <- symnum(significant_features$padj, 
                                             symbols = c("***", "**", "*", ""), cutpoints = c(0, 
                                                                                              0.001, 0.01, 0.05, 1), corr = FALSE)
        star_vec <- significant_features$stars[match(df_long$plot_taxa, 
                                                     significant_features[, sub_level])]
        star_vec[is.na(star_vec)] <- ""
        df_long$plot_taxa <- paste0(df_long$plot_taxa, 
                                    " ", star_vec)
      }
      df_long$legend_label <- ifelse(df_long$differential_abundance, 
                                     paste0("<b>", df_long$plot_taxa, "</b>"), as.character(df_long$plot_taxa))
      df_long$plot_taxa <- df_long$legend_label
      df_long$plot_taxa <- factor(df_long$plot_taxa, levels = unique(df_long$plot_taxa))
    }
  }
  if (differential_analysis && length(unique(meta[[exp_group]])) == 
      2) {
    main_glom <- tax_glom(ps_object, taxrank = main_level)
    if (sum(rowSums(otu_table(main_glom)) == 0) > 0) {
      main_glom <- prune_taxa(rowSums(otu_table(main_glom)) > 
                                0, main_glom)
    }
    diagdds_main = phyloseq_to_deseq2(main_glom, formula(paste("~", 
                                                               exp_group)))
    if (test == "Wald") {
      diag_main = DESeq(diagdds_main, test = test, fitType = fitType, 
                        sfType = sfType, betaPrior = betaPrior, quiet = quiet, 
                        minReplicatesForReplace = minReplicatesForReplace, 
                        modelMatrixType = modelMatrixType, useT = useT, 
                        minmu = minmu, parallel = parallel)
    }
    else {
      diag_main = DESeq(diagdds_main, test = test, fitType = fitType, 
                        sfType = sfType, betaPrior = betaPrior, reduced = formula(reduced), 
                        quiet = quiet, minReplicatesForReplace = minReplicatesForReplace, 
                        modelMatrixType = modelMatrixType, useT = useT, 
                        minmu = minmu, parallel = parallel)
    }
    results_main <- results(diag_main)
    significant_features_main <- subset(results_main, padj < 
                                          fdr_threshold)
    if (nrow(significant_features_main) == 0) {
      message("no significant feature detected at ", main_level, 
              " level.")
    }
    else {
      significant_features_main <- as.data.frame(cbind(significant_features_main, 
                                                       tax_table(main_glom)[rownames(tax_table(main_glom)) %in% 
                                                                              rownames(significant_features_main)]))
      df_long$differential_abundance_main <- FALSE
      df_long$differential_abundance_main[df_long[, main_level] %in% 
                                            significant_features_main[, main_level]] <- TRUE
      df_long$legend_label_main <- df_long[, main_level]
      significant_features_main$stars <- ""
      if (sig_lab == T) {
        significant_features_main$stars <- symnum(significant_features_main$padj, 
                                                  symbols = c("***", "**", "*", ""), cutpoints = c(0, 
                                                                                                   0.001, 0.01, 0.05, 1), corr = FALSE)
        star_vec_main <- significant_features_main$stars[match(df_long[, 
                                                                       main_level], significant_features_main[, main_level])]
        star_vec_main[is.na(star_vec_main)] <- ""
        df_long[, main_level] <- paste0(df_long[, main_level], 
                                        " ", star_vec_main)
        significant_features_main <- significant_features_main[, 
                                                               -c(ncol(significant_features_main))]
      }
      df_long[, main_level] <- ifelse(df_long$differential_abundance_main, 
                                      paste0("<b>", df_long[, main_level], "</b>"), 
                                      as.character(df_long[, main_level]))
      df_long[, main_level] <- factor(df_long[, main_level], 
                                      levels = unique(df_long[, main_level]))
    }
  }
  if (differential_analysis && mult_comp == T) {
    fam_glom <- tax_glom(ps_object, taxrank = sub_level)
    if (sum(rowSums(otu_table(fam_glom)) == 0) > 0) {
      fam_glom <- prune_taxa(rowSums(otu_table(fam_glom)) > 
                               0, fam_glom)
    }
    diagdds = phyloseq_to_deseq2(fam_glom, formula(paste("~", 
                                                         exp_group)))
    if (test == "Wald") {
      diag = DESeq(diagdds, test = test, fitType = fitType, 
                   sfType = sfType, betaPrior = betaPrior, quiet = quiet, 
                   minReplicatesForReplace = minReplicatesForReplace, 
                   modelMatrixType = modelMatrixType, useT = useT, 
                   minmu = minmu, parallel = parallel)
    }
    else {
      diag = DESeq(diagdds, test = test, fitType = fitType, 
                   sfType = sfType, betaPrior = betaPrior, reduced = formula(reduced), 
                   quiet = quiet, minReplicatesForReplace = minReplicatesForReplace, 
                   modelMatrixType = modelMatrixType, useT = useT, 
                   minmu = minmu, parallel = parallel)
    }
    comparisons <- combn(levels(meta[[exp_group]]), 2, simplify = FALSE)
    if (is.null(selected_comparisons) == F) {
      comparisons <- Filter(function(cmp) {
        any(sapply(selected_comparisons, function(sel) all(sel == 
                                                             cmp)))
      }, comparisons)
    }
    comparison_names <- sapply(comparisons, function(cmp) paste(cmp, 
                                                                collapse = "_vs_"))
    results_list <- list()
    results_list <- lapply(comparisons, function(cmp) {
      contrast_vec <- c(exp_group, cmp[1], cmp[2])
      tryCatch({
        res <- results(diag, contrast = contrast_vec)
        return(res)
      }, error = function(e) {
        message("Failed to compute results for comparison: ", 
                paste(cmp, collapse = " vs "), "\nError: ", 
                e$message)
        return(NULL)
      })
    })
    results_list <- setNames(results_list, comparison_names)
    significant_features_sub <- list()
    for (i in seq_along(results_list)) {
      res <- results_list[[i]]
      comparison_name <- names(results_list)[i]
      significant_features <- subset(res, padj < fdr_threshold)
      if (nrow(significant_features) == 0) {
        message("No significant feature detected for comparison ", 
                comparison_name, " at the ", sub_level, " level.")
      }
      else {
        significant_features <- cbind(significant_features, 
                                      tax_table(fam_glom)[rownames(tax_table(fam_glom)) %in% 
                                                            rownames(significant_features), , drop = FALSE])
        significant_features <- as.data.frame(significant_features)
        significant_features_sub[[comparison_name]] <- significant_features
      }
    }
  }
  if (differential_analysis && mult_comp == T) {
    main_glom <- tax_glom(ps_object, taxrank = main_level)
    if (sum(rowSums(otu_table(main_glom)) == 0) > 0) {
      fam_glom <- prune_taxa(rowSums(otu_table(main_glom)) > 
                               0, main_glom)
    }
    diagdds = phyloseq_to_deseq2(main_glom, formula(paste("~", 
                                                          exp_group)))
    if (test == "Wald") {
      diag = DESeq(diagdds, test = test, fitType = fitType, 
                   sfType = sfType, betaPrior = betaPrior, quiet = quiet, 
                   minReplicatesForReplace = minReplicatesForReplace, 
                   modelMatrixType = modelMatrixType, useT = useT, 
                   minmu = minmu, parallel = parallel)
    }
    else {
      diag = DESeq(diagdds, test = test, fitType = fitType, 
                   sfType = sfType, betaPrior = betaPrior, reduced = formula(reduced), 
                   quiet = quiet, minReplicatesForReplace = minReplicatesForReplace, 
                   modelMatrixType = modelMatrixType, useT = useT, 
                   minmu = minmu, parallel = parallel)
    }
    comparisons <- combn(levels(meta[[exp_group]]), 2, simplify = FALSE)
    if (is.null(selected_comparisons) == F) {
      comparisons <- Filter(function(cmp) {
        any(sapply(selected_comparisons, function(sel) all(sel == 
                                                             cmp)))
      }, comparisons)
    }
    comparison_names <- sapply(comparisons, function(cmp) paste(cmp, 
                                                                collapse = "_vs_"))
    results_list <- list()
    results_list <- lapply(comparisons, function(cmp) {
      contrast_vec <- c(exp_group, cmp[1], cmp[2])
      tryCatch({
        res <- results(diag, contrast = contrast_vec)
        return(res)
      }, error = function(e) {
        message("Failed to compute results for comparison: ", 
                paste(cmp, collapse = " vs "), "\nError: ", 
                e$message)
        return(NULL)
      })
    })
    results_list <- setNames(results_list, comparison_names)
    significant_features_main <- list()
    for (i in seq_along(results_list)) {
      res <- results_list[[i]]
      comparison_name <- names(results_list)[i]
      significant_features <- subset(res, padj < fdr_threshold)
      if (nrow(significant_features) == 0) {
        message("No significant feature detected for comparison ", 
                comparison_name, " at the ", main_level, " level.")
      }
      else {
        significant_features <- cbind(significant_features, 
                                      tax_table(main_glom)[rownames(tax_table(main_glom)) %in% 
                                                             rownames(significant_features), , drop = FALSE])
        significant_features <- as.data.frame(significant_features)
        significant_features_main[[comparison_name]] <- significant_features
      }
    }
  }
  MyColors <- df_long$MyColors
  names(MyColors) <- df_long$plot_taxa
  MyColors2 <- unique(df_long$MyColors)
  names(MyColors2) <- unique(df_long$plot_taxa)
  main_level_col[length(main_level_col) + 1] <- "#000000"
  df_long[, main_level] <- factor(df_long[, main_level], levels = unique(df_long[, 
                                                                                 main_level]))
  vec1 <- unique(df_long[, main_level])
  vec2 <- c(pull(topx[, main_level]), paste0("Others"))
  core_text_vec1 <- gsub("(<[^>]*>|\\*|\\s+$)", "", vec1)
  core_text_vec1 <- trimws(core_text_vec1)
  order_index <- match(core_text_vec1, vec2)
  main_level_col <- main_level_col[order_index]
  names(main_level_col) <- as.character(vec1)
  if (mean_group == F) {
    p <- ggplot(df_long, aes(x = !!as.name(sample_name), 
                             y = value, fill = plot_taxa)) + geom_bar(stat = "identity", 
                                                                      width = 0.85) + ylab("Relative abundance (%)\n") + 
      guides(fill = guide_legend(reverse = FALSE, title = sub_level, 
                                 order = 2)) + theme(line = element_line(colour = "black", 
                                                                         linewidth = 0.5), text = element_text(size = 9), 
                                                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                     panel.border = element_blank(), axis.line = element_line(colour = "black", 
                                                                                                              linewidth = 0.5)) + theme_bw() + theme(panel.grid.major = element_blank(), 
                                                                                                                                                     panel.grid.minor = element_blank(), panel.border = element_blank(), 
                                                                                                                                                     legend.position = "right", legend.box = "vertical", 
                                                                                                                                                     axis.text.x = element_text(angle = 45, hjust = 1, 
                                                                                                                                                                                size = x_axis_size), text = element_text(size = text_size), 
                                                                                                                                                     legend.text = element_markdown(size = legend_size)) + 
      geom_bar(aes(alpha = df_long[, main_level]), stat = "identity", 
               show.legend = TRUE) + scale_alpha_manual(values = rep(1, 
                                                                     length(unique(df_long[, main_level]))), guide = guide_legend(order = 1, 
                                                                                                                                  override.aes = list(fill = main_level_col)), name = main_level) + 
      scale_fill_manual("plot_taxa", values = MyColors2)
    p <- p + facet_wrap(~df_long[[exp_group]], scales = "free_x", 
                        nrow = n_row, ncol = n_col)
    p
  }
  else {
    df_long <- df_long %>% group_by(!!as.name(exp_group), 
                                    plot_taxa, !!as.name(main_level)) %>% reframe(n = n(), 
                                                                                  sum = sum(as.double(value)))
    df_long <- data.frame(df_long)
    i <- 1
    for (i in 1:length(unique(df_long[, exp_group]))) {
      df_long[, "sum"][df_long[, exp_group] == unique(df_long[, 
                                                              exp_group])[i]] <- df_long[, "sum"][df_long[, 
                                                                                                          exp_group] == unique(df_long[, exp_group])[i]]/count(meta[, 
                                                                                                                                                                    exp_group] == unique(meta[, exp_group])[i])
    }
    colnames(df_long)[colnames(df_long) == "sum"] <- "value"
    p <- ggplot(df_long, aes(x = !!as.name(exp_group), y = value, 
                             fill = plot_taxa)) + geom_bar(stat = "identity", 
                                                           width = 0.85) + ylab("Relative abundance (%)\n") + 
      guides(fill = guide_legend(reverse = FALSE, title = sub_level, 
                                 order = 2)) + theme(line = element_line(colour = "black", 
                                                                         linewidth = 0.5), text = element_text(size = 9), 
                                                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                     panel.border = element_blank(), axis.line = element_line(colour = "black", 
                                                                                                              linewidth = 0.5)) + theme_bw() + theme(panel.grid.major = element_blank(), 
                                                                                                                                                     panel.grid.minor = element_blank(), panel.border = element_blank(), 
                                                                                                                                                     legend.position = "right", legend.box = "vertical", 
                                                                                                                                                     axis.text.x = element_text(angle = 45, hjust = 1, 
                                                                                                                                                                                size = x_axis_size), text = element_text(size = text_size), 
                                                                                                                                                     legend.text = element_markdown(size = legend_size)) + 
      geom_bar(aes(alpha = df_long[, main_level]), stat = "identity", 
               show.legend = TRUE) + scale_alpha_manual(values = rep(1, 
                                                                     length(unique(df_long[, main_level]))), guide = guide_legend(order = 1, 
                                                                                                                                  override.aes = list(fill = main_level_col)), name = main_level) + 
      scale_fill_manual("plot_taxa", values = MyColors2)
    p <- p + facet_wrap(~df_long[[exp_group]], scales = "free_x", 
                        nrow = n_row, ncol = n_col)
    p
  }
  if (differential_analysis == T) {
    return(list(significant_table_main = significant_features_main, 
                significant_table_sub = significant_features_sub, 
                plot = p, main_names=unique(df_long[[main_level]]), sub_names=unique(df_long$plot_taxa)))
  }
  else {
    return(list(plot = p))
  }
}

# Function to write and save stackbarExtended sig_table
writeStackbarExtendedSigTable <- function(main_table, sub_table, filepath){
  
  # Initialize empty dataframe to append tables to
  table_to_write <- data.frame()
  
  # Iterate over the list of tables for main table
  for (i in seq_along(main_table)){
    
    # Extract the table
    table <- main_table[[i]]
    
    # Add a column with the name of the current table
    table$comparaison <- names(main_table)[i]
    
    # Add col indication if it is from main or sub table
    table$level <- "main"
    
    # Add signifiance symbols column
    table$significance <- cut(table$padj,
                              breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                              labels = c("***", "**", "*", "NS"))
    
    # Append to the final table
    table_to_write <- rbind(table_to_write, table)
    
  }
  
  # Iterate over the list of tables for sub table
  for (i in seq_along(sub_table)){
    
    # Extract the table
    table <- sub_table[[i]]
    
    # Add a column with the name of the current table
    table$comparaison <- names(sub_table)[i]
    
    # Add col indication if it is from main or sub table
    table$level <- "sub"
    
    # Add signifiance symbols column
    table$significance <- cut(table$padj,
                              breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                              labels = c("***", "**", "*", "NS"))
    
    # Append to the final table
    table_to_write <- rbind(table_to_write, table)
    
  }
  
  write_xlsx(x = table_to_write, path = filepath)
  
}

# Function that takes a stats excel file generated by the function above and generates
# heatmap for the pvalues

pvaluesHmap <- function(stats, selected_comparisons,
                        taxons, lvl, txn_lvl, group, path){
  
  stats=stats[stats$level==lvl,] # Table with only stats for main taxa
  stat_hmap=matrix(nrow = length(selected_comparisons), ncol = length(taxons)) # Initiate matrix with comparisons as rows and taxa as cols
  stat_hmap[]=1# Fill matrix with 1 (not significant)
  row.names(stat_hmap) = selected_comparisons
  colnames(stat_hmap) = c(taxons)
  
  # Iterate through dataframe and add p_adjusted values in the matrix
  for(i in 1:nrow(stat_hmap)){
    for(k in 1:ncol(stat_hmap)){
      g = row.names(stat_hmap)[i]
      subdf = stats[stats$comparaison==g,]
      if(length(subdf$padj[subdf[txn_lvl]==colnames(stat_hmap)[k]])>0){
        stat_hmap[i,k]=subdf$padj[subdf[txn_lvl]==colnames(stat_hmap)[k]]
      }
    }
  }
  
  # Reshape the data to long format
  stat_hmap <- melt(stat_hmap)
  colnames(stat_hmap)[1] <- "comparaison"
  colnames(stat_hmap)[2] <- txn_lvl
  stat_hmap$comparaison <- factor(stat_hmap$comparaison, levels = rev(selected_comparisons))
  
  # Plot heatmap with ggplot2
  ggplot(stat_hmap, aes(x = .data[[txn_lvl]], y = comparaison, fill = value)) +
    geom_tile(color = "black", lwd = 0.75, linetype = 1) +
    coord_fixed() + # Makes thing squared
    geom_text(aes(label = ifelse(value=="1","NS",round(value, 5))), color = "gray", size = 4) + # Show significance labels
    scale_fill_gradient(low = "black", high = "white") +
    theme_minimal() +
    labs(x = txn_lvl, y = "Comparison", fill = "") +
    guides(fill = "none")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  
  
}