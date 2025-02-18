library(phyloseq)
library(ggplot2)
library(ggtext)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(DESeq2)

# For multiple factors design (fac1 + fac2 + fac1:fac2). For now can only handle differences
# with 2 groups per factor
plot_microbiota_2Fac <- function(ps_object = ps,
                            exp_group = 'group',
                            subset_group = NULL,
                            twoFactor = FALSE,
                            fac1 = NULL,
                            refFac1 = NULL,
                            fac2 = NULL,
                            refFac2 = NULL,
                            sample_name = 'SampleID',
                            main_level = 'Phylum',
                            sub_level = 'Family',
                            threshold = 1,
                            n_phy = 4,
                            mean_group = FALSE,
                            hues = c("Oranges", "Greens", "Blues", "Purples"),
                            color_bias = 2,
                            n_row = 1,
                            n_col = NULL,
                            text_size = 9,
                            legend_size = 7,
                            x_axis_size = 8,
                            differential_analysis = FALSE,
                            mult_comp = FALSE,
                            selected_comparisons = NULL,
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
                            parallel = FALSE,
                            showOnlySubLegend = FALSE
) {
  
  
  
  # Validate inputs
  if (!("phyloseq" %in% class(ps_object)))
    stop("ps_object must be a phyloseq-class object.")
  if (!main_level %in% colnames(tax_table(ps_object)))
    stop("main_level column not found in the tax_table.")
  if (!sub_level %in% colnames(tax_table(ps_object)))
    stop("sub_level column not found in the tax_table.")
  if (!exp_group %in% names(sample_data(ps_object)))
    stop("exp_group column not found in the sample_data.")
  
  
  #Assure taxa are rows
  if (!taxa_are_rows(ps_object)) {
    ps_object <- t(ps_object)
  }
  
  #Subset the samples belonging to the groups selected
  if (is.null(subset_group) == F) {
    keep_samples = as.character(get_variable(ps_object, exp_group)) %in% subset_group
    ps_object = prune_samples(keep_samples, ps_object)
    
  }
  
  #transform to relative abundance
  ps_prop <-
    transform_sample_counts(ps_object, function(OTU)
      ((OTU / sum(OTU)) * 100))
  
  
  #store TAX, OTU and metadata tables
  otu <-  as.data.frame(otu_table(ps_prop))
  tax <-  as.data.frame(tax_table(ps_prop))
  meta <- data.frame(sample_data(ps_prop))
  
  #clean up OTU absent and taxa not present
  
  otu <- otu %>% filter_all(any_vars(. != 0))
  
  tax <- subset(tax, rownames(tax) %in% rownames(otu))
  
  #Agglomerate at the sub level defined
  ps_f <- tax_glom(ps_prop, taxrank = sub_level , NArm = FALSE)
  #Store OTU and TAX
  otu_f <-  as.data.frame(otu_table(ps_f))
  tax_f <-  as.data.frame(tax_table(ps_f))
  
  #Bind TAX and OTU
  
  otu_tax_f <- cbind(tax_f, otu_f)
  
  print(otu_tax_f[is.na(otu_tax_f$Phylum),])
  
  
  #create a vector, if True, the ASV mean among the samples is higher than the threshold
  
  #keep the position of the column storing the abundance of the taxa in the otu_tax_f object
  position <- ncol(tax) + 1:(ncol(otu_tax_f) - ncol(tax))
  
  message(paste0('\n', length(position), ' samples are analyzed \n'))
  
  #Create 2 vectors in the dataframe :
  #1 storing the abundance mean of the taxa
  #2 T or F this taxa is above the define threshold of filtering
  
  otu_tax_f <- otu_tax_f %>%
    rowwise() %>%
    mutate(Mean = mean(c_across(all_of(position)))) %>%
    mutate(high_abundance = case_when(Mean > threshold  ~ TRUE,
                                      TRUE  ~ FALSE))
  
  print(otu_tax_f$Phylum)
  
  #Store top x of the main_level phylogeny
  
  topx <- otu_tax_f %>%
    dplyr::group_by(!!as.name(main_level)) %>%
    dplyr::summarise(sum_top = sum(c_across(all_of(position)))) %>%
    dplyr::arrange(desc(sum_top)) %>%
    dplyr::slice_head(n = n_phy)
  
  
  #Create a vector, TRUE when the phylum belongs to the top X, the ones we want to plot
  otu_tax_f <- mutate(otu_tax_f,
                      selected_top = !!as.name(main_level) %in% pull(topx[, main_level]))
  
  #replace NA by "Unknown" 
  otu_tax_f[is.na(otu_tax_f)] <- "Unknown"
  
  #Add a column 'plot_taxa' to create the legend of the graphs :
  
  otu_tax_f <- otu_tax_f %>%
    mutate(
      plot_taxa = case_when(
        high_abundance == TRUE  &
          selected_top == TRUE &
          !!as.name(sub_level) == "Unknown"  ~
          paste0("Unknown ",!!as.name(main_level)),
        high_abundance == TRUE  &
          selected_top == TRUE  &
          !!as.name(sub_level) != "Unknown"  ~ !!as.name(sub_level),
        high_abundance == FALSE &
          selected_top == TRUE  ~ paste0("Others ", !!as.name(main_level)),
        selected_top   == FALSE ~ paste0("Others ")
      )
    )
  
  
  if (nrow(topx) != length(hues)) {
    message(
      'the number of colors choosen (',
      length(hues),
      ') is different from the defined number of features to plot (',
      nrow(topx),
      ')'
    )
  }
  
  print(topx)
  return(NULL)
  
  
  #initialize 
  df <-
    as_tibble(matrix(nrow = 0, ncol = length(colnames(otu_tax_f))),
              .name_repair = ~ colnames(otu_tax_f))
  main_level_col <- c()
  
  #loop through selected main_level to add color
  i <- 1
  
  for (i in 1:nrow(topx)) {
    top_loop <- pull(topx[, main_level])[i]
    # print(top_loop)
    
    # Subset the df with unique features and in function of the Means abundance
    
    temp <- otu_tax_f %>%
      dplyr::filter(!!as.name(main_level) == top_loop) %>%
      dplyr::arrange(desc(Mean))
    
    message(
      'Among the ',
      main_level,
      ' ' ,
      top_loop,
      ' ' ,
      length(unique(temp$plot_taxa)) ,
      ' features at ',
      sub_level,
      ' level will be plotted'
    )
    
    #define the color palette
    getPalette = colorRampPalette(brewer.pal(length(unique(temp$plot_taxa)), hues[i]), bias = color_bias)
    col <-
      getPalette(length(unique(temp$plot_taxa)) + 1)
    
    #store the darkestcolor of the palette to plot main_level legend graph
    if (length(col) >= 3) {
      main_level_col[i] <-   col[length(col)-1]
    } else {
      main_level_col[i] <-   col[length(col)]
    }  
    
    #loop among the sub_level feature and add color to it
    u <- 1
    for (u in 1:length(unique(temp$plot_taxa))) {
      print(unique(temp$plot_taxa)[u])
      
      #Add the color
      temp[which(temp$plot_taxa ==  unique(temp$plot_taxa)[u]), 'MyColors'] <-
        col[u + 1]
      
    }
    
    #create the dataframe storing the colors information
    df <- rbind(df, temp)
    
  }
  
  
  #Select features that don't belong to high abundance levels and assigned the black color to them
  unselescted <- otu_tax_f %>%
    filter(!(!!as.name(main_level) %in% pull(topx[, main_level])))
  
  unselescted$MyColors <- '#000000'
  
  #Merge the 2 dataframes
  df <- rbind(df, unselescted)
  
  #Order df to define the order that will be plotted, first the selected main_level, then the main_level in alphabetical order, then the mean abundance
  
  df <- df %>% ungroup %>%
    dplyr::arrange(desc(selected_top), !!as.name(main_level) , desc(Mean))
  
  
  #check if no features have been lost
  nrow(df) == nrow(otu_tax_f)
  
  #factor to keep order in the plot
  df$plot_taxa <-
    factor(df$plot_taxa, levels = unique(df$plot_taxa))
  
  
  #Transform df into long format by sample
  df_long <- melt(
    df,
    id = c("plot_taxa", "MyColors", main_level),
    measure.vars = meta[, sample_name],
    variable.name = sample_name
  )
  
  
  #Add exp_group to df_long
  df_long <- left_join(df_long, meta, by = sample_name)
  
  
  #Replace low abundance features at level X by "Others "
  df_long[,main_level] <- ifelse(df_long[,main_level] %in% pull(topx[, main_level]), 
                                 df_long[,main_level], 
                                 "Others ")
  
  
  if (differential_analysis &&
      length(unique(meta[[exp_group]])) != 2) {
    print(
      "The number of experimental group to test is not equal to 2, no statistical significance will appear in the legend"
    )
    mult_comp <- T
    
  }
  
  
  
  # Differential abundance analysis on sub_level using DESeq2 : 2 groups only
  if (differential_analysis &&
      length(unique(meta[[exp_group]])) == 2) {
    fam_glom <- tax_glom(ps_object, taxrank = sub_level)
    
    #remove OTU with 0 count across the dataset
    
    if (sum(rowSums(otu_table(fam_glom)) == 0) > 0) {
      fam_glom <-
        prune_taxa(rowSums(otu_table(fam_glom)) > 0, fam_glom)
      
    }
    
    diagdds = phyloseq_to_deseq2(fam_glom,  formula(paste("~", exp_group)))
    # Run DESeq2 analysis
    if (test == "Wald"){
      
      
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    } else {
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, reduced = formula(reduced), quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    }
    
    # Get the differentially abundant features
    results <- results(diag)
    
    # Extract features with adjusted p-value below the threshold
    significant_features <- subset(results, padj < fdr_threshold)
    
    #Check if significant features were detected
    if (nrow(significant_features) == 0) {
      message("no significant feature detected at ", sub_level, " level.")
    } else {
      #We have significantly different taxa, we link them to their sub_level
      
      
      significant_features <-
        as.data.frame(cbind(significant_features,
                            tax_table(fam_glom)[rownames(tax_table(fam_glom)) %in% rownames(significant_features)]))
      
      #store the sig table to return it 
      significant_features_sub <- significant_features
      
      # Add information to df_long to mark differentially abundant features
      df_long$differential_abundance <- FALSE
      df_long$differential_abundance[df_long$plot_taxa %in% significant_features[, sub_level]] <-
        TRUE
      
      #Add stars to legends if sig_lab is defined as TRUE
      significant_features$stars <- ""
      
      if (sig_lab == T) {
        significant_features$stars <- symnum(
          significant_features$`padj`,
          symbols   = c("***", "**", "*", ""),
          cutpoints = c(0,  .001, .01, .05, 1),
          corr      = FALSE
        )
        
        #Add stars to the name of the taxa
        star_vec <-
          significant_features$stars[match(df_long$plot_taxa , significant_features[, sub_level])]
        star_vec[is.na(star_vec)]  <- ""
        df_long$plot_taxa <- paste0(df_long$plot_taxa, " ", star_vec)
        
        
      }
      
      #put significant features in bold
      df_long$legend_label <-
        ifelse(
          df_long$differential_abundance,
          paste0("<b>", df_long$plot_taxa, "</b>"),
          as.character(df_long$plot_taxa)
        )
      df_long$plot_taxa <- df_long$legend_label
      df_long$plot_taxa <-
        factor(df_long$plot_taxa, levels = unique(df_long$plot_taxa))
      
    }
    
  }
  
  
  # Differential abundance analysis on main_level using DESeq2 : 2 groups only
  if (differential_analysis &&
      length(unique(meta[[exp_group]])) == 2) {
    main_glom <- tax_glom(ps_object, taxrank = main_level)
    
    #remove OTU with 0 count across the dataset
    
    if (sum(rowSums(otu_table(main_glom)) == 0) > 0) {
      main_glom <-
        prune_taxa(rowSums(otu_table(main_glom)) > 0, main_glom)
      
    }
    
    
    diagdds_main = phyloseq_to_deseq2(main_glom,  formula(paste("~", exp_group)))
    
    # Run DESeq2 analysis
    if (test == "Wald"){
      
      diag_main = DESeq(diagdds_main, test = test, fitType = fitType, sfType = sfType, 
                        betaPrior = betaPrior, quiet= quiet, 
                        minReplicatesForReplace = minReplicatesForReplace,
                        modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                        parallel = parallel)
      
    }else {
      
      diag_main = DESeq(diagdds_main, test = test, fitType = fitType, sfType = sfType, 
                        betaPrior = betaPrior, reduced = formula(reduced), quiet= quiet, 
                        minReplicatesForReplace = minReplicatesForReplace,
                        modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                        parallel = parallel)
      
    }
    
    # Get the differentially abundant features
    results_main <- results(diag_main)
    
    # Extract features with adjusted p-value below the threshold
    significant_features_main <- subset(results_main, padj < fdr_threshold)
    
    #Check if significant features were detected
    if (nrow(significant_features_main) == 0) {
      message("no significant feature detected at ", main_level, " level.")
    } else {
      #We have significantly different taxa, we link them to their sub_level
      
      significant_features_main <-
        as.data.frame(cbind(significant_features_main,
                            tax_table(main_glom)[rownames(tax_table(main_glom)) %in% rownames(significant_features_main)]))
      
      
      # Add information to df_long to mark differentially abundant features
      df_long$differential_abundance_main <- FALSE
      df_long$differential_abundance_main[df_long[,main_level] %in% significant_features_main[, main_level]] <-
        TRUE
      
      #store the initial main_level names
      df_long$legend_label_main <- df_long[,main_level]
      
      #Add stars to legends if sig_lab is defined as TRUE
      significant_features_main$stars <- ""
      
      if (sig_lab == T) {
        significant_features_main$stars <- symnum(
          significant_features_main$`padj`,
          symbols   = c("***", "**", "*", ""),
          cutpoints = c(0,  .001, .01, .05, 1),
          corr      = FALSE
        )
        
        #Add stars to the name of the taxa
        star_vec_main <-
          significant_features_main$stars[match(df_long[,main_level], significant_features_main[, main_level])]
        star_vec_main[is.na(star_vec_main)]  <- ""
        df_long[,main_level] <- paste0(df_long[,main_level], " ", star_vec_main)
        
        #remove the stars column from the datafram 
        significant_features_main <-  significant_features_main[ ,-c(ncol(significant_features_main)) ]
        
      }
      
      #put significant features in bold
      
      
      df_long[,main_level] <-
        ifelse(
          df_long$differential_abundance_main,
          paste0("<b>", df_long[,main_level], "</b>"),
          as.character(df_long[,main_level])
        )
      #df_long[,main_level] <- df_long$legend_label_main
      df_long[,main_level] <-
        factor(df_long[,main_level], levels = unique(df_long[,main_level]))
      
    }
    
  }
  
  
  # Differential abundance analysis on sub_level using DESeq2 : multiple comparisons
  if (differential_analysis && mult_comp == T) {
    
    fam_glom <- tax_glom(ps_object, taxrank = sub_level)
    
    #remove OTU with 0 count across the dataset
    if (sum(rowSums(otu_table(fam_glom)) == 0) > 0) {
      fam_glom <-
        prune_taxa(rowSums(otu_table(fam_glom)) > 0, fam_glom)
      
    }
    
    if(twoFactor){
      
      diagdds = phyloseq_to_deseq2(fam_glom, formula(paste("~", 
                                                           fac1, "+", fac2, "+", fac1, ":", fac2))) # Full formula including interaction between factors
      
      colData(diagdds)[[fac1]] <- relevel(colData(diagdds)[[fac1]], ref=refFac1) #Setting refFac1 as the baseline for fac1
      
      colData(diagdds)[[fac2]] <- relevel(colData(diagdds)[[fac2]], ref=refFac2) #Setting refFac2 as the baseline for fac2
      
      
    }else{
      diagdds = phyloseq_to_deseq2(fam_glom,  formula(paste("~", exp_group)))
    }
    
    # Run DESeq2 analysis
    if (test == "Wald"){
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    } else {
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, reduced = formula(reduced), quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    }
    
    comparisons <- combn(levels(meta[[exp_group]]), 2, simplify = FALSE)
    
    
    if (is.null(selected_comparisons) == F)  {
      
      # Filter only selected comparisons
      comparisons <- Filter(function(cmp) {
        any(sapply(selected_comparisons, function(sel) all(sel == cmp)))
      }, comparisons)
    }
    
    # Dummy fix, TODO: fix issue, where comparisons do not inherit selected_comparisons order which makes us misinterpret results
    comparisons[[1]] = selected_comparisons[[1]]
    comparisons[[2]] = selected_comparisons[[2]]
    comparisons[[3]] = selected_comparisons[[3]]
    comparisons[[4]] = selected_comparisons[[4]]
    message("performed comparisons will follow the order")
    print(comparisons)

    # Generate names for the results based on comparisons
    comparison_names <- sapply(comparisons, function(cmp) paste(cmp, collapse = "_vs_"))
    
    # Using correctly formatted contrasts based on results names
    results_list <- list()
    if(twoFactor){
      
      #Partition results for specific pairwise comparaisons
      res_subset1 <- results(diag, contrast = list(c(resultsNames(diag)[3]))) #wt putrescine vs vehicle
      res_subset2 <- results(diag, contrast = list(c(resultsNames(diag)[3], resultsNames(diag)[4]))) #il22 ko putrescine vs vehicle
      res_subset3 <- results(diag, contrast = list(c(resultsNames(diag)[2]))) #vehicle wt vs il22 ko
      res_subset4 <- results(diag, contrast = list(c(resultsNames(diag)[2], resultsNames(diag)[4]))) #putrescine wt vs il22 ko
      results_list <- append(results_list, list(res_subset1, res_subset2, res_subset3, res_subset4)) #Append elements
      
    }else{
      
      results_list <- lapply(comparisons, function(cmp) {
        contrast_vec <- c(exp_group, cmp[1], cmp[2])  # cmp[1] is the numerator, cmp[2] is the denominator
        tryCatch({
          res <- results(diag, contrast = contrast_vec)
          return(res)
        }, error = function(e) {
          message("Failed to compute results for comparison: ", paste(cmp, collapse = " vs "), "\nError: ", e$message)
          return(NULL)  
        })
      })
      
    }
    
    # Set names based on comparisons
    results_list <- setNames(results_list, comparison_names)
    
    # Initialize an empty list to store significant features for each comparison
    significant_features_sub <- list()
    
    # Iterate over the results_list to process each comparison
    for (i in seq_along(results_list)) {
      res <- results_list[[i]]
      comparison_name <- names(results_list)[i]
      
      # Apply the significance threshold
      significant_features <- subset(res, padj < fdr_threshold)
      
      # Check if significant features were detected
      if (nrow(significant_features) == 0) {
        message("No significant feature detected for comparison ", comparison_name ," at the ", sub_level, " level.")
      } else {
        # Link significant features with their taxonomy information
        significant_features <- cbind(
          significant_features,
          tax_table(fam_glom)[rownames(tax_table(fam_glom)) %in% rownames(significant_features), , drop = FALSE]
        )
        
        # Convert to data frame and add to the list
        significant_features <- as.data.frame(significant_features)
        significant_features_sub[[comparison_name]] <- significant_features
      }
      
    }
    
  }
  
  # Differential abundance analysis on main_level using DESeq2 : multiple comparisons
  if (differential_analysis && mult_comp == T) {
    
    main_glom <- tax_glom(ps_object, taxrank = main_level)
    
    #remove OTU with 0 count across the dataset
    
    if (sum(rowSums(otu_table(main_glom)) == 0) > 0) {
      fam_glom <-
        prune_taxa(rowSums(otu_table(main_glom)) > 0, main_glom)
      
    }
    
    if(twoFactor){
      
      diagdds = phyloseq_to_deseq2(main_glom, formula(paste("~", 
                                                           fac1, "+", fac2, "+", fac1, ":", fac2))) # Full formula including interaction between factors
      
      colData(diagdds)[[fac1]] <- relevel(colData(diagdds)[[fac1]], ref=refFac1) #Setting refFac1 as the baseline for fac1
      
      colData(diagdds)[[fac2]] <- relevel(colData(diagdds)[[fac2]], ref=refFac2) #Setting refFac2 as the baseline for fac2
      
      
    }else{
      
      diagdds = phyloseq_to_deseq2(main_glom, formula(paste("~", 
                                                           exp_group)))
    }
    
    # Run DESeq2 analysis
    if (test == "Wald"){
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    } else {
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, reduced = formula(reduced), quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    }
    
    comparisons <- combn(levels(meta[[exp_group]]), 2, simplify = FALSE)
    
    if (is.null(selected_comparisons) == F)  {
      
      # Filter only selected comparisons
      comparisons <- Filter(function(cmp) {
        any(sapply(selected_comparisons, function(sel) all(sel == cmp)))
      }, comparisons)
    }
    
    # Dummy fix, TODO: fix issue, where comparisons do not inherit selected_comparisons order which makes us misinterpret results
    comparisons[[1]] = selected_comparisons[[1]]
    comparisons[[2]] = selected_comparisons[[2]]
    comparisons[[3]] = selected_comparisons[[3]]
    comparisons[[4]] = selected_comparisons[[4]]
    message("performed comparisons will follow the order")
    print(comparisons)
    
    # Generate names for the results based on comparisons
    comparison_names <- sapply(comparisons, function(cmp) paste(cmp, collapse = "_vs_"))
    
    # Using correctly formatted contrasts based on results names
    results_list <- list()
    if(twoFactor){
      
      #Partition results for specific pairwise comparaisons
      res_subset1 <- results(diag, contrast = list(c(resultsNames(diag)[3]))) #wt putrescine vs vehicle
      res_subset2 <- results(diag, contrast = list(c(resultsNames(diag)[3], resultsNames(diag)[4]))) #il22 ko putrescine vs vehicle
      res_subset3 <- results(diag, contrast = list(c(resultsNames(diag)[2]))) #vehicle wt vs il22 ko
      res_subset4 <- results(diag, contrast = list(c(resultsNames(diag)[2], resultsNames(diag)[4]))) #putrescine wt vs il22 ko
      results_list <- append(results_list, list(res_subset1, res_subset2, res_subset3, res_subset4)) #Append elements
      
    }else{
      
      results_list <- lapply(comparisons, function(cmp) {
        contrast_vec <- c(exp_group, cmp[1], cmp[2])  # cmp[1] is the numerator, cmp[2] is the denominator
        tryCatch({
          res <- results(diag, contrast = contrast_vec)
          return(res)
        }, error = function(e) {
          message("Failed to compute results for comparison: ", paste(cmp, collapse = " vs "), "\nError: ", e$message)
          return(NULL)  
        })
      })
      
    }
    
    # Set names based on comparisons
    results_list <- setNames(results_list, comparison_names)
    
    # Initialize an empty list to store significant features for each comparison
    significant_features_main <- list()
    
    # Iterate over the results_list to process each comparison
    for (i in seq_along(results_list)) {
      res <- results_list[[i]]
      comparison_name <- names(results_list)[i]
      
      # Apply the significance threshold
      significant_features <- subset(res, padj < fdr_threshold)
      
      # Check if significant features were detected
      if (nrow(significant_features) == 0) {
        message("No significant feature detected for comparison ", comparison_name," at the ", main_level, " level.")
      } else {
        
        # print(dim(significant_features))
        print(significant_features)
        print(tax_table(main_glom)[rownames(tax_table(main_glom)) %in% rownames(significant_features),, drop = FALSE])
        # print(dim(tax_table(main_glom)[rownames(tax_table(main_glom)) %in% rownames(significant_features), , drop = FALSE]))
        
        
        # Link significant features with their taxonomy information
        
        significant_features <- cbind(
          significant_features,
          tax_table(main_glom)[rownames(tax_table(main_glom)) %in% rownames(significant_features), , drop = FALSE]
        )
        
        
        # Convert to data frame and add to the list
        significant_features <- as.data.frame(significant_features)
        significant_features_main[[comparison_name]] <- significant_features
        
      }
    }
    
    
    
    
    
  }

  
  
  #Prepare the color vector
  MyColors <- df_long$MyColors
  names(MyColors) <- df_long$plot_taxa
  
  MyColors2 <- unique(df_long$MyColors)
  names(MyColors2) <- unique(df_long$plot_taxa)
  
  #add the black color to "Other_main_level"
  main_level_col[length(main_level_col)+1] <- '#000000'
  
  #Order the colors
  df_long[,main_level] <- factor(df_long[,main_level], levels = unique(df_long[,main_level]))
  
  vec1 <- unique(df_long[,main_level])
  vec2 <- c(pull(topx[,main_level]), paste0("Others"))
  core_text_vec1 <- gsub("(<[^>]*>|\\*|\\s+$)", "", vec1)
  core_text_vec1 <- trimws(core_text_vec1)
  order_index <- match( core_text_vec1, vec2)
  main_level_col <- main_level_col[order_index]
  names(main_level_col) <- as.character(vec1)
  
  #plot
  
  if(mean_group == F) {  
    
    
    p <-
      ggplot(df_long, aes(
        x = !!as.name(sample_name),
        y = value,
        fill = plot_taxa
      )) +
      geom_bar(stat = "identity", width = 0.85) +
      ylab("Relative abundance (%)\n") +
      guides(fill = guide_legend(reverse = FALSE, title = sub_level, order = 2)) +
      theme(
        line = element_line(colour = "black", linewidth = .5),
        text = element_text(size = 9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = .5)
      ) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = "right",
        legend.box = "vertical",
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = x_axis_size
        ),
        text = element_text(size = text_size),
        legend.text = element_markdown(size = legend_size)
      ) +
      geom_bar(aes(alpha = df_long[, main_level]), stat = "identity", show.legend = ifelse(showOnlySubLegend, FALSE, TRUE)) +
      scale_alpha_manual(
        values = rep(1, length(unique(df_long[, main_level]))),
        guide = guide_legend(order = 1,override.aes = list(fill = main_level_col)),
        name = main_level
      ) +
      scale_fill_manual('plot_taxa', values = MyColors2)
    p <- p + facet_wrap( ~ df_long[[exp_group]] , scales  = "free_x", nrow = n_row, ncol = n_col)
    
    p
    
  } else {
    
    
    df_long <- df_long %>% 
      group_by(!!as.name(exp_group), plot_taxa, !!as.name(main_level)) %>%
      reframe(
        n = n(),
        sum = sum(as.double(value)))
    df_long <- data.frame(df_long)
    
    i<-1
    for (i in 1: length(unique(df_long[,exp_group]))) {
      
      df_long[,"sum"][df_long[,exp_group] == unique(df_long[,exp_group] )[i]] <- df_long[,"sum"][df_long[,exp_group] == unique(df_long[,exp_group] )[i]]/count(meta[,exp_group] == unique(meta[,exp_group])[i])
      
    }
    
    colnames(df_long)[colnames(df_long) == "sum"] <- "value"  
    
    
    
    #plot
    p <-
      ggplot(df_long, aes(
        x = !!as.name(exp_group),
        y = value,
        fill = plot_taxa
      )) +
      geom_bar(stat = "identity", width = 0.85) +
      ylab("Relative abundance (%)\n") +
      guides(fill = guide_legend(reverse = FALSE, title = sub_level, order = 2)) +
      theme(
        line = element_line(colour = "black", linewidth = .5),
        text = element_text(size = 9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = .5)
      ) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = "right",
        legend.box = "vertical",
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = x_axis_size
        ),
        text = element_text(size = text_size),
        legend.text = element_markdown(size = legend_size)
      ) +
      geom_bar(aes(alpha = df_long[, main_level]), stat = "identity", show.legend = ifelse(showOnlySubLegend, FALSE, TRUE)) +
      scale_alpha_manual(
        values = rep(1, length(unique(df_long[, main_level]))),
        guide = guide_legend(order = 1,override.aes = list(fill = main_level_col)),
        name = main_level
      ) +
      scale_fill_manual('plot_taxa', values = MyColors2)
    p <- p + facet_wrap( ~ df_long[[exp_group]] , scales  = "free_x", nrow = n_row, ncol = n_col)
    
    p
    
    
    
  }
  
  if (differential_analysis == T ) {
    
    return(
      list(
        significant_table_main = significant_features_main,
        significant_table_sub = significant_features_sub,
        plot = p,
        main_names=unique(df_long[[main_level]]),
        sub_names=unique(df_long$plot_taxa)
      )
    )
    
  } else {
    return(list(plot = p))
  }
  
}  


# For timepoints design, does not compare between timepoints but compares groups
# of interest within timepoints 
plot_microbiota_timepoints <- function(ps_object = ps,
                            exp_group = 'group',
                            timePoints = TRUE,
                            time_variable = 'week',
                            combined_group = 'gg_group',
                            subset_group = NULL,
                            sample_name = 'SampleID',
                            main_level = 'Phylum',
                            sub_level = 'Family',
                            threshold = 1,
                            n_phy = 4,
                            mean_group = FALSE,
                            hues = c("Oranges", "Greens", "Blues", "Purples"),
                            color_bias = 2,
                            n_row = 1,
                            n_col = NULL,
                            text_size = 9,
                            legend_size = 7,
                            x_axis_size = 8,
                            differential_analysis = FALSE,
                            mult_comp = FALSE,
                            selected_comparisons = NULL,
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
                            parallel = FALSE,
                            showOnlySubLegend = FALSE
) {
  
  
  
  # Validate inputs
  if (!("phyloseq" %in% class(ps_object)))
    stop("ps_object must be a phyloseq-class object.")
  if (!main_level %in% colnames(tax_table(ps_object)))
    stop("main_level column not found in the tax_table.")
  if (!sub_level %in% colnames(tax_table(ps_object)))
    stop("sub_level column not found in the tax_table.")
  if (!exp_group %in% names(sample_data(ps_object)))
    stop("exp_group column not found in the sample_data.")
  
  
  #Assure taxa are rows
  if (!taxa_are_rows(ps_object)) {
    ps_object <- t(ps_object)
  }
  
  #Subset the samples belonging to the groups selected
  if (is.null(subset_group) == F) {
    keep_samples = as.character(get_variable(ps_object, exp_group)) %in% subset_group
    ps_object = prune_samples(keep_samples, ps_object)
    
  }
  
  #transform to relative abundance
  ps_prop <-
    transform_sample_counts(ps_object, function(OTU)
      ((OTU / sum(OTU)) * 100))
  
  
  #store TAX, OTU and metadata tables
  otu <-  as.data.frame(otu_table(ps_prop))
  tax <-  as.data.frame(tax_table(ps_prop))
  meta <- data.frame(sample_data(ps_prop))
  
  #clean up OTU absent and taxa not present
  
  otu <- otu %>% filter_all(any_vars(. != 0))
  
  tax <- subset(tax, rownames(tax) %in% rownames(otu))
  
  #Agglomerate at the sub level defined
  ps_f <- tax_glom(ps_prop, taxrank = sub_level , NArm = FALSE)
  #Store OTU and TAX
  otu_f <-  as.data.frame(otu_table(ps_f))
  tax_f <-  as.data.frame(tax_table(ps_f))
  
  #Bind TAX and OTU
  
  otu_tax_f <- cbind(tax_f, otu_f)
  
  
  #create a vector, if True, the ASV mean among the samples is higher than the threshold
  
  #keep the position of the column storing the abundance of the taxa in the otu_tax_f object
  position <- ncol(tax) + 1:(ncol(otu_tax_f) - ncol(tax))
  
  message(paste0('\n', length(position), ' samples are analyzed \n'))
  
  #Create 2 vectors in the dataframe :
  #1 storing the abundance mean of the taxa
  #2 T or F this taxa is above the define threshold of filtering
  
  otu_tax_f <- otu_tax_f %>%
    rowwise() %>%
    mutate(Mean = mean(c_across(all_of(position)))) %>%
    mutate(high_abundance = case_when(Mean > threshold  ~ TRUE,
                                      TRUE  ~ FALSE))
  
  #Store top x of the main_level phylogeny
  
  topx <- otu_tax_f %>%
    dplyr::group_by(!!as.name(main_level)) %>%
    dplyr::summarise(sum_top = sum(c_across(all_of(position)))) %>%
    dplyr::arrange(desc(sum_top)) %>%
    dplyr::slice_head(n = n_phy)
  
  
  #Create a vector, TRUE when the phylum belongs to the top X, the ones we want to plot
  otu_tax_f <- mutate(otu_tax_f,
                      selected_top = !!as.name(main_level) %in% pull(topx[, main_level]))
  
  #replace NA by "Unknown" 
  otu_tax_f[is.na(otu_tax_f)] <- "Unknown"
  
  #Add a column 'plot_taxa' to create the legend of the graphs :
  
  otu_tax_f <- otu_tax_f %>%
    mutate(
      plot_taxa = case_when(
        high_abundance == TRUE  &
          selected_top == TRUE &
          !!as.name(sub_level) == "Unknown"  ~
          paste0("Unknown ",!!as.name(main_level)),
        high_abundance == TRUE  &
          selected_top == TRUE  &
          !!as.name(sub_level) != "Unknown"  ~ !!as.name(sub_level),
        high_abundance == FALSE &
          selected_top == TRUE  ~ paste0("Others ", !!as.name(main_level)),
        selected_top   == FALSE ~ paste0("Others ")
      )
    )
  
  
  if (nrow(topx) != length(hues)) {
    message(
      'the number of colors choosen (',
      length(hues),
      ') is different from the defined number of features to plot (',
      nrow(topx),
      ')'
    )
  }
  
  
  #initialize 
  df <-
    as_tibble(matrix(nrow = 0, ncol = length(colnames(otu_tax_f))),
              .name_repair = ~ colnames(otu_tax_f))
  main_level_col <- c()
  
  #loop through selected main_level to add color
  i <- 1
  
  for (i in 1:nrow(topx)) {
    top_loop <- pull(topx[, main_level])[i]
    # print(top_loop)
    
    # Subset the df with unique features and in function of the Means abundance
    
    temp <- otu_tax_f %>%
      dplyr::filter(!!as.name(main_level) == top_loop) %>%
      dplyr::arrange(desc(Mean))
    
    message(
      'Among the ',
      main_level,
      ' ' ,
      top_loop,
      ' ' ,
      length(unique(temp$plot_taxa)) ,
      ' features at ',
      sub_level,
      ' level will be plotted'
    )
    
    #define the color palette
    getPalette = colorRampPalette(brewer.pal(length(unique(temp$plot_taxa)), hues[i]), bias = color_bias)
    col <-
      getPalette(length(unique(temp$plot_taxa)) + 1)
    
    #store the darkestcolor of the palette to plot main_level legend graph
    if (length(col) >= 3) {
      main_level_col[i] <-   col[length(col)-1]
    } else {
      main_level_col[i] <-   col[length(col)]
    }  
    
    #loop among the sub_level feature and add color to it
    u <- 1
    for (u in 1:length(unique(temp$plot_taxa))) {
      print(unique(temp$plot_taxa)[u])
      
      #Add the color
      temp[which(temp$plot_taxa ==  unique(temp$plot_taxa)[u]), 'MyColors'] <-
        col[u + 1]
      
    }
    
    #create the dataframe storing the colors information
    df <- rbind(df, temp)
    
  }
  
  
  #Select features that don't belong to high abundance levels and assigned the black color to them
  unselescted <- otu_tax_f %>%
    filter(!(!!as.name(main_level) %in% pull(topx[, main_level])))
  
  unselescted$MyColors <- '#000000'
  
  #Merge the 2 dataframes
  df <- rbind(df, unselescted)
  
  #Order df to define the order that will be plotted, first the selected main_level, then the main_level in alphabetical order, then the mean abundance
  
  df <- df %>% ungroup %>%
    dplyr::arrange(desc(selected_top), !!as.name(main_level) , desc(Mean))
  
  
  #check if no features have been lost
  nrow(df) == nrow(otu_tax_f)
  
  #factor to keep order in the plot
  df$plot_taxa <-
    factor(df$plot_taxa, levels = unique(df$plot_taxa))
  
  
  #Transform df into long format by sample
  df_long <- melt(
    df,
    id = c("plot_taxa", "MyColors", main_level),
    measure.vars = meta[, sample_name],
    variable.name = sample_name
  )
  
  
  #Add exp_group to df_long
  df_long <- left_join(df_long, meta, by = sample_name)
  
  
  #Replace low abundance features at level X by "Others "
  
  df_long[,main_level] <- ifelse(df_long[,main_level] %in% pull(topx[, main_level]), 
                                 df_long[,main_level], 
                                 "Others ")
  
  
  if (differential_analysis &&
      length(unique(meta[[exp_group]])) != 2) {
    print(
      "The number of experimental group to test is not equal to 2, no statistical significance will appear in the legend"
    )
    mult_comp <- T
    
  }
  
  
  
  # # Differential abundance analysis on sub_level using DESeq2 : 2 groups only
  # if (differential_analysis &&
  #     length(unique(meta[[exp_group]])) == 2) {
  #   fam_glom <- tax_glom(ps_object, taxrank = sub_level)
  #   
  #   #remove OTU with 0 count across the dataset
  #   
  #   if (sum(rowSums(otu_table(fam_glom)) == 0) > 0) {
  #     fam_glom <-
  #       prune_taxa(rowSums(otu_table(fam_glom)) > 0, fam_glom)
  #     
  #   }
  #   
  #   diagdds = phyloseq_to_deseq2(fam_glom,  formula(paste("~", exp_group)))
  #   
  #   # Run DESeq2 analysis
  #   if (test == "Wald"){
  #     
  #     
  #     
  #     diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
  #                  betaPrior = betaPrior, quiet= quiet, 
  #                  minReplicatesForReplace = minReplicatesForReplace,
  #                  modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
  #                  parallel = parallel)
  #     
  #   } else {
  #     
  #     diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
  #                  betaPrior = betaPrior, reduced = formula(reduced), quiet= quiet, 
  #                  minReplicatesForReplace = minReplicatesForReplace,
  #                  modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
  #                  parallel = parallel)
  #     
  #   }
  #   
  #   # Get the differentially abundant features
  #   results <- results(diag)
  #   
  #   # Extract features with adjusted p-value below the threshold
  #   significant_features <- subset(results, padj < fdr_threshold)
  #   
  #   #Check if significant features were detected
  #   if (nrow(significant_features) == 0) {
  #     message("no significant feature detected at ", sub_level, " level.")
  #   } else {
  #     #We have significantly different taxa, we link them to their sub_level
  #     
  #     
  #     significant_features <-
  #       as.data.frame(cbind(significant_features,
  #                           tax_table(fam_glom)[rownames(tax_table(fam_glom)) %in% rownames(significant_features)]))
  #     
  #     #store the sig table to return it 
  #     significant_features_sub <- significant_features
  #     
  #     # Add information to df_long to mark differentially abundant features
  #     df_long$differential_abundance <- FALSE
  #     df_long$differential_abundance[df_long$plot_taxa %in% significant_features[, sub_level]] <-
  #       TRUE
  #     
  #     #Add stars to legends if sig_lab is defined as TRUE
  #     significant_features$stars <- ""
  #     
  #     if (sig_lab == T) {
  #       significant_features$stars <- symnum(
  #         significant_features$`padj`,
  #         symbols   = c("***", "**", "*", ""),
  #         cutpoints = c(0,  .001, .01, .05, 1),
  #         corr      = FALSE
  #       )
  #       
  #       #Add stars to the name of the taxa
  #       star_vec <-
  #         significant_features$stars[match(df_long$plot_taxa , significant_features[, sub_level])]
  #       star_vec[is.na(star_vec)]  <- ""
  #       df_long$plot_taxa <- paste0(df_long$plot_taxa, " ", star_vec)
  #       
  #       
  #     }
  #     
  #     #put significant features in bold
  #     df_long$legend_label <-
  #       ifelse(
  #         df_long$differential_abundance,
  #         paste0("<b>", df_long$plot_taxa, "</b>"),
  #         as.character(df_long$plot_taxa)
  #       )
  #     df_long$plot_taxa <- df_long$legend_label
  #     df_long$plot_taxa <-
  #       factor(df_long$plot_taxa, levels = unique(df_long$plot_taxa))
  #     
  #   }
  #   
  # }
  # 
  # 
  # # Differential abundance analysis on main_level using DESeq2 : 2 groups only
  # if (differential_analysis &&
  #     length(unique(meta[[exp_group]])) == 2) {
  #   main_glom <- tax_glom(ps_object, taxrank = main_level)
  #   
  #   #remove OTU with 0 count across the dataset
  #   
  #   if (sum(rowSums(otu_table(main_glom)) == 0) > 0) {
  #     main_glom <-
  #       prune_taxa(rowSums(otu_table(main_glom)) > 0, main_glom)
  #     
  #   }
  #   
  #   
  #   diagdds_main = phyloseq_to_deseq2(main_glom,  formula(paste("~", exp_group)))
  #   
  #   # Run DESeq2 analysis
  #   if (test == "Wald"){
  #     
  #     diag_main = DESeq(diagdds_main, test = test, fitType = fitType, sfType = sfType, 
  #                       betaPrior = betaPrior, quiet= quiet, 
  #                       minReplicatesForReplace = minReplicatesForReplace,
  #                       modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
  #                       parallel = parallel)
  #     
  #   }else {
  #     
  #     diag_main = DESeq(diagdds_main, test = test, fitType = fitType, sfType = sfType, 
  #                       betaPrior = betaPrior, reduced = formula(reduced), quiet= quiet, 
  #                       minReplicatesForReplace = minReplicatesForReplace,
  #                       modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
  #                       parallel = parallel)
  #     
  #   }
  #   
  #   # Get the differentially abundant features
  #   results_main <- results(diag_main)
  #   
  #   # Extract features with adjusted p-value below the threshold
  #   significant_features_main <- subset(results_main, padj < fdr_threshold)
  #   
  #   #Check if significant features were detected
  #   if (nrow(significant_features_main) == 0) {
  #     message("no significant feature detected at ", main_level, " level.")
  #   } else {
  #     #We have significantly different taxa, we link them to their sub_level
  #     
  #     significant_features_main <-
  #       as.data.frame(cbind(significant_features_main,
  #                           tax_table(main_glom)[rownames(tax_table(main_glom)) %in% rownames(significant_features_main)]))
  #     
  #     
  #     # Add information to df_long to mark differentially abundant features
  #     df_long$differential_abundance_main <- FALSE
  #     df_long$differential_abundance_main[df_long[,main_level] %in% significant_features_main[, main_level]] <-
  #       TRUE
  #     
  #     #store the initial main_level names
  #     df_long$legend_label_main <- df_long[,main_level]
  #     
  #     #Add stars to legends if sig_lab is defined as TRUE
  #     significant_features_main$stars <- ""
  #     
  #     if (sig_lab == T) {
  #       significant_features_main$stars <- symnum(
  #         significant_features_main$`padj`,
  #         symbols   = c("***", "**", "*", ""),
  #         cutpoints = c(0,  .001, .01, .05, 1),
  #         corr      = FALSE
  #       )
  #       
  #       #Add stars to the name of the taxa
  #       star_vec_main <-
  #         significant_features_main$stars[match(df_long[,main_level], significant_features_main[, main_level])]
  #       star_vec_main[is.na(star_vec_main)]  <- ""
  #       df_long[,main_level] <- paste0(df_long[,main_level], " ", star_vec_main)
  #       
  #       #remove the stars column from the datafram 
  #       significant_features_main <-  significant_features_main[ ,-c(ncol(significant_features_main)) ]
  #       
  #     }
  #     
  #     #put significant features in bold
  #     
  #     
  #     df_long[,main_level] <-
  #       ifelse(
  #         df_long$differential_abundance_main,
  #         paste0("<b>", df_long[,main_level], "</b>"),
  #         as.character(df_long[,main_level])
  #       )
  #     #df_long[,main_level] <- df_long$legend_label_main
  #     df_long[,main_level] <-
  #       factor(df_long[,main_level], levels = unique(df_long[,main_level]))
  #     
  #   }
  #   
  # }
  
  
  # Differential abundance analysis on sub_level using DESeq2 : multiple comparisons
  if (differential_analysis && mult_comp == T) {
    
    fam_glom <- tax_glom(ps_object, taxrank = sub_level)
    
    #remove OTU with 0 count across the dataset
    
    if (sum(rowSums(otu_table(fam_glom)) == 0) > 0) {
      fam_glom <-
        prune_taxa(rowSums(otu_table(fam_glom)) > 0, fam_glom)
      
    }
    
    diagdds = phyloseq_to_deseq2(fam_glom,  formula(paste("~", exp_group)))
    
    # Run DESeq2 analysis
    if (test == "Wald"){
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    } else {
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, reduced = formula(reduced), quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    }
    
    comparisons <- combn(levels(meta[[exp_group]]), 2, simplify = FALSE)
    
    if (is.null(selected_comparisons) == F)  {
      
      # Filter only selected comparisons
      comparisons <- Filter(function(cmp) {
        any(sapply(selected_comparisons, function(sel) all(sel == cmp)))
      }, comparisons)
    }
    
    
    # Generate names for the results based on comparisons
    comparison_names <- sapply(comparisons, function(cmp) paste(cmp, collapse = "_vs_"))
    
    # Using correctly formatted contrasts based on results names
    results_list <- list()
    results_list <- lapply(comparisons, function(cmp) {
      contrast_vec <- c(exp_group, cmp[1], cmp[2])  # cmp[1] is the numerator, cmp[2] is the denominator
      tryCatch({
        res <- results(diag, contrast = contrast_vec)
        return(res)
      }, error = function(e) {
        message("Failed to compute results for comparison: ", paste(cmp, collapse = " vs "), "\nError: ", e$message)
        return(NULL)  
      })
    })
    
    # Set names based on comparisons
    results_list <- setNames(results_list, comparison_names)
    
    
    # Initialize an empty list to store significant features for each comparison
    significant_features_sub <- list()
    
    # Iterate over the results_list to process each comparison
    for (i in seq_along(results_list)) {
      res <- results_list[[i]]
      comparison_name <- names(results_list)[i]
      
      # Apply the significance threshold
      significant_features <- subset(res, padj < fdr_threshold)
      
      # Check if significant features were detected
      if (nrow(significant_features) == 0) {
        message("No significant feature detected for comparison ", comparison_name ," at the ", sub_level, " level.")
      } else {
        # Link significant features with their taxonomy information
        significant_features <- cbind(
          significant_features,
          tax_table(fam_glom)[rownames(tax_table(fam_glom)) %in% rownames(significant_features), , drop = FALSE]
        )
        
        # Convert to data frame and add to the list
        significant_features <- as.data.frame(significant_features)
        significant_features_sub[[comparison_name]] <- significant_features
      }
    }
    
    
    
    
    
  }
  
  # Differential abundance analysis on main_level using DESeq2 : multiple comparisons
  if (differential_analysis && mult_comp == T) {
    
    main_glom <- tax_glom(ps_object, taxrank = main_level)
    
    #remove OTU with 0 count across the dataset
    
    if (sum(rowSums(otu_table(main_glom)) == 0) > 0) {
      fam_glom <-
        prune_taxa(rowSums(otu_table(main_glom)) > 0, main_glom)
      
    }
    
    diagdds = phyloseq_to_deseq2(main_glom,  formula(paste("~", exp_group)))
    
    # Run DESeq2 analysis
    if (test == "Wald"){
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    } else {
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, reduced = formula(reduced), quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    }
    
    
    comparisons <- combn(levels(meta[[exp_group]]), 2, simplify = FALSE)
    
    if (is.null(selected_comparisons) == F)  {
      
      # Filter only selected comparisons
      comparisons <- Filter(function(cmp) {
        any(sapply(selected_comparisons, function(sel) all(sel == cmp)))
      }, comparisons)
    }
    
    # Generate names for the results based on comparisons
    comparison_names <- sapply(comparisons, function(cmp) paste(cmp, collapse = "_vs_"))
    
    # Using correctly formatted contrasts based on results names
    results_list <- list()
    results_list <- lapply(comparisons, function(cmp) {
      contrast_vec <- c(exp_group, cmp[1], cmp[2])  # cmp[1] is the numerator, cmp[2] is the denominator
      tryCatch({
        res <- results(diag, contrast = contrast_vec)
        return(res)
      }, error = function(e) {
        message("Failed to compute results for comparison: ", paste(cmp, collapse = " vs "), "\nError: ", e$message)
        return(NULL)  
      })
    })
    
    # Set names based on comparisons
    results_list <- setNames(results_list, comparison_names)
    
    
    # Initialize an empty list to store significant features for each comparison
    significant_features_main <- list()
    
    # Iterate over the results_list to process each comparison
    for (i in seq_along(results_list)) {
      res <- results_list[[i]]
      comparison_name <- names(results_list)[i]
      
      # Apply the significance threshold
      significant_features <- subset(res, padj < fdr_threshold)
      
      # Check if significant features were detected
      if (nrow(significant_features) == 0) {
        message("No significant feature detected for comparison ", comparison_name," at the ", main_level, " level.")
      } else {
        # Link significant features with their taxonomy information
        significant_features <- cbind(
          significant_features,
          tax_table(main_glom)[rownames(tax_table(main_glom)) %in% rownames(significant_features), , drop = FALSE]
        )
        
        # Convert to data frame and add to the list
        significant_features <- as.data.frame(significant_features)
        significant_features_main[[comparison_name]] <- significant_features
      }
    }
    
    
    
    
    
  }
  
  # Differential abundance analysis on sub_level using DESeq2 : multiple timepoints
  if (differential_analysis && timePoints == T) {
    
    fam_glom <- tax_glom(ps_object, taxrank = sub_level)
    
    #remove OTU with 0 count across the dataset
    
    if (sum(rowSums(otu_table(fam_glom)) == 0) > 0) {
      fam_glom <-
        prune_taxa(rowSums(otu_table(fam_glom)) > 0, fam_glom)
      
    }
    
    comparisons <- combn(levels(meta[[combined_group]]), 2, simplify = FALSE)
    
    if (is.null(selected_comparisons) == F)  {
      
      # Filter only selected comparisons
      comparisons <- Filter(function(cmp) {
        any(sapply(selected_comparisons, function(sel) all(sel == cmp)))
      }, comparisons)
    }
    
    # Generate names for the results based on comparisons
    comparison_names <- sapply(comparisons, function(cmp) paste(cmp, collapse = "_vs_"))
    
    # Initialize empty list for results
    results_list <- list()
    
    # Iterate through timepoints and create a separate deseq object for each
    for(i in seq_along(levels(sample_data(ps_object)[[time_variable]]))){
      
      #Creating phyloseq objects for each timepoint
      fam_subset <- prune_samples(sample_data(ps_object)[[time_variable]] == levels(sample_data(ps_object)[[time_variable]])[i], fam_glom)
      
      diagdds = phyloseq_to_deseq2(fam_subset,  formula(paste("~", exp_group)))
      
      # Run DESeq2 analysis
      if (test == "Wald"){
        
        diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                     betaPrior = betaPrior, quiet= quiet, 
                     minReplicatesForReplace = minReplicatesForReplace,
                     modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                     parallel = parallel)
        
      } else {
        
        diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                     betaPrior = betaPrior, reduced = formula(reduced), quiet= quiet, 
                     minReplicatesForReplace = minReplicatesForReplace,
                     modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                     parallel = parallel)
        
      }
      
      results_list[[comparison_names[i]]] <- results(diag, contrast = list(c(resultsNames(diag)[2])))
    }
    
    # # Using correctly formatted contrasts based on results names
    # 
    # results_list <- lapply(comparisons, function(cmp) {
    #   contrast_vec <- c(exp_group, cmp[1], cmp[2])  # cmp[1] is the numerator, cmp[2] is the denominator
    #   tryCatch({
    #     res <- results(diag, contrast = contrast_vec)
    #     return(res)
    #   }, error = function(e) {
    #     message("Failed to compute results for comparison: ", paste(cmp, collapse = " vs "), "\nError: ", e$message)
    #     return(NULL)  
    #   })
    # })
    # 
    # # Set names based on comparisons
    # results_list <- setNames(results_list, comparison_names)
    
    # Initialize an empty list to store significant features for each comparison
    significant_features_sub <- list()
    
    # Iterate over the results_list to process each comparison
    for (i in seq_along(results_list)) {
      res <- results_list[[i]]
      comparison_name <- names(results_list)[i]
      
      # Apply the significance threshold
      significant_features <- subset(res, padj < fdr_threshold)
      
      # Check if significant features were detected
      if (nrow(significant_features) == 0) {
        message("No significant feature detected for comparison ", comparison_name ," at the ", sub_level, " level.")
      } else {
        # Link significant features with their taxonomy information
        significant_features <- cbind(
          significant_features,
          tax_table(fam_subset)[rownames(tax_table(fam_subset)) %in% rownames(significant_features), , drop = FALSE]
        )
        
        # Convert to data frame and add to the list
        significant_features <- as.data.frame(significant_features)
        significant_features_sub[[comparison_name]] <- significant_features
      }
    }
    
      
  }
  
  # Differential abundance analysis on main_level using DESeq2 : multiple timepoints
  if (differential_analysis && timePoints == T) {
    
    main_glom <- tax_glom(ps_object, taxrank = main_level)
    
    #remove OTU with 0 count across the dataset
    
    if (sum(rowSums(otu_table(main_glom)) == 0) > 0) {
      main_glom <-
        prune_taxa(rowSums(otu_table(main_glom)) > 0, main_glom)
      
    }
    
    comparisons <- combn(levels(meta[[combined_group]]), 2, simplify = FALSE)
    
    if (is.null(selected_comparisons) == F)  {
      
      # Filter only selected comparisons
      comparisons <- Filter(function(cmp) {
        any(sapply(selected_comparisons, function(sel) all(sel == cmp)))
      }, comparisons)
    }
    
    # Generate names for the results based on comparisons
    comparison_names <- sapply(comparisons, function(cmp) paste(cmp, collapse = "_vs_"))
    
    # Initialize empty list for results
    results_list <- list()
    
    # Iterate through timepoints and create a separate deseq object for each
    for(i in seq_along(levels(sample_data(ps_object)[[time_variable]]))){
      
      #Creating phyloseq objects for each timepoint
      main_subset <- prune_samples(sample_data(ps_object)[[time_variable]] == levels(sample_data(ps_object)[[time_variable]])[i], main_glom)
      
      diagdds = phyloseq_to_deseq2(main_subset,  formula(paste("~", exp_group)))
      
      # Run DESeq2 analysis
      if (test == "Wald"){
        
        diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                     betaPrior = betaPrior, quiet= quiet, 
                     minReplicatesForReplace = minReplicatesForReplace,
                     modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                     parallel = parallel)
        
      } else {
        
        diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                     betaPrior = betaPrior, reduced = formula(reduced), quiet= quiet, 
                     minReplicatesForReplace = minReplicatesForReplace,
                     modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                     parallel = parallel)
        
      }
      
      results_list[[comparison_names[i]]] <- results(diag, contrast = list(c(resultsNames(diag)[2])))
    }
    
    # # Using correctly formatted contrasts based on results names
    # 
    # results_list <- lapply(comparisons, function(cmp) {
    #   contrast_vec <- c(exp_group, cmp[1], cmp[2])  # cmp[1] is the numerator, cmp[2] is the denominator
    #   tryCatch({
    #     res <- results(diag, contrast = contrast_vec)
    #     return(res)
    #   }, error = function(e) {
    #     message("Failed to compute results for comparison: ", paste(cmp, collapse = " vs "), "\nError: ", e$message)
    #     return(NULL)  
    #   })
    # })
    # 
    # # Set names based on comparisons
    # results_list <- setNames(results_list, comparison_names)
    
    
    # Initialize an empty list to store significant features for each comparison
    significant_features_main <- list()
    
    # Iterate over the results_list to process each comparison
    for (i in seq_along(results_list)) {
      res <- results_list[[i]]
      comparison_name <- names(results_list)[i]
      
      # Apply the significance threshold
      significant_features <- subset(res, padj < fdr_threshold)
      
      # Check if significant features were detected
      if (nrow(significant_features) == 0) {
        message("No significant feature detected for comparison ", comparison_name ," at the ", main_level, " level.")
      } else {
        # Link significant features with their taxonomy information
        significant_features <- cbind(
          significant_features,
          tax_table(main_subset)[rownames(tax_table(main_subset)) %in% rownames(significant_features), , drop = FALSE]
        )
        
        # Convert to data frame and add to the list
        significant_features <- as.data.frame(significant_features)
        significant_features_main[[comparison_name]] <- significant_features
      }
    }
    
    
  }
  
  
  #Prepare the color vector
  MyColors <- df_long$MyColors
  names(MyColors) <- df_long$plot_taxa
  
  MyColors2 <- unique(df_long$MyColors)
  names(MyColors2) <- unique(df_long$plot_taxa)
  
  #add the black color to "Other_main_level"
  main_level_col[length(main_level_col)+1] <- '#000000'
  
  #Order the colors
  df_long[,main_level] <- factor(df_long[,main_level], levels = unique(df_long[,main_level]))
  
  vec1 <- unique(df_long[,main_level])
  vec2 <- c(pull(topx[,main_level]), paste0("Others"))
  core_text_vec1 <- gsub("(<[^>]*>|\\*|\\s+$)", "", vec1)
  core_text_vec1 <- trimws(core_text_vec1)
  order_index <- match( core_text_vec1, vec2)
  main_level_col <- main_level_col[order_index]
  names(main_level_col) <- as.character(vec1)
  
  #plot
  
  if(mean_group == F) {  
    
    
    p <-
      ggplot(df_long, aes(
        x = !!as.name(sample_name),
        y = value,
        fill = plot_taxa
      )) +
      geom_bar(stat = "identity", width = 0.85) +
      ylab("Relative abundance (%)\n") +
      guides(fill = guide_legend(reverse = FALSE, title = sub_level, order = 2)) +
      theme(
        line = element_line(colour = "black", linewidth = .5),
        text = element_text(size = 9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = .5)
      ) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = "right",
        legend.box = "vertical",
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = x_axis_size
        ),
        text = element_text(size = text_size),
        legend.text = element_markdown(size = legend_size)
      ) +
      geom_bar(aes(alpha = df_long[, main_level]), stat = "identity", show.legend = ifelse(showOnlySubLegend, FALSE, TRUE)) +
      scale_alpha_manual(
        values = rep(1, length(unique(df_long[, main_level]))),
        guide = guide_legend(order = 1,override.aes = list(fill = main_level_col)),
        name = main_level
      ) +
      scale_fill_manual('plot_taxa', values = MyColors2) 
    p <- p + facet_wrap( ~ df_long[[combined_group]] , scales  = "free_x", nrow = n_row, ncol = n_col) # Replaced "exp_group" by "combined_group"
    
    p
    
  } else {
    
    
    df_long <- df_long %>% 
      group_by(!!as.name(exp_group), plot_taxa, !!as.name(main_level)) %>%
      reframe(
        n = n(),
        sum = sum(as.double(value)))
    df_long <- data.frame(df_long)
    
    i<-1
    for (i in 1: length(unique(df_long[,exp_group]))) {
      
      df_long[,"sum"][df_long[,exp_group] == unique(df_long[,exp_group] )[i]] <- df_long[,"sum"][df_long[,exp_group] == unique(df_long[,exp_group] )[i]]/count(meta[,exp_group] == unique(meta[,exp_group])[i])
      
    }
    
    colnames(df_long)[colnames(df_long) == "sum"] <- "value"  
    
    
    
    #plot
    p <-
      ggplot(df_long, aes(
        x = !!as.name(exp_group),
        y = value,
        fill = plot_taxa
      )) +
      geom_bar(stat = "identity", width = 0.85) +
      ylab("Relative abundance (%)\n") +
      guides(fill = guide_legend(reverse = FALSE, title = sub_level, order = 2)) +
      theme(
        line = element_line(colour = "black", linewidth = .5),
        text = element_text(size = 9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = .5)
      ) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = "right",
        legend.box = "vertical",
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = x_axis_size
        ),
        text = element_text(size = text_size),
        legend.text = element_markdown(size = legend_size)
      ) +
      geom_bar(aes(alpha = df_long[, main_level]), stat = "identity", show.legend = TRUE) +
      scale_alpha_manual(
        values = rep(1, length(unique(df_long[, main_level]))),
        guide = guide_legend(order = 1,override.aes = list(fill = main_level_col)),
        name = main_level
      ) +
      scale_fill_manual('plot_taxa', values = MyColors2)
    p <- p + facet_wrap( ~ df_long[[exp_group]] , scales  = "free_x", nrow = n_row, ncol = n_col)
    
    p
    
    
    
  }
  
  
  
  if (differential_analysis == T ) {
    
    return(
      list(
        significant_table_main = significant_features_main,
        significant_table_sub = significant_features_sub,
        plot = p,
        main_names = unique(df_long[[main_level]]),
        sub_names = unique(df_long$plot_taxa)
      )
    )
    
  } else {
    return(list(plot = p))
  }
  
}  


plot_microbiota_multiFac_timepoints <- function(
    ps_object = ps,
    exp_group = 'group',
    subset_group = NULL,
    timePoints = TRUE,
    time_variable = 'week',
    combined_group = 'gg_group',
    twoFactor = FALSE,
    fac1 = NULL,
    refFac1 = NULL,
    fac2 = NULL,
    refFac2 = NULL,
    sample_name = 'SampleID',
    main_level = 'Phylum',
    sub_level = 'Family',
    threshold = 1,
    n_phy = 4,
    mean_group = FALSE,
    hues = c("Oranges", "Greens", "Blues", "Purples"),
    color_bias = 2,
    n_row = 1,
    n_col = NULL,
    text_size = 9,
    legend_size = 7,
    x_axis_size = 8,
    differential_analysis = FALSE,
    mult_comp = FALSE,
    selected_comparisons = NULL,
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
    parallel = FALSE,
    showOnlySubLegend = FALSE
) {
  #### 1. Data Validation & Preparation ####
  if (!("phyloseq" %in% class(ps_object)))
    stop("ps_object must be a phyloseq-class object.")
  if (!main_level %in% colnames(tax_table(ps_object)))
    stop("main_level column not found in the tax_table.")
  if (!sub_level %in% colnames(tax_table(ps_object)))
    stop("sub_level column not found in the tax_table.")
  if (!exp_group %in% names(sample_data(ps_object)))
    stop("exp_group column not found in the sample_data.")
  
  if(timePoints){
    if(!combined_group %in% names(sample_data(ps_object)))
      stop("combined_group column not found in sample_data (needed for timepoints design).")
    if(!time_variable %in% names(sample_data(ps_object)))
      stop("time_variable column not found in sample_data (needed for timepoints design).")
  }
  
  # Ensure taxa are rows
  if (!taxa_are_rows(ps_object)) {
    ps_object <- t(ps_object)
  }
  
  # Subset samples if subset_group is provided
  if (!is.null(subset_group)) {
    keep_samples <- as.character(get_variable(ps_object, exp_group)) %in% subset_group
    ps_object <- prune_samples(keep_samples, ps_object)
  }
  
  # Transform counts to relative abundance (%)
  ps_prop <- transform_sample_counts(ps_object, function(OTU) ((OTU / sum(OTU)) * 100))
  
  # Extract OTU, taxonomy, and metadata tables
  otu <- as.data.frame(otu_table(ps_prop))
  tax <- as.data.frame(tax_table(ps_prop))
  meta <- data.frame(sample_data(ps_prop))
  
  # Remove OTUs with zero counts across all samples
  otu <- otu %>% filter_all(any_vars(. != 0))
  tax <- subset(tax, rownames(tax) %in% rownames(otu))
  
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
  
  if (nrow(topx) != length(hues)) {
    message('the number of colors chosen (', length(hues),
            ') is different from the defined number of features to plot (', nrow(topx), ')')
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
  df_long <- melt(df,
                  id = c("plot_taxa", "MyColors", main_level),
                  measure.vars = meta[, sample_name],
                  variable.name = sample_name)
  df_long <- left_join(df_long, meta, by = sample_name)
  df_long[, main_level] <- ifelse(df_long[, main_level] %in% pull(topx[, main_level]),
                                  df_long[, main_level],
                                  "Others ")
  
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
  if (mean_group == FALSE) {
    p <- ggplot(df_long, aes(x = !!as.name(sample_name), y = value, fill = plot_taxa)) +
      geom_bar(stat = "identity", width = 0.85) +
      ylab("Relative abundance (%)\n") +
      guides(fill = guide_legend(reverse = FALSE, title = sub_level, order = 2)) +
      theme(line = element_line(colour = "black", linewidth = 0.5),
            text = element_text(size = 9),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black", linewidth = 0.5)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            legend.position = "right",
            legend.box = "vertical",
            axis.text.x = element_text(angle = 45, hjust = 1, size = x_axis_size),
            text = element_text(size = text_size),
            legend.text = element_markdown(size = legend_size)) +
      geom_bar(aes(alpha = df_long[, main_level]), stat = "identity",
               show.legend = ifelse(showOnlySubLegend, FALSE, TRUE)) +
      scale_alpha_manual(values = rep(1, length(unique(df_long[, main_level]))),
                         guide = guide_legend(order = 1, override.aes = list(fill = main_level_col)),
                         name = main_level) +
      scale_fill_manual('plot_taxa', values = MyColors2)
    if(timePoints){
      p <- p + facet_wrap(~ df_long[[combined_group]], scales = "free_x", nrow = n_row, ncol = n_col)
    } else {
      p <- p + facet_wrap(~ df_long[[exp_group]], scales = "free_x", nrow = n_row, ncol = n_col)
    }
  } else {
    df_long <- df_long %>% group_by(!!as.name(exp_group), plot_taxa, !!as.name(main_level)) %>%
      reframe(n = n(), sum = sum(as.double(value)))
    df_long <- data.frame(df_long)
    for (i in 1:length(unique(df_long[, exp_group]))) {
      df_long[,"sum"][df_long[, exp_group] == unique(df_long[, exp_group])[i]] <-
        df_long[,"sum"][df_long[, exp_group] == unique(df_long[, exp_group])[i]] /
        count(meta[, exp_group] == unique(meta[, exp_group])[i])
    }
    colnames(df_long)[colnames(df_long) == "sum"] <- "value"
    p <- ggplot(df_long, aes(x = !!as.name(exp_group), y = value, fill = plot_taxa)) +
      geom_bar(stat = "identity", width = 0.85) +
      ylab("Relative abundance (%)\n") +
      guides(fill = guide_legend(reverse = FALSE, title = sub_level, order = 2)) +
      theme(line = element_line(colour = "black", linewidth = 0.5),
            text = element_text(size = 9),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = "black", linewidth = 0.5)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            legend.position = "right",
            legend.box = "vertical",
            axis.text.x = element_text(angle = 45, hjust = 1, size = x_axis_size),
            text = element_text(size = text_size),
            legend.text = element_markdown(size = legend_size)) +
      geom_bar(aes(alpha = df_long[, main_level]), stat = "identity", show.legend = TRUE) +
      scale_alpha_manual(values = rep(1, length(unique(df_long[, main_level]))),
                         guide = guide_legend(order = 1, override.aes = list(fill = main_level_col)),
                         name = main_level) +
      scale_fill_manual('plot_taxa', values = MyColors2)
    p <- p + facet_wrap(~ df_long[[exp_group]], scales = "free_x", nrow = n_row, ncol = n_col)
  }
  
  #### 6. Return ####
  if (differential_analysis == TRUE) {
    return(list(
      significant_table_main = significant_features_main,
      significant_table_sub = significant_features_sub,
      plot = p,
      main_names = unique(df_long[[main_level]]),
      sub_names = unique(df_long$plot_taxa)
    ))
  } else {
    return(list(plot = p))
  }
}




# Function to write and save stackbarExtended sig_table
writeStackbarExtendedSigTable <- function(main_table, includeSubTable = FALSE, sub_table = NULL, filepath){
  
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
    
    # Append to the final table
    table_to_write <- rbind(table_to_write, table)
    
  }
  
  if(includeSubTable){
    
    # Iterate over the list of tables for sub table
    for (i in seq_along(sub_table)){
      
      # Extract the table
      table <- sub_table[[i]]
      
      # Add a column with the name of the current table
      table$comparaison <- names(sub_table)[i]
      
      # Add col indication if it is from main or sub table
      table$level <- "sub"
      
      # Append to the final table
      table_to_write <- rbind(table_to_write, table)
      
    }
  }
  
  write_xlsx(x = table_to_write, path = filepath)
  
}

# Function that takes a stats excel file generated by the function above and generates
# heatmap for the pvalues

pvaluesHmap <- function(stats, selected_comparisons,
                        taxons, lvl, txn_lvl, group, displayPValues = TRUE, displayChangeArrows = FALSE, path){
  
  stats=stats[stats$level==lvl,] # Table with only stats for taxon level of interest
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
  stat_hmap$comparaison <- factor(stat_hmap$comparaison, levels = selected_comparisons)
  
  # Define color breaks
  stat_hmap$significance <- cut(stat_hmap$value,
                         breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                         labels = c("<0.001", "<0.01", "<0.05", "n.s."),
                         right = FALSE)
  
  if(displayChangeArrows){
    # Add log2 fold change information
    for(i in 1:nrow(stat_hmap)){
      matched_value <- stats$log2FoldChange[stats[[txn_lvl]] == stat_hmap[[txn_lvl]][i] & stats$comparaison == stat_hmap$comparaison[i]]
      stat_hmap$log2FoldChange[i] <- ifelse(length(matched_value) > 0, matched_value, NA)    
      }
  }
  
  # Plot heatmap with ggplot2
  p <- ggplot(stat_hmap, aes(y = .data[[txn_lvl]], x = comparaison, fill = significance)) +
    geom_tile(color = "black", lwd = 0.5, linetype = 1) +
    coord_fixed() + # Makes thing squared
    scale_fill_manual(
      values = c("n.s." = "white", "<0.05" = "#F4A3A8", "<0.01" = "#E04B54", "<0.001" = "#A40000"))+
    theme_minimal() +
    labs(x = "", y = "", fill = "P value") +
    scale_y_discrete(limits = rev(unique(stat_hmap[[txn_lvl]])))+  # Reverse y-axis order
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))  
  
  if(displayPValues){
    p <- p + geom_text(aes(label = ifelse(value=="1","NS",round(value, 5))), color = "gray", size = 4) # Show significance labels
  }
  
  if(displayChangeArrows){
    p <- p + geom_text(aes(label = ifelse(is.na(log2FoldChange),"",ifelse(log2FoldChange>0, "↑", "↓"))), color = "black", size = 6) # Show significance labels
  }

  return(p)
  
}



