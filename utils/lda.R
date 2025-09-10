library(microbiomeMarker)

multipleGroupsLDA <- function(ps, deseq, taxon, gg_group, pairs, customColors){
  
  asvList <- row.names(tax_table(ps_samuel)[!is.na(tax_table(ps_samuel)[,"Species"]),]) # List of ASVs which were identified at the taxonomic level of interest
  
  # Perform LEfSe analysis
  lefse_results <- run_lefse(
    ps = ps_samuel,
    group = "gg_group",   # Replace with your grouping variable
    kw_cutoff = 0.05,                # Kruskal-Wallis test p-value cutoff
    lda_cutoff = 2,                  # LDA score cutoff
    multigrp_strat = TRUE            # Enable multi-group comparison
  )
  
  # Plot LDA scores
  plot_lda(lefse_results)
  plot_ef_bar(lefse_results)
  microbiomeMarker::plot
    
  # Apply rlog transformation to the count data
  rlog_data <- rlog(deseq_samuel, blind = TRUE)
  
  # Extract the log-transformed counts for the significant features
  lda_data <- assay(rlog_data)[asvList, ]
  
  # Combine the data with the group labels (e.g., condition)
  lda_input <- merge(t(lda_data), sample_data(ps_samuel)[,"gg_group"], by = "row.names")
  row.names(lda_input) <- lda_input$Row.names
  lda_input <- lda_input[-1]
  
  
  train_control <- trainControl(method = "cv", number = 10)
  # Train the LDA model using caret (this returns a train object)
  lda_model <- train(as.formula(paste(varToCompare, "~ .")),
                     data = lda_input,
                     method = "lda",
                     trControl = train_control)  # ensure you have defined train_control
  
  # Extract ASV-level LDA coefficients (effect sizes)
  lda_effects <- data.frame(
    ASV = rownames(lda_model$finalModel$scaling),
    LDA_Score = lda_model$finalModel$scaling[,1],  # Use first discriminant
    taxa = paste0(substring(tax_table(ps)[rownames(lda_model$finalModel$scaling),"Genus"], 1, 1), ". ",
                  tax_table(ps)[rownames(lda_model$finalModel$scaling),"Species"])
  )
  
  lda_summary <- lda_effects %>%
    group_by(taxa) %>%
    summarise(mean_LDA = mean(LDA_Score)) %>%
    mutate(Group = ifelse(mean_LDA > 0, levels(sample_data(ps)[[varToCompare]])[2], levels(sample_data(ps)[[varToCompare]])[1]))  # Separate step for clarity
  
  print(lda_summary)
  
  # Plot species-level mean LDA scores
  p <- ggplot(lda_summary, aes(x = reorder(taxa, mean_LDA), y = mean_LDA, fill = Group)) +
    geom_col() +
    coord_flip() +
    geom_hline(yintercept = 0, color = "black") +
    scale_fill_manual(values = customColors) +
    labs(x = "", y = "LDA Score") +
    theme_minimal() +
    theme(legend.position = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(face = "bold"),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold.italic")
    )
  
  ggsave(plot = p, filename = paste0(dir, "/lda.png"), dpi = 300, bg = "white", width = 4, height = 6)
}


sample_data(ps_samuel)$gg_group <- factor(sample_data(ps_samuel)$gg_group, levels = c("Wt:Vehicle", "Wt:Putrescine", "IL-22ra1-/-:Vehicle", "IL-22ra1-/-:Putrescine"))  # Vehicle as reference

#Creates ps subset for taxonomical level of interest
ps_subset <- tax_glom(ps_samuel, taxrank = "Species")

#Differential abundance Samuel
deseq_samuel <- phyloseq_to_deseq2(ps_samuel, ~ genotype + treatment+ genotype : treatment) # Full formula with interaction term

#Setting "Wt" as the baseline for genotype
colData(deseq_samuel)$genotype <- relevel(colData(deseq_samuel)$genotype, ref="Wt")

#Setting "Vehicle" as the baseline for treatment
colData(deseq_samuel)$treatment <- relevel(colData(deseq_samuel)$treatment, ref="Vehicle")

deseq_samuel <- DESeq(deseq_samuel, test="Wald", fitType = "parametric")

resultsNames(deseq_samuel)

customColors = list('black','#A22004',"#AB8F23","#04208D")
pairs <- list(list("Wt:Vehicle","Wt:Putrescine"), list("IL-22ra1-/-:Vehicle","IL-22ra1-/-:Putrescine"), list("Wt:Vehicle","IL-22ra1-/-:Vehicle"), list("Wt:Putrescine","IL-22ra1-/-:Putrescine"))
