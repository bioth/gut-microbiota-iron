library(ggplot2)
library(ggsignif)
library(readxl)



# Generating graphs for 16s analysis with qPCR
setwd("I:/Chercheurs/Santos_Manuela/Thibault M/")
f_rodentium <- read_excel("RT-PCR F rodentium.xlsx")
colnames(f_rodentium) <- f_rodentium[3,]
f_rodentium <- f_rodentium[-c(1:3),]
meta <- read_excel("gut-microbiota-iron/experiments/ongoing exp/young-abx-exp6/dissection.xlsx")
meta$ID <- substring(meta$ID, first = 1, last = 5)
meta$gg_group <- factor(paste(meta$diet, meta$treatment, sep = ":"), levels = c("50:water", "500:water", "50:abx", "500:abx"))


data <- merge(f_rodentium, meta, by.x = "Mouse ID", by.y = "ID")
colnames(data)[4] <- "per16s"
data$per16s <- as.numeric(data$per16s)
customColors <- c("blue", "red", "darkblue", "darkred")

ggplot(data = data, aes(x = gg_group, y = per16s, color = gg_group)) +
  geom_point(size = 1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = -0.75)) + 
  
  #Error bars
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
               aes(color = gg_group),
               width = 0.2, size = 0.7,
               position = position_dodge(-0.75)) +
  
  #Mean lines
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar",
               aes(ymin = ..y.., ymax = ..y.., group = gg_group),
               color = "black", linewidth = 0.5, width = 0.5,
               position = position_dodge(-0.75))+
  
  
  labs(title = "F. rodentium 16S",
       y = "16S", color = "Groups", x = "Groups") +
  scale_color_manual(values = customColors)+
  
  # #Add significance bars
  # geom_signif(comparisons = list(c(groups[1],groups[2])),
  #             annotations = ifelse(displayPvalue, paste("p = ", 
  #                                                       format(sigtab_taxon[sigtab_taxon$comparaison == 1, "padj"], digits = 2, scientific = TRUE)),
  #                                  sigtab_taxon[sigtab_taxon$comparaison == 1, "significance"]),
  #             tip_length = 0.02,
  #             y_position =  max(relative_abundance$rel_ab)+1/12*max(relative_abundance$rel_ab),
  #             size = 1.2,  # Make the bar wider
  #             color = "black") +
  # 
  # geom_signif(comparisons = list(c(groups[3],groups[4])),
  #             annotations = ifelse(displayPvalue, paste("p = ", 
  #                                                       format(sigtab_taxon[sigtab_taxon$comparaison == 2, "padj"], digits = 2, scientific = TRUE)),
  #                                  sigtab_taxon[sigtab_taxon$comparaison == 2, "significance"]),
  #             tip_length = 0.02,
  #             y_position =  max(relative_abundance$rel_ab)+1/12*max(relative_abundance$rel_ab),
  #             size = 1.2,  # Make the bar wider
  #             color = "black") +
  # 
  # geom_signif(comparisons = list(c(groups[1],groups[3])),
  #             annotations = ifelse(displayPvalue, paste("p = ", 
  #                                                       format(sigtab_taxon[sigtab_taxon$comparaison == 3, "padj"], digits = 2, scientific = TRUE)),
  #                                  sigtab_taxon[sigtab_taxon$comparaison == 3, "significance"]),
  #             tip_length = 0.02,
  #             y_position =  max(relative_abundance$rel_ab)+2/12*max(relative_abundance$rel_ab),
  #             size = 1.2,  # Make the bar wider
  #             color = "black") +
  # 
  # geom_signif(comparisons = list(c(groups[2],groups[4])),
  #             annotations = ifelse(displayPvalue, paste("p = ", 
  #                                                       format(sigtab_taxon[sigtab_taxon$comparaison == 4, "padj"], digits = 2, scientific = TRUE)),
  #                                  sigtab_taxon[sigtab_taxon$comparaison == 4, "significance"]),
  #             tip_length = 0.02,
  #             y_position =  max(relative_abundance$rel_ab)+3/12*max(relative_abundance$rel_ab),
  #             size = 1.2,  # Make the bar wider
  #             color = "black") +
  
  theme_minimal()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Adjust title font size and style
    axis.title.x = element_text(size = 14, face = "bold"),  # Adjust x-axis label font size and style
    axis.title.y = element_text(size = 14, face = "bold"),  # Adjust y-axis label font size and style
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust x-axis tick label font size
    axis.text.y = element_text(size = 12),  # Adjust y-axis tick label font size
    legend.title = element_text(size = 12, face = "bold"),  # Remove legend title
    legend.text = element_text(size = 12),  # Adjust legend font size
    panel.grid.major = element_blank(),  # Add major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black", size = 1)) # Include axis lines  # Include axis bar





ggsave(plot = p, filename = paste(dir_taxon,"/",taxonName,"_relab.png", sep = ""), dpi = 300, height = 6, width = 6, bg = 'white')