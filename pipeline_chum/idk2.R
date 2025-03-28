library(readxl)
library(ggplot2)
library(car)

setwd("I:/Chercheurs/Santos_Manuela/Thibault M")
ferrozineT35 <- readxl::read_excel("Ferrozine T35.xlsx")
ferrozineT35 = ferrozineT35[-c(1,2,27),]
colnames(ferrozineT35)[14:16] = list("iron_concentration","treatment","diet")
ferrozineT35$gg_group = factor(paste(ferrozineT35$diet, ferrozineT35$treatment, sep = "_"), levels =  c("50_water", "50_dss", "500_water", "500_dss"))
levels(ferrozineT35$gg_group)
custom_colors = list("blue","red")
ferrozineT35$iron_concentration = as.numeric(ferrozineT35$iron_concentration)

max(ferrozineT35$iron_concentration[ferrozineT35$gg_group == "50_water"])

ggplot(data = ferrozineT35, aes(x = factor(gg_group),  y = iron_concentration, color = as.character(diet))) +
  
  stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..-0.25, yend=..y.., color = diet), size =1)+ #adding horizontal bars representing means
  stat_summary(fun="mean", geom = "segment", mapping=aes(xend=..x..+0.25, yend=..y.., color = diet), size =1)+ 
  stat_summary(aes(color = diet), fun.data="mean_cl_normal", geom="errorbar", width=0.2, size = 1) + #adding SEM error bars
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), alpha = 0.5) +
  
  # geom_signif(comparisons =   desired_order1,
  #             map_signif_level = TRUE, # Pour afficher l'étoile
  #             y_position = 9, # Ajuste cette valeur selon ton graphique
  #             annotations = ifelse(p_value < 0.05, "*", "ns")) + # Affiche * si significatif +
  
  
  labs(title = "Iron concentration after 5 weeks of iron exposure",
       x = "Diet + Treatment",
       y = "Iron concentration in stools",
       color = "Diet")+ 
  scale_color_discrete(labels = c("50 ppm FeSO4", "500 ppm FeSO4"))+
  scale_color_manual(values = custom_colors)+
  scale_x_discrete(labels = c("50 ppm + water", "50 ppm + DSS", "500 ppm + water", "500 ppm + DSS"))+
  theme_minimal() +  
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1)
  )

#stats

leveneTest()

group1 = ferrozineT35[ferrozineT35$gg_group == "50_water",]
group2 = ferrozineT35[ferrozineT35$gg_group == "50_dss",]

print(leveneTest(group1$iron_concentration, group2$iron_concentration))
print(shapiro.test(group2$iron_concentration))
