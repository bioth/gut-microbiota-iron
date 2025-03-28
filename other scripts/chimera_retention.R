# Analysis of chimeras distribution across Microbiota_18 dataset

library(readxl)
library(ggplot2)

df <- read_excel("Downloads/Microbiota_18.chimeraMetrics.xlsx") # Load data from genome quebec
df$id <- substring(df$Sample, first = 1, last = 5)
df$timepoint <- substring(df$Sample, first = 7, last = length(df$Sample))
meta <- read.csv("Documents/CHUM_git/Microbiota_18/metadata/metadata.csv", sep = ";")
meta <- meta[-46,]
meta$id <- substring(meta$id, first = 1, last = 5)

merged_df <- merge(df, meta, by = "id") # Merge metadata with chimera metrics data
merged_df$timepoint <- factor(merged_df$timepoint, levels = c("T0","T35","T49","d53","T54","Tfinal"))
merged_df$gg_group <- paste(merged_df$diet, merged_df$treatment, sep = ":")
merged_df$gg_group <- factor(merged_df$gg_group, levels = c("50:water","500:water","50:dss","500:dss"))

ggplot(data = merged_df, aes(x = timepoint, y = `Overall Retention %`, fill = gg_group)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.8), alpha = 0.5) +
  geom_point(aes(colour = gg_group), position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8), alpha = 0.5)+
  scale_color_manual(values = c("blue", "red", "darkblue", "darkred"))+
  scale_fill_manual(values = c("blue", "red", "darkblue", "darkred"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

df <- read_excel("Documents/CHUM_git/Microbiota_18/r_console_output/retention stats/microbiota_18_own_retention_stats.xlsx") # Load my data
df$id <- substring(df$Sample, first = 1, last = 5)
df$timepoint <- substring(df$Sample, first = 7, last = length(df$Sample))
merged_df <- merge(df, meta, by = "id") # Merge metadata with chimera metrics data
merged_df$timepoint <- factor(merged_df$timepoint, levels = c("T0","T35","T49","d53","T54","Tfinal"))
merged_df$gg_group <- paste(merged_df$diet, merged_df$treatment, sep = ":")
merged_df$gg_group <- factor(merged_df$gg_group, levels = c("50:water","500:water","50:dss","500:dss"))
merged_df$`Overall Retention %` <- merged_df$nonchim/merged_df$merged*100
merged_df$`Merged Retention %` <- merged_df$merged/merged_df$FilteredF*100

ggplot(data = merged_df, aes(x = timepoint, y = `Overall Retention %`, fill = gg_group)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.8), alpha = 0.5) +
  geom_point(aes(colour = gg_group), position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8), alpha = 0.5)+
  scale_color_manual(values = c("blue", "red", "darkblue", "darkred"))+
  scale_fill_manual(values = c("blue", "red", "darkblue", "darkred"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggplot(data = merged_df, aes(x = timepoint, y = `Merged Retention %`, fill = gg_group)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.8), alpha = 0.5) +
  geom_point(aes(colour = gg_group), position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8), alpha = 0.5)+
  scale_color_manual(values = c("blue", "red", "darkblue", "darkred"))+
  scale_fill_manual(values = c("blue", "red", "darkblue", "darkred"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
