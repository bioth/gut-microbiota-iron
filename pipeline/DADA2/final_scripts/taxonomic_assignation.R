library("dada2"); packageVersion("dada2")
library(data.table)

# Load asv table for each microbiota 18 and microbiota 19 runs
asv_table_M18 <- as.data.frame(fread("~/Documents/CHUM_git/Microbiota_18_final/Microbiota18_final_data2/asv_table/asv_table.csv", sep = ";"))
rownames(asv_table_M18) <- asv_table_M18[,1]  # Use the first column as row names
asv_table_M18 <- asv_table_M18[,-1]  # Drop the first column
asv_table_M18 <- as.matrix(asv_table_M18)
table(nchar(getSequences(asv_table_M18)))
ncol(asv_table_M18)

# Assign taxonomy
taxa <- assignTaxonomy(asv_table_M18, "~/Documents/CHUM_git/training_set/silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)
taxa_w_species <- addSpecies(taxtab = taxa, refFasta = "~/Documents/CHUM_git/training_set/silva_species_assignment_v138.1.fa.gz")
write.table(taxa_w_species, sep = ";", file = "~/Documents/CHUM_git/Microbiota_18_final/taxonomy/taxa_annotation_final.csv", col.names = TRUE)

asv_table_M19 <- as.data.frame(fread("~/Documents/CHUM_git/Microbiota_19/Microbiota_19_final_data2/asv_table/asv_table.csv", sep = ";"))
rownames(asv_table_M19) <- asv_table_M19[,1]  # Use the first column as row names
asv_table_M19 <- asv_table_M19[,-1]  # Drop the first column
asv_table_M19 <- as.matrix(asv_table_M19)
ncol(asv_table_M19)
table(nchar(getSequences(asv_table_M19)))

taxa <- assignTaxonomy(asv_table_M19, "~/Documents/CHUM_git/training_set/silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)
taxa_w_species <- addSpecies(taxtab = taxa, refFasta = "~/Documents/CHUM_git/training_set/silva_species_assignment_v138.1.fa.gz")
write.table(taxa_w_species, sep = ";", file = "~/Documents/CHUM_git/Microbiota_19/taxonomy/taxa_annotation_final.csv", col.names = TRUE)
