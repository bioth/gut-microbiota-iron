library("dada2"); packageVersion("dada2")
library(ShortRead)
library(data.table)

# Load asv table for each microbiota 18 and microbiota 19 runs
asv_table_M18 <- as.matrix(read.csv("~/Documents/CHUM_git/Microbiota_18_final/asv_table/asv_table.csv", sep = ";"))
asv_table_M19 <- as.matrix(read.csv("~/Documents/CHUM_git/Microbiota_19/asv_table/asv_table.csv", sep = ";"))
ncol(asv_table_M18) # Check number of ASVs
ncol(asv_table_M19)

# Merge the ASV table by summing similar sequences ids
asv_table <- mergeSequenceTables(asv_table_M18, asv_table_M19, repeats="sum")

# Assign taxonomy
taxa <- assignTaxonomy(asv_table, "~/Documents/CHUM_git/training_set/silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)
taxa_w_species <- addSpecies(taxtab = taxa, refFasta = "~/Documents/CHUM_git/training_set/silva_species_assignment_v138.1.fa.gz")

existingDirCheck("~/Documents/CHUM_git/Microbiota_18_19_merged")
# Save taxa matrix and merged asv table
write.table(asv_table, sep = ";", file = "~/Documents/CHUM_git/Microbiota_18_19_merged/asv_table.csv", col.names = TRUE)
write.table(taxa_w_species, sep = ";", file = "~/Documents/CHUM_git/Microbiota_18_19_merged/taxa_annotation.csv", col.names = TRUE)