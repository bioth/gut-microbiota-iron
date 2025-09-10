library("dada2"); packageVersion("dada2")
library(ShortRead)
library(data.table)

# Load asv table for each microbiota 18 and microbiota 19 runs
asv_table_M18 <- as.data.frame(fread("~/projects/def-santosmm/jr34106/Microbiota_18_19_merged/asv_table/asv_table_M18.csv", sep = ";"))
rownames(asv_table_M18) <- asv_table_M18[,1]  # Use the first column as row names
asv_table_M18 <- asv_table_M18[,-1]  # Drop the first column

asv_table_M19 <- as.data.frame(fread("~/projects/def-santosmm/jr34106/Microbiota_18_19_merged/asv_table/asv_table_M19_2.csv", sep = ";"))
rownames(asv_table_M19) <- asv_table_M19[,1]  # Use the first column as row names
asv_table_M19 <- asv_table_M19[,-1]  # Drop the first column

asv_table_M18 <- as.matrix(asv_table_M18)
asv_table_M19 <- as.matrix(asv_table_M19)

# Merge the ASV table by summing similar sequences ids
asv_table <- mergeSequenceTables(asv_table_M18, asv_table_M19, repeats="sum")

# Removing chimeras
seqtab.nochim <- removeBimeraDenovo(asv_table, method="consensus", multithread=TRUE, verbose=TRUE,
                                    allowOneOff=FALSE, minFoldParentOverAbundance=8)

# What percentage of chimeras over the total dataset
print("Percentage of chimeras over the total dataset:")
1-sum(seqtab.nochim)/sum(asv_table)

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "~/projects/def-santosmm/jr34106/data/training_set/silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)
taxa_w_species <- addSpecies(taxtab = taxa, refFasta = "~/projects/def-santosmm/jr34106/data/training_set/silva_species_assignment_v138.1.fa.gz")

write.table(seqtab.nochim, sep = ";", file = "~/projects/def-santosmm/jr34106/Microbiota_18_19_merged/asv_table/merged_asv_table.csv", col.names = TRUE)
write.table(taxa_w_species, sep = ";", file = "~/projects/def-santosmm/jr34106/Microbiota_18_19_merged/taxonomy/merged_taxa_annotation.csv", col.names = TRUE)