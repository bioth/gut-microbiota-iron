require("dada2"); packageVersion("dada2")
require(data.table)

# Load ASV table
asv_table <- as.data.frame(fread("path/to/asv/table", sep = ";"))
rownames(asv_table) <- asv_table[,1]  # Use the first column as row names
asv_table <- asv_table[,-1]  # Drop the first column
asv_table <- as.matrix(asv_table)
table(nchar(getSequences(asv_table)))
ncol(asv_table)

# Assign taxonomy
taxa <- assignTaxonomy(asv_table, "path/to/silva/training/set/silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)
taxa_w_species <- addSpecies(taxtab = taxa, refFasta = "path/to/silva/species/training/set/silva_species_assignment_v138.1.fa.gz")