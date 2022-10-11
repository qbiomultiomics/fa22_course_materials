##### NOTE: You must have rna_se loaded in to run this file

##### RNA_CLINICAL
# Create rna_clinical as a data frame
rna_clinical <-  as.data.frame(rna_se@colData)

# Subset out any pertinent NA patients (in this example, pertinent refers to the age_at_index column)
age_mask <-  is.na(rna_se@colData$age_at_index) # if you are exploring a different clinical variable, remove NAs for that variable (ie race, sex, stage, etc.)
rna_clinical <- rna_clinical[!age_mask, ]

# Add descriptive row names to rna_clinical
row.names(rna_clinical) <- rna_clinical$barcode

# Subset out the normal samples
# These can be useful as control data for your final project, so there are situations in which you would want to keep this data
tissue_type_mask <- ifelse(rna_clinical$definition == "Solid Tissue Normal", F, T)
rna_clinical <- rna_clinical[tissue_type_mask, ]

# Add any calculated columns you desire:
rna_clinical$age_category <- ifelse(rna_clinical$age_at_index < 50, "Young", "Old")
rna_clinical$five_year_surv <- ifelse(is.na(ifelse(rna_clinical$days_to_death == "NA", NA,
                                                   ifelse(rna_clinical$days_to_death > 365.25*5, TRUE, FALSE))),
                                      ifelse(rna_clinical$days_to_last_follow_up > 365.25*5, TRUE, FALSE),
                                      ifelse(rna_clinical$days_to_death == "NA", NA,
                                             ifelse(rna_clinical$days_to_death > 365.25*5, TRUE, FALSE)))

##### RNA_GENES

# Create rna_genes as a data frame
rna_genes <- as.data.frame(rna_se@rowRanges@elementMetadata)

# Add descriptive row names to rna_genes
row.names(rna_genes) <- rna_genes$gene_id

##### RNA_COUNTS

# Create rna_counts as a data frame
rna_counts <- as.data.frame(rna_se@assays@data$unstranded)

# Subset out the same NA patients as done in rna_clinical using the same mask
rna_counts <- rna_counts[, !age_mask]

# If you removed normal samples from rna_clinical, do the same to rna_counts
rna_counts <- rna_counts[ , tissue_type_mask]

# Add descriptive row and column names to rna_counts
row.names(rna_counts) <- rna_genes$gene_id
colnames(rna_counts) <- rna_clinical$barcode

##### SAVING FILES AS .RDATA

# before you do this, clean up your environment with remove()
remove(age_mask, tissue_type_mask)
save.image("/Users/nicoleblack/Desktop/QBIO/qbio490_nicole/analysis_data/rnaclinical_rnagenes_rnacounts.RData") # this takes a while

##### READING IN FILES FROM .RDATA
load("/Users/nicoleblack/Desktop/QBIO/qbio490_nicole/analysis_data/rnaclinical_rnagenes_rnacounts.RData")

##### SAVING FILES AS .CSVs

# Subset out treatments, disease_type, and primary_site columns
# NOTE: We only do this so we can save as a .csv file!!
treatments_mask <- ifelse(colnames(rna_clinical) == "treatments" | colnames(rna_clinical) == "disease_type" | colnames(rna_clinical) == "primary_site", FALSE, TRUE)
rna_clinical <- rna_clinical[, treatments_mask] 

write.csv(rna_clinical, "/Users/nicoleblack/Desktop/QBIO/qbio490_nicole/brca_rna_clincial_data.csv", row.names = TRUE, col.names = TRUE)
write.csv(rna_genes, "/Users/nicoleblack/Desktop/QBIO/qbio490_nicole/brca_rna_gene_data.csv", row.names = TRUE, col.names = TRUE)
write.csv(rna_counts, "/Users/nicoleblack/Desktop/QBIO/qbio490_nicole/brca_rna_count_data.csv", row.names = TRUE, col.names = TRUE)

##### READING IN FILES FROM .CSVs
rna_clinical <- read.csv("/Users/nicoleblack/Desktop/QBIO/qbio490_nicole/brca_rna_clincial_data.csv", header = TRUE, row.names = 1)
rna_counts <- read.csv("/Users/nicoleblack/Desktop/QBIO/qbio490_nicole/brca_rna_count_data.csv", header = TRUE, row.names = 1)
rna_genes <- read.csv("/Users/nicoleblack/Desktop/QBIO/qbio490_nicole/brca_rna_gene_data.csv", header = TRUE, row.names = 1)

