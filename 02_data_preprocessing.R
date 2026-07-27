# ==========================================================
# Project: Differential Gene Expression Analysis of COVID-19
# Dataset: GSE152418
# Platform: GPL24676 (Illumina NovaSeq 6000)
# Script: 02_Data_Preprocessing.R
# Author: Lopamudra Basu
# ==========================================================

#----------------------------------------------------------
# Load Required Libraries
#----------------------------------------------------------

library(DESeq2)

#----------------------------------------------------------
# Load Metadata
#----------------------------------------------------------

pheno <- read.csv(
  "data/GSE152418_Metadata.csv",
  row.names = 1,
  check.names = FALSE
)

#----------------------------------------------------------
# Load Raw Count Matrix
#----------------------------------------------------------

counts <- read.delim(
  "data/GSE152418_p20047_Study1_RawCounts.txt.gz",
  row.names = 1,
  check.names = FALSE
)

#----------------------------------------------------------
# Match Sample Names
#----------------------------------------------------------

rownames(pheno) <- pheno$title

stopifnot(
  all(colnames(counts) == rownames(pheno))
)

#----------------------------------------------------------
# Remove Convalescent Sample
#----------------------------------------------------------

pheno_filtered <- pheno[
  pheno$`disease state:ch1` != "Convalescent",
]

counts_filtered <- counts[
  ,
  rownames(pheno_filtered)
]

stopifnot(
  all(colnames(counts_filtered) == rownames(pheno_filtered))
)

#----------------------------------------------------------
# Create Experimental Groups
#----------------------------------------------------------

pheno_filtered$condition <- factor(
  pheno_filtered$`disease state:ch1`,
  levels = c("Healthy", "COVID-19")
)

#----------------------------------------------------------
# Create DESeq2 Dataset
#----------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = pheno_filtered,
  design = ~ condition
)

#----------------------------------------------------------
# Filter Lowly Expressed Genes
# Keep genes with at least 10 counts in at least 3 samples
#----------------------------------------------------------

keep <- rowSums(counts(dds) >= 10) >= 3

dds <- dds[keep, ]

#----------------------------------------------------------
# Save Processed Objects
#----------------------------------------------------------

saveRDS(
  dds,
  file = "results/dds_filtered.rds"
)

write.csv(
  pheno_filtered,
  file = "results/Filtered_Metadata.csv",
  row.names = TRUE
)

write.csv(
  counts_filtered,
  file = "results/Filtered_Count_Matrix.csv"
)

#----------------------------------------------------------
# Summary
#----------------------------------------------------------

cat("=====================================\n")
cat("Data Preprocessing Completed\n")
cat("=====================================\n")
cat("Samples retained :", ncol(dds), "\n")
cat("Genes retained   :", nrow(dds), "\n")
cat("=====================================\n")