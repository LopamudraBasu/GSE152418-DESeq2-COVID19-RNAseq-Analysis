# ==========================================================
# Project: Differential Gene Expression Analysis of COVID-19
# Dataset: GSE152418
# Platform: GPL24676 (Illumina NovaSeq 6000)
# Script: 03_Differential_Expression_Analysis.R
# Author: Lopamudra Basu
# ==========================================================

#----------------------------------------------------------
# Load Required Libraries
#----------------------------------------------------------

library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)

#----------------------------------------------------------
# Load Preprocessed DESeq2 Dataset
#----------------------------------------------------------

dds <- readRDS(
  "results/dds_filtered.rds"
)

#----------------------------------------------------------
# Run Differential Expression Analysis
#----------------------------------------------------------

dds <- DESeq(dds)

#----------------------------------------------------------
# Extract Differential Expression Results
#----------------------------------------------------------

res <- results(dds)

#----------------------------------------------------------
# Save Complete Results
#----------------------------------------------------------

saveRDS(
  res,
  file = "results/res.rds"
)

write.csv(
  as.data.frame(res),
  file = "results/All_DEGs.csv"
)

#----------------------------------------------------------
# Filter Significant Differentially Expressed Genes (DEGs)
#----------------------------------------------------------

sig_res <- subset(
  res,
  padj < 0.05 &
    abs(log2FoldChange) > 1
)

# Sort by adjusted p-value
sig_res <- sig_res[
  order(sig_res$padj),
]

#----------------------------------------------------------
#----------------------------------------------------------
# Annotate All Genes with Gene Symbols
#----------------------------------------------------------

res$GeneSymbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Replace missing Gene Symbols with Ensembl IDs
res$GeneSymbol[
  is.na(res$GeneSymbol)
] <- rownames(res)[
  is.na(res$GeneSymbol)
]

#----------------------------------------------------------
# Save Complete Annotated Results
#----------------------------------------------------------

saveRDS(
  res,
  file = "results/res.rds"
)

#----------------------------------------------------------
# Extract Significant DEGs
#----------------------------------------------------------

sig_res <- subset(
  res,
  padj < 0.05 &
    abs(log2FoldChange) > 1
)

sig_res <- sig_res[
  order(sig_res$padj),
]

#----------------------------------------------------------
# Save Significant DEGs
#----------------------------------------------------------

write.csv(
  as.data.frame(sig_res),
  file = "results/Significant_DEGs_Annotated.csv"
)
#----------------------------------------------------------
# Summary
#----------------------------------------------------------

cat("=====================================\n")
cat("Differential Expression Analysis Completed\n")
cat("=====================================\n")
cat("Total genes analyzed      :", nrow(res), "\n")
cat("Significant DEGs detected :", nrow(sig_res), "\n")
cat("=====================================\n")