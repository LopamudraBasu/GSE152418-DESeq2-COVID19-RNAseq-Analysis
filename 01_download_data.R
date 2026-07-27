# ==========================================================
# Project: Differential Gene Expression Analysis of COVID-19
# Dataset: GSE152418
# Platform: GPL24676 (Illumina NovaSeq 6000)
# Script: 01_Data_Download.R
# Author: Lopamudra Basu
# ==========================================================

#----------------------------------------------------------
# Load Required Library
#----------------------------------------------------------

library(GEOquery)

#----------------------------------------------------------
# Download GEO Dataset
#----------------------------------------------------------

gse <- getGEO(
  "GSE152418",
  GSEMatrix = TRUE
)

#----------------------------------------------------------
# Extract Sample Metadata
#----------------------------------------------------------

pheno <- pData(gse[[1]])

#----------------------------------------------------------
# Save Metadata
#----------------------------------------------------------

write.csv(
  pheno,
  file = "data/GSE152418_Metadata.csv",
  row.names = TRUE
)

#----------------------------------------------------------
# Completion Message
#----------------------------------------------------------

message("GSE152418 metadata downloaded and saved successfully.")