# ==========================================================
# Project: Differential Gene Expression Analysis of COVID-19
# Dataset: GSE152418
# Platform: GPL24676 (Illumina NovaSeq 6000)
# Script: 06_Heatmap.R
# Author: Lopamudra Basu
# ==========================================================

#----------------------------------------------------------
# Load Required Libraries
#----------------------------------------------------------

library(DESeq2)
library(pheatmap)
library(RColorBrewer)

#----------------------------------------------------------
# Load DESeq2 Dataset
#----------------------------------------------------------

dds <- readRDS(
  "results/dds_filtered.rds"
)

#----------------------------------------------------------
# Load Annotated Differential Expression Results
#----------------------------------------------------------

res <- readRDS(
  "results/res.rds"
)

#----------------------------------------------------------
# Variance Stabilizing Transformation (VST)
#----------------------------------------------------------

vsd <- vst(
  dds,
  blind = FALSE
)

#----------------------------------------------------------
# Extract Normalized Expression Matrix
#----------------------------------------------------------

expression_matrix <- assay(vsd)

#----------------------------------------------------------
# Remove Genes with Missing Adjusted P-values
#----------------------------------------------------------

res <- subset(
  res,
  !is.na(padj)
)

#----------------------------------------------------------
# Sort by Adjusted P-value
#----------------------------------------------------------

res <- res[
  order(res$padj),
]

#----------------------------------------------------------
# Select Top 50 Differentially Expressed Genes
#----------------------------------------------------------

top50_genes <- head(
  rownames(res),
  50
)

#----------------------------------------------------------
# Subset Expression Matrix
#----------------------------------------------------------

heatmap_matrix <- expression_matrix[
  top50_genes,
]

#----------------------------------------------------------
# Replace Ensembl IDs with Gene Symbols
#----------------------------------------------------------

gene_symbols <- res$GeneSymbol[
  match(
    top50_genes,
    rownames(res)
  )
]

gene_symbols[
  is.na(gene_symbols)
] <- top50_genes[
  is.na(gene_symbols)
]

rownames(heatmap_matrix) <- gene_symbols

#----------------------------------------------------------
# Scale Expression (Row-wise Z-score)
#----------------------------------------------------------

heatmap_matrix <- t(
  scale(
    t(heatmap_matrix)
  )
)

#----------------------------------------------------------
# Sample Annotation
#----------------------------------------------------------

annotation_col <- data.frame(
  
  Condition = colData(dds)$condition
  
)

rownames(annotation_col) <- colnames(heatmap_matrix)

#----------------------------------------------------------
# Annotation Colors
#----------------------------------------------------------

annotation_colors <- list(
  
  Condition = c(
    
    Healthy = "#377EB8",
    
    `COVID-19` = "#E41A1C"
    
  )
  
)

#----------------------------------------------------------
# Create Heatmap
#----------------------------------------------------------

heatmap_plot <- pheatmap(
  
  heatmap_matrix,
  
  color = colorRampPalette(
    
    rev(
      brewer.pal(11, "RdBu")
    )
    
  )(100),
  
  cluster_rows = TRUE,
  
  cluster_cols = TRUE,
  
  scale = "none",
  
  annotation_col = annotation_col,
  
  annotation_colors = annotation_colors,
  
  show_colnames = FALSE,
  
  show_rownames = TRUE,
  
  fontsize_row = 5.5,
  
  fontsize = 10,
  
  border_color = NA,
  
  main = "Expression Profile of Top 50 Differentially Expressed Genes
  COVID-19 vs Healthy PBMC"
  
)

#----------------------------------------------------------
# Save Heatmap
#----------------------------------------------------------

png(
  filename = "figures/Figure3_Heatmap_Top50_DEGs.png",
  width = 3600,
  height = 4200,
  res = 300
)

pheatmap(
  
  heatmap_matrix,
  
  color = colorRampPalette(
    
    rev(
      brewer.pal(11, "RdBu")
    )
    
  )(100),
  
  cluster_rows = TRUE,
  
  cluster_cols = TRUE,
  
  scale = "none",
  
  annotation_col = annotation_col,
  
  annotation_colors = annotation_colors,
  
  show_colnames = FALSE,
  
  show_rownames = TRUE,
  
  fontsize_row = 7,
  
  fontsize = 10,
  
  border_color = NA,
  
  main = "Expression Profile of Top 50 Differentially Expressed Genes
  COVID-19 vs Healthy PBMC"
  
  
)

dev.off()

#----------------------------------------------------------
# Summary
#----------------------------------------------------------

cat("========================================\n")
cat("Heatmap Generated Successfully\n")
cat("========================================\n")
cat("Top DEGs Displayed : 50\n")
cat("Figure Saved As:\n")
cat("figures/Figure3_Heatmap_Top50_DEGs.png\n")
cat("========================================\n")