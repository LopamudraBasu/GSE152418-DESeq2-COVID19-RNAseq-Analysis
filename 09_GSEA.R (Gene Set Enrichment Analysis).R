# ==========================================================
# Project: Differential Gene Expression Analysis of COVID-19
# Dataset: GSE152418
# Script: 09_GSEA.R
# Author: Lopamudra Basu
# ==========================================================

#----------------------------------------------------------
# Load Libraries
#----------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

#----------------------------------------------------------
# Load Annotated DESeq2 Results
#----------------------------------------------------------

res <- readRDS("results/res.rds")

#----------------------------------------------------------
# Remove Missing Values
#----------------------------------------------------------

res <- subset(
  res,
  !is.na(log2FoldChange) &
    !is.na(GeneSymbol)
)

#----------------------------------------------------------
# Remove Duplicate Gene Symbols
#----------------------------------------------------------

res <- res[
  !duplicated(res$GeneSymbol),
]

#----------------------------------------------------------
# Convert Gene Symbols to Entrez IDs
#----------------------------------------------------------

gene_df <- bitr(
  
  res$GeneSymbol,
  
  fromType = "SYMBOL",
  
  toType = "ENTREZID",
  
  OrgDb = org.Hs.eg.db
  
)

#----------------------------------------------------------
# Merge Annotation
#----------------------------------------------------------

gene_table <- merge(
  
  gene_df,
  
  as.data.frame(res),
  
  by.x = "SYMBOL",
  
  by.y = "GeneSymbol"
  
)

#----------------------------------------------------------
# Remove Duplicate Entrez IDs
#----------------------------------------------------------

gene_table <- gene_table[
  !duplicated(gene_table$ENTREZID),
]

#----------------------------------------------------------
# Remove Missing Fold Change
#----------------------------------------------------------

gene_table <- gene_table[
  !is.na(gene_table$log2FoldChange),
]

#----------------------------------------------------------
# Create Ranked Gene List
#----------------------------------------------------------

gene_list <- gene_table$log2FoldChange

names(gene_list) <- gene_table$ENTREZID

gene_list <- sort(
  gene_list,
  decreasing = TRUE
)

#----------------------------------------------------------
# Verify Unique IDs
#----------------------------------------------------------

gene_list <- gene_list[
  !duplicated(names(gene_list))
]

#----------------------------------------------------------
# Run GSEA
#----------------------------------------------------------

gsea <- gseKEGG(
  
  geneList = gene_list,
  
  organism = "hsa",
  
  keyType = "ncbi-geneid",
  
  by = "fgsea",
  
  pvalueCutoff = 0.05,
  
  pAdjustMethod = "BH",
  
  verbose = FALSE
  
)

#----------------------------------------------------------
# Check Results
#----------------------------------------------------------

head(gsea)

summary(gsea)

#----------------------------------------------------------
# Save Results
#----------------------------------------------------------

write.csv(
  
  as.data.frame(gsea),
  
  "results/GSEA_KEGG_Results.csv",
  
  row.names = FALSE
  
)

#----------------------------------------------------------
# Dot Plot
#----------------------------------------------------------

gsea_dot <- dotplot(
  
  gsea,
  
  showCategory = 10,
  
  title = "Top Enriched KEGG Gene Sets"
  
) +
  
  labs(
    
    x = "Gene Ratio",
    
    y = NULL
    
  ) +
  
  theme_bw(base_size = 15) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold"
    )
    
  )

print(gsea_dot)

ggsave(
  
  "figures/Figure6_GSEA_Dotplot.png",
  
  gsea_dot,
  
  width = 10,
  
  height = 8,
  
  dpi = 300
  
)

#----------------------------------------------------------
# Enrichment Plot
#----------------------------------------------------------

gsea_curve <- gseaplot2(
  
  gsea,
  
  geneSetID = 1,
  
  title = gsea$Description[1]
  
)

print(gsea_curve)

ggsave(
  
  "figures/Figure7_GSEA_EnrichmentPlot.png",
  
  gsea_curve,
  
  width = 10,
  
  height = 7,
  
  dpi = 300
  
)

#----------------------------------------------------------
# Summary
#----------------------------------------------------------

cat("========================================\n")
cat("Gene Set Enrichment Analysis Completed\n")
cat("========================================\n")
cat("Results saved to:\n")
cat("results/GSEA_KEGG_Results.csv\n\n")
cat("Figures saved to:\n")
cat("Figure6_GSEA_Dotplot.png\n")
cat("Figure7_GSEA_EnrichmentPlot.png\n")
cat("========================================\n")