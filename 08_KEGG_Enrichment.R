# ==========================================================
# Project: Differential Gene Expression Analysis of COVID-19
# Dataset: GSE152418
# Platform: GPL24676 (Illumina NovaSeq 6000)
# Script: 08_KEGG_Pathway_Analysis.R
# Author: Lopamudra Basu
# ==========================================================

#----------------------------------------------------------
# Load Required Libraries
#----------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

#----------------------------------------------------------
# Load Significant Differentially Expressed Genes
#----------------------------------------------------------

sig_res <- read.csv(
  "results/Significant_DEGs_Annotated.csv",
  row.names = 1
)

#----------------------------------------------------------
# Remove Missing Gene Symbols
#----------------------------------------------------------

sig_res <- subset(
  sig_res,
  !is.na(GeneSymbol)
)

#----------------------------------------------------------
# Convert Gene Symbols to Entrez IDs
#----------------------------------------------------------

gene_ids <- bitr(
  sig_res$GeneSymbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

#----------------------------------------------------------
# Perform KEGG Pathway Enrichment Analysis
#----------------------------------------------------------

kegg <- enrichKEGG(
  
  gene = gene_ids$ENTREZID,
  
  organism = "hsa",
  
  keyType = "ncbi-geneid",
  
  pvalueCutoff = 0.05,
  
  pAdjustMethod = "BH"
  
)

#----------------------------------------------------------
# Convert Entrez IDs to Gene Symbols
#----------------------------------------------------------

kegg <- setReadable(
  kegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

#----------------------------------------------------------
# View Results
#----------------------------------------------------------

head(kegg)

summary(kegg)

#----------------------------------------------------------
# Save KEGG Results
#----------------------------------------------------------

write.csv(
  
  as.data.frame(kegg),
  
  file = "results/KEGG_Pathway_Enrichment.csv",
  
  row.names = FALSE
  
)

#----------------------------------------------------------
# KEGG Dot Plot
#----------------------------------------------------------
kegg_dotplot <- dotplot(
  
  kegg,
  
  showCategory = 10,
  
  title = "Top Enriched KEGG Pathways"
  
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
    
    axis.title = element_text(face = "bold")
    
  )

  

print(kegg_dotplot)

ggsave(
  
  filename = "figures/Figure5_KEGG_Dotplot.png",
  
  plot = kegg_dotplot,
  
  width = 10,
  
  height = 8,
  
  dpi = 300
  
)

#----------------------------------------------------------
# KEGG Bar Plot
#----------------------------------------------------------

kegg_barplot <- barplot(
  
  kegg,
  
  showCategory = 10,
  
  title = "Top KEGG Pathways"
  
) +
  
  theme_bw(base_size = 15) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    )
    
  )

print(kegg_barplot)

ggsave(
  
  filename = "figures/Figure6_KEGG_Barplot.png",
  
  plot = kegg_barplot,
  
  width = 10,
  
  height = 8,
  
  dpi = 300
  
)

#----------------------------------------------------------
# Summary
#----------------------------------------------------------

cat("========================================\n")
cat("KEGG Pathway Analysis Completed\n")
cat("========================================\n")
cat("Results saved to:\n")
cat("results/KEGG_Pathway_Enrichment.csv\n")
cat("\n")
cat("Figures saved to:\n")
cat("figures/Figure5_KEGG_Dotplot.png\n")
cat("figures/Figure6_KEGG_Barplot.png\n")
cat("========================================\n")