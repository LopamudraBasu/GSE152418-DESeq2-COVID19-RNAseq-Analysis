# ==========================================================
# Project: Differential Gene Expression Analysis of COVID-19
# Dataset: GSE152418
# Platform: GPL24676 (Illumina NovaSeq 6000)
# Script: 07_GO_Enrichment.R
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
# Perform GO Enrichment Analysis
#----------------------------------------------------------

go_bp <- enrichGO(
  
  gene = gene_ids$ENTREZID,
  
  OrgDb = org.Hs.eg.db,
  
  keyType = "ENTREZID",
  
  ont = "BP",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff = 0.05,
  
  qvalueCutoff = 0.05,
  
  readable = TRUE
  
)

#----------------------------------------------------------
# View Results
#----------------------------------------------------------

head(go_bp)

summary(go_bp)

#----------------------------------------------------------
# Save GO Enrichment Results
#----------------------------------------------------------

write.csv(
  
  as.data.frame(go_bp),
  
  file = "results/GO_BP_Enrichment.csv",
  
  row.names = FALSE
  
)

#----------------------------------------------------------
# Dot Plot
#----------------------------------------------------------

go_dotplot <- dotplot(
  
  go_bp,
  
  showCategory = 10,
  
  title = "GO Biological Process Enrichment"
  
) +
  
  theme_bw(base_size = 14) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 18
    ),
    
    axis.title = element_text(
      face = "bold"
    )
    
  )

print(go_dotplot)

ggsave(
  
  filename = "figures/Figure4_GO_BP_Dotplot.png",
  
  plot = go_dotplot,
  
  width = 10,
  
  height = 8,
  
  dpi = 300
  
)


#----------------------------------------------------------
# Summary
#----------------------------------------------------------

cat("========================================\n")
cat("GO Enrichment Analysis Completed\n")
cat("========================================\n")
cat("Results saved to:\n")
cat("results/GO_BP_Enrichment.csv\n")
cat("\n")
cat("Figures saved to:\n")
cat("figures/Figure4_GO_BP_Dotplot.png\n")
cat("========================================\n")