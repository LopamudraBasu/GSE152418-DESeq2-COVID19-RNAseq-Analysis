# ==========================================================
# Project: Differential Gene Expression Analysis of COVID-19
# Dataset: GSE152418
# Platform: GPL24676 (Illumina NovaSeq 6000)
# Script: 05B_Volcano_Plot_COVID_Genes.R
# Author: Lopamudra Basu
# ==========================================================

#----------------------------------------------------------
# Load Required Libraries
#----------------------------------------------------------

library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(AnnotationDbi)

#----------------------------------------------------------
# Load DESeq2 Results
#----------------------------------------------------------

res <- readRDS(
  "results/res.rds"
)

#----------------------------------------------------------
# Convert Results to Data Frame
#----------------------------------------------------------

res_df <- as.data.frame(res)

#----------------------------------------------------------
# Annotate Ensembl IDs with Gene Symbols
#----------------------------------------------------------

res_df$GeneSymbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res_df),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Replace missing Gene Symbols with Ensembl IDs
res_df$GeneSymbol[
  is.na(res_df$GeneSymbol)
] <- rownames(res_df)[
  is.na(res_df$GeneSymbol)
]

#----------------------------------------------------------
# Remove Missing Values
#----------------------------------------------------------

res_df <- subset(
  res_df,
  !is.na(log2FoldChange) &
    !is.na(padj)
)

#----------------------------------------------------------
# Classify Genes
#----------------------------------------------------------

res_df$Category <- "Not Significant"

res_df$Category[
  res_df$padj < 0.05 &
    res_df$log2FoldChange > 1
] <- "Upregulated (COVID-19)"

res_df$Category[
  res_df$padj < 0.05 &
    res_df$log2FoldChange < -1
] <- "Downregulated (COVID-19)"

#----------------------------------------------------------
# COVID-19 Signature Genes
#----------------------------------------------------------

genes_to_label <- c(
  "IFI27",
  "ISG15",
  "MX1",
  "IFI44L"
)

label_genes <- subset(
  res_df,
  GeneSymbol %in% genes_to_label &
    padj < 0.05 &
    abs(log2FoldChange) > 1
)

#----------------------------------------------------------
# Create Volcano Plot
#----------------------------------------------------------

volcano_plot <- ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj)
  )
) +
  
  geom_point(
    aes(
      color = Category,
      size = Category
    ),
    alpha = 0.65
  ) +
  
  scale_color_manual(
    values = c(
      "Downregulated (COVID-19)" = "#377EB8",
      "Not Significant" = "grey75",
      "Upregulated (COVID-19)" = "#E41A1C"
    )
  ) +
  
  scale_size_manual(
    values = c(
      "Downregulated (COVID-19)" = 2.5,
      "Not Significant" = 1.8,
      "Upregulated (COVID-19)" = 2.5
    ),
    guide = "none"
  ) +
  
  geom_vline(
    xintercept = c(-1,1),
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey50"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey50"
  ) +
  
  geom_text_repel(
    data = label_genes,
    aes(label = GeneSymbol),
    color = "black",
    size = 5,
    force = 6,
    box.padding = 0.8,
    point.padding = 0.6,
    segment.color = "grey40",
    segment.size = 0.4,
    max.overlaps = Inf
  ) +
  
  labs(
    title = "COVID-19 Signature Genes",
    subtitle = "Differential Expression in PBMC (GSE152418)",
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~Adjusted~P~Value),
    color = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 20,
      hjust = 0.5,
      margin = margin(b = 8)
    ),
    
    plot.subtitle = element_text(
      size = 14,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 16
    ),
    
    axis.text = element_text(
      size = 13
    ),
    
    legend.position = "right",
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 12
    )
  )

#----------------------------------------------------------
# Display Plot
#----------------------------------------------------------

print(volcano_plot)

#----------------------------------------------------------
# Save Figure
#----------------------------------------------------------

ggsave(
  filename = "figures/Figure2B_COVID_Genes_Volcano.png",
  plot = volcano_plot,
  width = 9,
  height = 6.5,
  dpi = 300
)

#----------------------------------------------------------
# Summary
#----------------------------------------------------------

cat("========================================\n")
cat("COVID-19 Signature Volcano Plot Complete\n")
cat("========================================\n")
cat("Output : figures/Figure2B_COVID_Genes_Volcano.png\n")
cat("========================================\n")