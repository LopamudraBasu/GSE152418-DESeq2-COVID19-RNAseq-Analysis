# ==========================================================
# Project: Differential Gene Expression Analysis of COVID-19
# Dataset: GSE152418
# Platform: GPL24676 (Illumina NovaSeq 6000)
# Script: 05_Volcano_Plot.R
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
# Convert DESeq2 Results to Data Frame
#----------------------------------------------------------

res_df <- as.data.frame(res)

#----------------------------------------------------------
# Annotate Gene Symbols
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

res_df <- na.omit(res_df)

#----------------------------------------------------------
# Categorize Genes
#----------------------------------------------------------

res_df$Category <- "Not Significant"

res_df$Category[
  res_df$padj < 0.05 &
    res_df$log2FoldChange > 1
] <- "Upregulated"

res_df$Category[
  res_df$padj < 0.05 &
    res_df$log2FoldChange < -1
] <- "Downregulated"

#----------------------------------------------------------
# Select Top 10 Significant Genes
#----------------------------------------------------------

top_genes <- head(
  res_df[
    order(res_df$padj),
  ],
  10
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
      "Downregulated" = "#377EB8",
      "Not Significant" = "grey75",
      "Upregulated" = "#E41A1C"
    )
  ) +
  
  scale_size_manual(
    values = c(
      "Downregulated" = 2.5,
      "Not Significant" = 1.8,
      "Upregulated" = 2.5
    ),
    guide = "none"
  ) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.6
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.6
  ) +
  
  geom_text_repel(
    data = top_genes,
    aes(label = GeneSymbol),
    size = 4,
    color = "black",
    box.padding = 0.6,
    point.padding = 0.5,
    force = 2,
    min.segment.length = 0,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  
  labs(
    title = "Differential Gene Expression",
    subtitle = "COVID-19 vs Healthy PBMC (GSE152418)",
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~Adjusted~P~Value),
    color = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 20,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 14,
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 16
    ),
    
    axis.text = element_text(
      size = 13
    ),
    
    legend.position = "right",
    
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
  filename = "figures/Figure2A_Volcano.png",
  plot = volcano_plot,
  width = 9,
  height = 6.5,
  dpi = 300
)

#----------------------------------------------------------
# Summary
#----------------------------------------------------------

cat("=====================================\n")
cat("Volcano Plot Generated Successfully\n")
cat("=====================================\n")
cat("Figure saved to: figures/Figure2A_Volcano.png\n")
cat("=====================================\n")