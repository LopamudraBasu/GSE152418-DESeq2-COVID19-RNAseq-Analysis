# ==========================================================
# Project: Differential Gene Expression Analysis of COVID-19
# Dataset: GSE152418
# Platform: GPL24676 (Illumina NovaSeq 6000)
# Script: 04_PCA_Visualization.R
# Author: Lopamudra Basu
# ==========================================================

#----------------------------------------------------------
# Load Required Libraries
#----------------------------------------------------------

library(DESeq2)
library(ggplot2)

#----------------------------------------------------------
# Load Preprocessed DESeq2 Dataset
#----------------------------------------------------------

dds <- readRDS(
  "results/dds_filtered.rds"
)

#----------------------------------------------------------
# Variance Stabilizing Transformation (VST)
#----------------------------------------------------------

vsd <- vst(
  dds,
  blind = FALSE
)

#----------------------------------------------------------
# Principal Component Analysis (PCA)
#----------------------------------------------------------

pcaData <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pcaData, "percentVar")
)

#----------------------------------------------------------
# Create Publication-Quality PCA Plot
#----------------------------------------------------------

p <- ggplot(
  pcaData,
  aes(
    PC1,
    PC2,
    color = condition
  )
) +
  
  geom_point(
    size = 5
  ) +
  
  stat_ellipse(
    level = 0.95,
    linewidth = 0.8
  ) +
  
  labs(
    title = "Principal Component Analysis of RNA-seq Samples",
    subtitle = "COVID-19 vs Healthy PBMC (GSE152418)",
    x = paste0(
      "PC1 (",
      percentVar[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      percentVar[2],
      "%)"
    ),
    color = "Condition"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 12
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    legend.title = element_text(
      face = "bold"
    )
  )

#----------------------------------------------------------
# Display Plot
#----------------------------------------------------------

print(p)

#----------------------------------------------------------
# Save Figure
#----------------------------------------------------------

ggsave(
  filename = "figures/Figure1_PCA.png",
  plot = p,
  width = 9,
  height = 6,
  dpi = 300
)

#----------------------------------------------------------
# Summary
#----------------------------------------------------------

cat("=====================================\n")
cat("PCA Visualization Completed\n")
cat("=====================================\n")
cat("PC1 Variance :", percentVar[1], "%\n")
cat("PC2 Variance :", percentVar[2], "%\n")
cat("=====================================\n")