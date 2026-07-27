<div align="center">

# 🧬 Differential Gene Expression Analysis of COVID-19 PBMC using DESeq2

### Comprehensive RNA-seq Transcriptomic Analysis of GEO Dataset **GSE152418**

**DESeq2 • Functional Enrichment • GSEA • STRINGdb • Cytoscape**

![R](https://img.shields.io/badge/R-4.4+-276DC3?style=for-the-badge&logo=r)
![Bioconductor](https://img.shields.io/badge/Bioconductor-DESeq2-green?style=for-the-badge)
![RNAseq](https://img.shields.io/badge/RNA--seq-Transcriptomics-blue?style=for-the-badge)
![GEO](https://img.shields.io/badge/Database-GEO-orange?style=for-the-badge)
![Status](https://img.shields.io/badge/Project-Completed-success?style=for-the-badge)

---

*A reproducible end-to-end RNA-seq analysis workflow for identifying differentially expressed genes, enriched biological pathways, and protein interaction networks associated with COVID-19.*

</div>

---

# 📖 Overview

This project presents a complete RNA-seq differential gene expression analysis of **Peripheral Blood Mononuclear Cells (PBMCs)** collected from COVID-19 patients and healthy controls using the publicly available **GEO dataset GSE152418**.

The analysis was performed entirely in **R** using **Bioconductor** packages and follows a reproducible workflow commonly used in transcriptomics research.

The project includes:

- RNA-seq preprocessing
- Differential gene expression analysis
- Principal Component Analysis (PCA)
- Volcano plot visualization
- Heatmap generation
- GO Biological Process enrichment
- KEGG pathway enrichment
- Gene Set Enrichment Analysis (GSEA)
- Protein–Protein Interaction (PPI) network construction
- Cytoscape network visualization
- Hub gene identification

---

# ⭐ Project Highlights

✅ Complete RNA-seq Analysis Pipeline

✅ Differential Gene Expression using DESeq2

✅ Functional Enrichment Analysis

✅ Gene Set Enrichment Analysis (GSEA)

✅ STRINGdb Protein Interaction Network

✅ Cytoscape Network Visualization

✅ Publication-quality Figures

✅ Fully Reproducible R Workflow

---

# 📂 Dataset

| Feature | Description |
|----------|-------------|
| Database | NCBI GEO |
| Dataset | GSE152418 |
| Organism | Homo sapiens |
| Platform | GPL24676 |
| Technology | Illumina NovaSeq 6000 |
| Sample Type | Peripheral Blood Mononuclear Cells (PBMCs) |

---

# 🎯 Objectives

- Identify differentially expressed genes (DEGs)
- Compare transcriptomic profiles between COVID-19 and healthy individuals
- Identify enriched biological processes and signaling pathways
- Investigate immune-related molecular mechanisms
- Construct protein interaction networks
- Identify highly connected hub genes

---

# 🔬 Analysis Workflow

```text
GEO Dataset Download
          │
          ▼
Metadata Processing
          │
          ▼
Raw Count Matrix Filtering
          │
          ▼
DESeq2 Differential Expression Analysis
          │
          ▼
Principal Component Analysis
          │
          ▼
Volcano Plot
          │
          ▼
Heatmap
          │
          ▼
GO Enrichment
          │
          ▼
KEGG Enrichment
          │
          ▼
Gene Set Enrichment Analysis
          │
          ▼
STRINGdb Protein Network
          │
          ▼
Cytoscape Visualization
          │
          ▼
Hub Gene Identification
```

---

# 💻 Software & Packages

### Programming

- R

### Bioconductor

- GEOquery
- DESeq2
- clusterProfiler
- enrichplot
- org.Hs.eg.db
- STRINGdb

### CRAN

- ggplot2
- ggrepel
- pheatmap
- igraph
- dplyr
- RColorBrewer

### Network Analysis

- Cytoscape

---

# 📁 Repository Structure

```
GSE152418-DESeq2-COVID19-RNAseq-Analysis

│── data
│── scripts
│── figures
│── results
└── README.md
```

---

# 📊 Results

## Principal Component Analysis

![](figures/Figure1_PCA.png)

The PCA demonstrates clear transcriptomic separation between COVID-19 patients and healthy controls.

---

## Differential Expression Volcano Plot

![](figures/Figure2A_Volcano.png)

Genes with significant differential expression are highlighted using adjusted p-values and log₂ fold-change thresholds.

---

## Heatmap of Top Differentially Expressed Genes

![](figures/Figure3_Heatmap_Top50_DEGs.png)

Expression patterns of the top 50 differentially expressed genes clearly distinguish COVID-19 samples from healthy controls.

---

## GO Biological Process Enrichment

![](figures/Figure4_GO_BP_Dotplot.png)

Functional enrichment analysis identified immune-related biological processes significantly associated with COVID-19.

---

## KEGG Pathway Enrichment

![](figures/Figure5_KEGG_Dotplot.png)

KEGG analysis revealed pathways involved in antiviral response and immune signaling.

---

## Protein–Protein Interaction Network

![](figures/Figure8_PPI_Network.png)

STRINGdb and Cytoscape were used to construct a high-confidence interaction network and identify hub genes.

---

# 🧬 Biological Insights

This analysis revealed substantial transcriptomic alterations in COVID-19 PBMC samples.

Key observations include:

- Strong activation of interferon-mediated antiviral responses.
- Significant enrichment of immune-related biological processes.
- Identification of highly connected hub genes within the protein interaction network.
- Clear separation of disease and healthy samples based on global gene expression profiles.

Representative genes identified include:

- IFI27
- ISG15
- MX1
- IFI44L

---

# 🧠 Skills Demonstrated

- RNA-seq Data Analysis
- Differential Gene Expression
- Transcriptomics
- Functional Genomics
- Pathway Enrichment Analysis
- Network Biology
- Protein–Protein Interaction Analysis
- Cytoscape Visualization
- Reproducible Research
- Scientific Data Visualization

---

# 🚀 Reproducibility

Run the scripts sequentially:

```
01_Data_Download.R

02_Data_Preprocessing.R

03_Differential_Expression_Analysis.R

04_PCA_Visualization.R

05A_Top_DEGs_Volcano.R

05B_COVID_Genes_Volcano.R

06_Heatmap.R

07_GO_Enrichment.R

08_KEGG_Pathway_Analysis.R

09_GSEA.R

10_Protein_Protein_Interaction.R
```

---

# 📦 Output

The project generates:

- Differentially Expressed Gene Tables
- PCA Plot
- Volcano Plots
- Heatmap
- GO Enrichment Results
- KEGG Enrichment Results
- GSEA Results
- STRING Protein Interaction Network
- Cytoscape Files
- Hub Gene Table

---

# 🔮 Future Improvements

- WGCNA
- Immune Cell Infiltration Analysis
- Machine Learning Classification
- Multi-omics Integration
- Survival Analysis
- Biomarker Validation

---

# 👩‍💻 Author

**Lopamudra Basu**

**M.Sc. Biotechnology**

 Research Assistant

**Areas of Interest**

- Transcriptomics
- RNA-seq
- Cancer Genomics
- Functional Genomics
- Multi-omics
- Computational Biology

GitHub: https://github.com/LopamudraBasu

---

<div align="center">

⭐ **If you found this repository useful, please consider giving it a star!**

</div>
