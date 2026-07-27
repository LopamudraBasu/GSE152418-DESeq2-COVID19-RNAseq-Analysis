# ==========================================================
# Project: Differential Gene Expression Analysis of COVID-19
# Dataset: GSE152418
# Script: 10_Protein_Protein_Interaction.R
# Author: Lopamudra Basu
# ==========================================================

#----------------------------------------------------------
# Load Libraries
#----------------------------------------------------------

library(STRINGdb)
library(org.Hs.eg.db)
library(dplyr)
library(igraph)

#----------------------------------------------------------
# Load Significant DEGs
#----------------------------------------------------------

deg <- read.csv(
  "results/Significant_DEGs_Annotated.csv",
  stringsAsFactors = FALSE
)

#----------------------------------------------------------
# Remove Missing Gene Symbols
#----------------------------------------------------------

deg <- deg %>%
  filter(!is.na(GeneSymbol))

#----------------------------------------------------------
# Remove Duplicate Gene Symbols
#----------------------------------------------------------

deg <- deg %>%
  distinct(GeneSymbol, .keep_all = TRUE)

#----------------------------------------------------------
# Keep Highly Significant Genes
#----------------------------------------------------------

deg <- deg %>%
  filter(
    padj < 0.05,
    abs(log2FoldChange) > 2
  )

#----------------------------------------------------------
# Select Top 100 Genes
#----------------------------------------------------------

deg <- deg %>%
  arrange(padj)

top100 <- head(deg, 100)

cat("Top genes selected:", nrow(top100), "\n")

#----------------------------------------------------------
# Initialize STRING Database
#----------------------------------------------------------

string_db <- STRINGdb(
  
  version = "12",
  
  species = 9606,
  
  score_threshold = 700,
  
  input_directory = ""
  
)

#----------------------------------------------------------
# Map Gene Symbols to STRING IDs
#----------------------------------------------------------

mapped <- string_db$map(
  
  top100,
  
  "GeneSymbol",
  
  removeUnmappedRows = TRUE
  
)

cat("Mapped genes:", nrow(mapped), "\n")

#----------------------------------------------------------
# Retrieve Protein Interactions
#----------------------------------------------------------

ppi <- string_db$get_interactions(
  
  mapped$STRING_id
  
)

cat("Protein interactions:", nrow(ppi), "\n")

#----------------------------------------------------------
# Keep Only Interactions Between Selected Genes
#----------------------------------------------------------

ppi <- ppi %>%
  
  filter(
    
    from %in% mapped$STRING_id,
    
    to %in% mapped$STRING_id
    
  )

cat("Filtered interactions:", nrow(ppi), "\n")

#----------------------------------------------------------
# Save Interaction Table
#----------------------------------------------------------

write.csv(
  
  ppi,
  
  "results/STRING_PPI_Network.csv",
  
  row.names = FALSE
  
)

#----------------------------------------------------------
# Build Network
#----------------------------------------------------------

network <- graph_from_data_frame(
  
  ppi,
  
  directed = FALSE
  
)

#----------------------------------------------------------
# Calculate Degree
#----------------------------------------------------------

degree_values <- degree(network)

hub_genes <- data.frame(
  
  STRING_ID = names(degree_values),
  
  Degree = degree_values
  
)

#----------------------------------------------------------
# Add Gene Symbols
#----------------------------------------------------------

hub_genes <- merge(
  
  hub_genes,
  
  mapped[, c("STRING_id","GeneSymbol")],
  
  by.x = "STRING_ID",
  
  by.y = "STRING_id"
  
)

hub_genes <- hub_genes %>%
  
  arrange(desc(Degree))

#----------------------------------------------------------
# Save Hub Genes
#----------------------------------------------------------

write.csv(
  
  hub_genes,
  
  "results/Hub_Genes.csv",
  
  row.names = FALSE
  
)

#----------------------------------------------------------
# Export Cytoscape Edge Table
#----------------------------------------------------------

edges <- data.frame(
  
  Source = ppi$from,
  
  Target = ppi$to
  
)

write.table(
  
  edges,
  
  "results/Cytoscape_Edges.tsv",
  
  sep = "\t",
  
  quote = FALSE,
  
  row.names = FALSE
  
)

#----------------------------------------------------------
# Export Cytoscape Node Table
#----------------------------------------------------------

nodes <- mapped %>%
  
  select(
    
    STRING_id,
    
    GeneSymbol,
    
    log2FoldChange,
    
    padj
    
  )

nodes$Degree <- degree_values[
  
  match(
    
    nodes$STRING_id,
    
    names(degree_values)
    
  )
  
]

write.table(
  
  nodes,
  
  "results/Cytoscape_Nodes.tsv",
  
  sep = "\t",
  
  quote = FALSE,
  
  row.names = FALSE
  
)

#----------------------------------------------------------
# Display Top Hub Genes
#----------------------------------------------------------

cat("\n")

cat("========== TOP HUB GENES ==========\n")

print(
  
  head(
    
    hub_genes,
    
    20
    
  )
  
)

cat("\n")

cat("Analysis Completed Successfully\n")