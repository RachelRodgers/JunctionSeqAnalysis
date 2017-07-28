#JS_VolcanoPlot.R

# source file containing functions
source('C:/users/rache/Desktop/repos/JunctionSeqAnalysis/JSFunctions2.R')

# Load libraries
library( ggplot2 )
library( ggrepel )
library( biomaRt )
library( dplyr )


JSFile.contents <- readRawJunctionSeq("PT_v_PPallGenes.results.txt.gz")

# From the full dataset, keep only features that were deemed testable, 
# are junctions, and have padjust < 0.95
JSResults.keepFeatures <- filter(JSFile.contents, 
                                 testable != FALSE,
                                 featureType != "exonic_part",
                                 padjust < 0.95)

# Extract only those columns needed for plotting
# featureID, geneID, pvalue, padjust, log2FC...
JSResults.forVolcano <- JSResults.keepFeatures[c(1,2,12,13,21)]

# Make list of genes for querying Biomart
geneList <- as.vector(JSResults.forVolcano$geneID)

# Query Biomart to get gene symbols
ensemblQuery <- getGeneSymbols(geneList)?
# For ensemblQuery df, change any blank cells to NA
ensemblQuery[ensemblQuery == "df"] <- NA

# Merge the PPvNN.forVolcano df with the ensemblQuery df on the geneID column
output.df <- merge(JSResults.forVolcano, ensemblQuery, by.x = "geneID", by.y = "ensembl_gene_id")



