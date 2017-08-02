#JS_VolcanoPlot.R

# source file containing functions:  JSFunctions2.R
# source('C:/users/rache/Desktop/repos/JunctionSeqAnalysis/JSFunctions2.R')

# Load libraries
library( ggplot2 )
library( ggrepel )
library( biomaRt )
library( dplyr )

#-------------------------------------------------------------------------------------#
# Build initial data frame from raw output file

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

#-------------------------------------------------------------------------------------#
# Query BioMart to get gene symbols

# Make list of genes for querying Biomart
geneList <- as.vector(JSResults.forVolcano$geneID)

# Query Biomart to get gene symbols
ensemblQuery <- getGeneSymbols(geneList)
# For ensemblQuery df, change any blank cells to NA
ensemblQuery[ensemblQuery == "df"] <- NA

#-------------------------------------------------------------------------------------#
# Modify the data frame to prepare for plotting

# Merge the PPvNN.forVolcano df with the ensemblQuery df on the geneID column
output.df <- merge(JSResults.forVolcano, ensemblQuery, by.x = "geneID", by.y = "ensembl_gene_id")

# Change class of the featureID column for further editing
output.df$featureID <- as.character(output.df$featureID)
# Create a new column in the output.df that is a combo of both the
# hgnc_symbol and the feature ID (J###) for easier identification on plot:

# Split featureID column on the colon
featureList <- strsplit(output.df$featureID, ":")

# Add the second element of each list within featureList to a new list
newFeatureList <- list()
for (i in 1:length(featureList)) {
        newFeatureList[i] <- featureList[[i]][[2]]
}

newFeatureList <- unlist(newFeatureList)

# Add this list to the dataframe
output.df <- cbind(output.df, "featureSymbol" = newFeatureList)

# Create new column for the final name using the newly-create "featureSymbol" column
# and the hgnc_symbol column
output.df$geneFeature <- paste(output.df$hgnc_symbol, output.df$featureSymbol, sep = "_")

#-------------------------------------------------------------------------------------#
# Create volcano plotting recipe and generate the plot

## ggplot2 formatting recipes
gg_bigger_texts = theme(
        axis.title = element_text( size=22 ),
        axis.text = element_text( size=20 ),
        legend.text = element_text(size=14 ),
        legend.title = element_text(size=15 ),
        plot.title = element_text( size=22 ),
        strip.text = element_text( size=15 )
)

gg_no_legend = theme(
        legend.position='none'
)

gg_no_grid = theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
)

gg_no_x_grid = theme(
        panel.grid.major.x = element_blank() )

gg_no_y_grid = theme(
        panel.grid.major.y = element_blank() )

gg_center_title = theme(
        plot.title = element_text( hjust = 0.5 )
)


volcano_plot = output.df %>%
        mutate( ., Significant = ifelse( abs(log2FC.Pre.Treated_Psoriasis.Post.Treated_Psoriasis.) > 1 & padjust < 0.05 , "Yes", "No" ) ) %>%
        ggplot( ., aes( x=log2FC.Pre.Treated_Psoriasis.Post.Treated_Psoriasis., y=-log10( padjust ), colour=Significant ) ) +
        scale_colour_manual( values = c( "gray", "red" ) ) +
        theme_bw() +
        gg_bigger_texts +
        gg_no_legend +
        gg_no_grid +
        gg_center_title +
        geom_point( size=1 ) +
        geom_hline( yintercept=1.30, linetype=2 ) +
        geom_vline( xintercept=c( -1, 1 ), linetype=2 ) +
        geom_text_repel( data=subset( output.df, padjust < 0.05 & abs( log2FC.Pre.Treated_Psoriasis.Post.Treated_Psoriasis. ) > 1.0 )[c(1:30),], colour="black", aes( label=geneFeature ), size=3 ) +
        xlab("\nlog2( Fold change )") +
        ylab("-log10( padjust )\n")

volcano_plot


