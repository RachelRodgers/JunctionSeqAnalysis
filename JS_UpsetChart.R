# JS_UpsetChart.R

# Load libraries(assuming they're already installed)
library(plyr)
library(dplyr)
library(magrittr)
library(UpSetR)

# PP_v_NN
PPvNN.unzipped <- readRawJunctionSeq("PP_v_NNsigGenes.results.txt.gz")
PPvNN.genes <- subsetSigGenes(PPvNN.unzipped)
PPvNN.genes <- as.vector(PPvNN.genes[ , "geneID"])

# PT_v_NN
PTvNN.unzipped <- readRawJunctionSeq("PT_v_NNsigGenes.results.txt.gz")
PTvNN.genes <- subsetSigGenes(PTvNN.unzipped)
PTvNN.genes <- as.vector(PTvNN.genes[ , "geneID"])

# PT_v_PP
PTvPP.unzipped <- readRawJunctionSeq("PT_v_PPsigGenes.results.txt.gz")
PTvPP.genes <- subsetSigGenes(PTvPP.unzipped)
PTvPP.genes <- as.vector(PTvPP.genes[ , "geneID"])

universal.gene.list <- unique( c(PPvNN.genes, PTvNN.genes, PTvPP.genes))

# Can't get mutate() to work with piping, need to figure out...
upset.df <- data.frame(gene_id = universal.gene.list) 
upset.df <- mutate(upset.df, PPvNN = ifelse(gene_id %in% PPvNN.genes, 1, 0)) 
upset.df <- mutate(upset.df, PTvNN = ifelse(gene_id %in% PTvNN.genes, 1, 0)) 
upset.df <- mutate(upset.df, PTvPP = ifelse(gene_id %in% PTvPP.genes, 1, 0))

row.names(upset.df) = upset.df$gene_id

upset.df <- upset.df[ , -c(1)]

pdf("JunctionSeq_UpsetChart.pdf", 11, 8.5)
upset( upset.df, text.scale=c( 3, 3, 3, 2, 3, 3), point.size = 5, line.size = 2 )
dev.off()