# JSFunctions2.R

# Read in the .gz file
readRawJunctionSeq <- function(fileName) {
        zippedFile <- gzfile(fileName, open = 'r')
        fileContents<- read.table(zippedFile, header = TRUE) #Read in as data frame
        close(zippedFile)
        return(fileContents)
}


# Exclude records from the raw data (previously read in via readRawJunctionSeq)
# based on testable, featureType, and padjust columns
subsetSigGenes <- function(junctionSeq_rawDF) {
        junctionSeq_sigDF <- subset(junctionSeq_rawDF, 
                                    subset = (testable != FALSE & 
                                              featureType != "exonic_part" & 
                                              padjust < 0.05))
        return(junctionSeq_sigDF)
}

# Get hgnc symbol and ensembl id for given list of ensembl ids 
# for human, version 85
getGeneSymbols <- function(gene.list) {
  ensembl <- useMart(host = 'jul2016.archive.ensembl.org',
                     biomart = 'ENSEMBL_MART_ENSEMBL',
                     dataset = 'hsapiens_gene_ensembl')
  query <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                 filters = 'ensembl_gene_id', values = gene.list,
                 mart = ensembl)
  return(query)
}