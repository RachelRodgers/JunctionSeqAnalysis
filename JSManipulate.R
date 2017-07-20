# JSManipulate.R

# Conditions tested in this analysis
conditionList <- list("PP_v_NN", "PT_v_NN", "PT_v_PP")

# Used to store the significant genes for each condition
rawGeneList <- list()

for(i in 1:length(conditionList)) {
        # Read in and unzip the .gz file
        rawFile <- paste(conditionList[[i]], "sigGenes.results.txt.gz", sep = "")
        unzippedDF <- readRawJunctionSeq(rawFile)
        
        # Subset based on splice junctions and "testable" status being TRUE
        sigGenesSubset <- subsetSigSJ(unzippedDF)
        
        rawGeneList[[i]] <- as.character(sigGenesSubset[ , "geneID"])
        
}

# Remove duplicates from the rawGeneList
universalGeneList <- unique(unlist(rawGeneList))

# Create data frame to plot the upset chart
upset.df = data.frame(gene_id = universalGeneList)
upset.df = mutate(upset.df, PPvNN = ifelse(gene_id %in% rawGeneList[[1]], 1, 0))
upset.df = mutate(upset.df, PTvNN = ifelse(gene_id %in% rawGeneList[[2]], 1, 0))
upset.df = mutate(upset.df, PTvPP = ifelse(gene_id %in% rawGeneList[[3]], 1, 0))
#row.names(upset.df) = upset.df$gene_id

#upset.df <- upset.df[ , -c(1)]