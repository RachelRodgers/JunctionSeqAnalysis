# JSFunctions.R

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

