# JSFunctions.R

# Read in the .gz file
readRawJunctionSeq <- function(f, mode = "r", header = TRUE) {
        f_zipped <- gzfile(f, open = mode)
        f_unzipped <- read.table(f_zipped, header = header) #Read in as data frame
        close(f_zipped)
        return(f_unzipped)
}

# Exclude records from the raw data (read in via readRawJunctionSeq)
# based on testable and featureType columns
subsetSigSJ <- function(x) {
        x_subset <- subset(x, subset = x$testable != FALSE & x$featureType != "exonic_part")
        return(x_subset)
}

