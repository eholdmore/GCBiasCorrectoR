#### GATK-style GC Bias Correction ####
### Erica M. Holdmore ###
## based on GCBiasCorrector.java by David Benjamin and Samuel Lee
## https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/copynumber/denoising/GCBiasCorrector.java

# Use simple regression to learn multiplicative correction factors as a function of GC content.
# Obtained regression curve by filling 1% GC-content bins with the coverages from genomic intervals 
# that correspond to each GC bin and take the median to get a robust estimate of the curve. 
# Convolve medians with exponential curve to smooth out bins with few data (i.e. extreme GC values that occur rarely).

# GC bins are 0%, 1%, ..., 100%
NUMBER_OF_GC_BINS <- 101

# scale (in units of GC content from 0 to 1) over which GC bias correlation decays
correlationLength <- 0.02
correlationDecayRatePerBin <- 1 / (correlationLength * NUMBER_OF_GC_BINS)

# Apache Commons median doesn't work on empty arrays; this value is a placeholder to avoid exceptions
DUMMY_VALUE_NEVER_USED <- 1.0

# Function to learn multiplicative correction factors as a function of GC content
GCBiasCorrector <- function(readCounts, intervalGCContent) {
  stopifnot(length(intervalGCContent) > 0)
  stopifnot(length(intervalGCContent) == length(readCounts))
  
  # Initialize a list to hold read counts for each GC bin
  readCountsByGC <- vector("list", length = NUMBER_OF_GC_BINS)
  for (i in seq_along(readCountsByGC)) {
    readCountsByGC[[i]] <- list()
  }
  
  # Assign read counts to appropriate GC bins
  for (i in seq_along(intervalGCContent)) {
    binIndex <- gcContentToBinIndex(intervalGCContent[i])
    readCountsByGC[[binIndex + 1]] <- c(readCountsByGC[[binIndex + 1]], readCounts[i])
  }
  
  # Function to calculate correction factors
  calculateCorrectionFactors <- function(readCountsByGC) {
    medians <- sapply(readCountsByGC, function(x) {
      if (length(x) > 0) {
        median(unlist(x))
      } else {
        DUMMY_VALUE_NEVER_USED
      }
    })
    
    correctionFactors <- sapply(1:NUMBER_OF_GC_BINS, function(bin) {
      weights <- sapply(1:NUMBER_OF_GC_BINS, function(n) {
        length(readCountsByGC[[n]]) * exp(-abs(bin - n) * correlationDecayRatePerBin)
      })
      weightedMedian <- sum(medians * weights) / sum(abs(weights))
      return(1 / weightedMedian)
    })
    
    return(correctionFactors)
  }
  
  # Calculate correction factors
  gcCorrectionFactors <- calculateCorrectionFactors(readCountsByGC)
  
  return(gcCorrectionFactors)
}

# Function to map GC content to bin index
gcContentToBinIndex <- function(gcContent) {
  return(round(gcContent * (NUMBER_OF_GC_BINS - 1)))
}

# Function to calculate the corrected coverage
#correctedCoverage <- function(readCount, gcContent, gcCorrectionFactors) {
#  binIndex <- gcContentToBinIndex(gcContent)
#  return(gcCorrectionFactors[binIndex + 1] * readCount)
#}

# Function to correct GC bias of read counts
correctGCBias <- function(readCounts, intervalGCContent) {
  stopifnot(ncol(readCounts) == length(intervalGCContent))
  
  totalCoveragePerSample <- rowSums(readCounts)
  gcBiasCorrectors <- lapply(seq_len(nrow(readCounts)), function(i) {
    GCBiasCorrector(readCounts[i, ], intervalGCContent)
  })
  
  # Correct the input counts in-place
  for (i in seq_len(nrow(readCounts))) {
    for (j in seq_len(ncol(readCounts))) {
      readCounts[i, j] <- gcBiasCorrectors[[i]][gcContentToBinIndex(intervalGCContent[j]) + 1] * readCounts[i, j]
    }
  }
  
  # Normalize the corrected counts
  sampleNormalizationFactors <- totalCoveragePerSample / rowSums(readCounts)
  readCounts <- t(t(readCounts) * sampleNormalizationFactors)
  
  return(readCounts)
}