\name{CSA_alignedPeaksTanimotoCoefficientCalculator}
\alias{CSA_alignedPeaksTanimotoCoefficientCalculator}
\title{
CSA Aligned Peaks Tanimoto Coefficient Calculator
}
\description{
This function groups co-detected peaks on the aligned table.
}
\usage{
CSA_alignedPeaksTanimotoCoefficientCalculator(address_input_msp, peakXcol,
minPercenetageDetection = 5, minNumberFragments = 2, minTanimotoCoefficient = 0.1,
RTtolerance = 0.05, number_processing_threads = 1)
}
\arguments{
  \item{address_input_msp}{
address of the .msp files generated via IDSL.IPA CSA aggregation
}
  \item{peakXcol}{
aligned indexed peak table
}
  \item{minPercenetageDetection}{
minimum CSA frequency detection
}
  \item{minNumberFragments}{
minimum frequency
}
  \item{minTanimotoCoefficient}{
minimum Tanimoto coefficient
}
  \item{RTtolerance}{
retention time tolerance to detect common peaks
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A 2-column dataframe of `Co-detectedIDs` and `TanimotoCoefficients`
}