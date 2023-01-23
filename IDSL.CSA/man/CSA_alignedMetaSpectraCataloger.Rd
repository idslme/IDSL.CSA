\name{CSA_alignedMetaSpectraCataloger}
\alias{CSA_alignedMetaSpectraCataloger}
\title{
CSA Aligned Meta-Spectra Cataloger
}
\description{
This function generates integrated and most abundant aligned spectra from the aligned spectra
}
\usage{
CSA_alignedMetaSpectraCataloger(address_input_msp, peakXcol, peak_height,
CSA_aligned_property_table, groupedID, minTanimotoCoefficient = 0.5,
number_processing_threads = 1)
}
\arguments{
  \item{address_input_msp}{
address of the .msp files generated via IDSL.IPA DIA grouping
}
  \item{peakXcol}{
aligned indexed peak table
}
  \item{peak_height}{
aligned peak height table
}
  \item{CSA_aligned_property_table}{
a matrix of three columns of "IPA detection frequency", "median_height", and "median_R13C" for the aligned peak table
}
  \item{groupedID}{
A 2-column dataframe of `Co-detectedIDs` and `TanimotoCoefficients` from the `CSA_alignedPeaksTanimotoCoefficientCalculator` module
}
  \item{minTanimotoCoefficient}{
minimum Tanimoto coefficient
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A list of two objects for `MSP_integrated_aligned_spectra` and `MSP_most_abundant_aligned_spectra`
}