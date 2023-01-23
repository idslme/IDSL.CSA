\name{IDSL.CSA_MSPgenerator}
\alias{IDSL.CSA_MSPgenerator}
\title{
IDSL.CSA MSP Generator
}
\description{
This function creates standard .msp files that can also be used for Pepsearch.
}
\usage{
IDSL.CSA_MSPgenerator(CSA_peaklist, msLevel, spectral_search_mode = "dda", 
spectral_search_mode_option = NA, number_processing_threads = 1)
}
\arguments{
  \item{CSA_peaklist}{
A dataframe peaklist of co-detected peaks
}
  \item{spectral_search_mode}{
Type of analysis. spectral_search_mode = c("dda", "dia", "csa")
}
  \item{msLevel}{
MS level = c(1, 2)
}
  \item{spectral_search_mode_option}{
Secondary type of analysis. spectral_search_mode_option = c(NA, "rawddaspectra", "alignedtable")
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A string of standard .msp file
}