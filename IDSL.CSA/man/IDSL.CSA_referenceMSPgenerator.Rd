\name{IDSL.CSA_referenceMSPgenerator}
\alias{IDSL.CSA_referenceMSPgenerator}
\title{
IDSL.CSA Reference MSP Generator
}
\description{
This function creates reference standard .msp files.
}
\usage{
IDSL.CSA_referenceMSPgenerator(REF_peaklist, refTable, selectedIPApeaks_IDref, msLevel,
spectral_search_mode = "dda", spectral_search_mode_option = NA)
}
\arguments{
  \item{REF_peaklist}{
A dataframe peaklist of co-detected peaks
}
  \item{refTable}{
reference CSA table
}
  \item{selectedIPApeaks_IDref}{
selectedIPApeaks_IDref
}
  \item{msLevel}{
MS level = c(1, 2)
}
  \item{spectral_search_mode}{
Type of analysis. spectral_search_mode = c("dda", "dia", "csa")
}
  \item{spectral_search_mode_option}{
Secondary type of analysis. spectral_search_mode_option = c(NA, "rawddaspectra", "alignedtable")
}
}
\value{
A string of standard .msp file
}