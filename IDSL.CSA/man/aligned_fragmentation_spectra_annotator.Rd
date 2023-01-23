\name{aligned_fragmentation_spectra_annotator}
\alias{aligned_fragmentation_spectra_annotator}
\title{
Aligned Fragmentation Spectra Annotator
}
\description{
This function detects frequent matched compounds across multiple samples on the aligned peak table matrix.
}
\usage{
aligned_fragmentation_spectra_annotator(PARAM_AT, output_path)
}
\arguments{
  \item{PARAM_AT}{
a parameter driven from the `CSA_AlignedTable_xlsxAnalyzer` module.
}
  \item{output_path}{
output path
}
}
\value{
This function stores `.Rdata` and `.csv` data from dataframe of aligned fragmentation spectra.
}