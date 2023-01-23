\name{DDA_rawSpectraDeconvolution}
\alias{DDA_rawSpectraDeconvolution}
\title{
DDA Raw Spectra Deconvolution
}
\description{
This function stacks all DDA scans.
}
\usage{
DDA_rawSpectraDeconvolution(DDA_hrms_address, DDA_hrms_file, rawDDAspectraVar = NULL,
number_processing_threads = 1)
}
\arguments{
  \item{DDA_hrms_address}{
path to the HRMS file
}
  \item{DDA_hrms_file}{
DDA HRMS file
}
  \item{rawDDAspectraVar}{
c(NULL, list(precursorMZvec, precursorRTvec, massError, RTtolerance)). When NULL, all scans with precursor values are used for DDA peaklist generation. When the list is provided, it filters the scans with respect to predefined `precursorMZvec` and `precursorRTvec` values.
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A dataframe stacked DDA scans.
}