\name{DDA_fragmentationPeakDetection}
\alias{DDA_fragmentationPeakDetection}
\title{
DDA Fragmentation Peaks Detection
}
\description{
This function detects fragmentation peaks for the Data-Dependent Acquisition (DDA) analysis.
}
\usage{
DDA_fragmentationPeakDetection(DDA_hrms_address, DDA_hrms_file, peaklist,
selectedIPApeaks, massErrorPrecursor, DDAprocessingMode = 'MostIntenseDDAspectra',
outputDDAspectra = NULL, number_processing_threads = 1)
}
\arguments{
  \item{DDA_hrms_address}{
path to the HRMS file
}
  \item{DDA_hrms_file}{
DDA HRMS file
}
  \item{peaklist}{
IDSL.IPA peaklist
}
  \item{selectedIPApeaks}{
A vector of selected IDSL.IPA peaks only when a number of IDSL.IPA peaks from one peaklist is processed.
}
  \item{massErrorPrecursor}{
Mass accuracy (Da) to find precursor m/z in IDSL.IPA peaklists
}
  \item{DDAprocessingMode}{
c('MostIntenseDDAspectra', c('DDAspectraIntegration', massErrorIntegration), c('IonFiltering', massErrorIonFiltering, minPercentageDetectedScans, rsdCutoff, pearsonRHOthreshold)). Required variables for each DDA processing mode should be provided in this vector.
}
  \item{outputDDAspectra}{
When `NULL` DDA spectra are not plotted. `outputDDAspectra` represents an address to save DDA spectra figures.
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A dataframe peaklist of co-detected DDA analysis.
}