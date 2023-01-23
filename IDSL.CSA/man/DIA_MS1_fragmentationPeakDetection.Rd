\name{DIA_MS1_fragmentationPeakDetection}
\alias{DIA_MS1_fragmentationPeakDetection}
\title{
CSA DIA MS1 Fragmentation Peaks Detection
}
\description{
This function detects fragmentation peaks for the Data-Independent Acquisition (DIA) analysis at  ms level 1.
}
\usage{
DIA_MS1_fragmentationPeakDetection(DIA_hrms_address, DIA_hrms_file, peaklist,
selectedIPApeaks, massError, smoothingWindowMS1, scanTolerance, nSpline,
topRatioPeakHeight, intensityThresholdFragment, pearsonRHOthreshold, outputDIAeic = NULL,
number_processing_threads = 1)
}
\arguments{
  \item{DIA_hrms_address}{
path to the HRMS file
}
  \item{DIA_hrms_file}{
DIA HRMS file
}
  \item{peaklist}{
IDSL.IPA peaklist
}
  \item{selectedIPApeaks}{
A vector of selected IDSL.IPA peaks only when a number of IDSL.IPA peaks from one peaklist is processed.
}
  \item{massError}{
Mass accuracy in Da
}
  \item{smoothingWindowMS1}{
number of scans for peak smoothing.
}
  \item{scanTolerance}{
a scan tolerance to extend the chromatogram for better calculations.
}
  \item{nSpline}{
number of points for further smoothing using a cubic spline smoothing method to add more points to calculate Pearson correlation rho values
}
  \item{topRatioPeakHeight}{
The top percentage of the chromatographic peak to calculate Pearson correlation rho values
}
  \item{intensityThresholdFragment}{
a value to represent intensity threshold for the fragment at the apex chromatogram scan
}
  \item{pearsonRHOthreshold}{
Minimum threshold for Pearson correlation rho values
}
  \item{outputDIAeic}{
When `NULL` DIA EICs are not plotted. `outputDIAeic` represents an address to save DIA EICs figures.
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A dataframe peaklist of co-detected DIA analysis.
}
\references{
Fakouri Baygi, S., Fernando, S., Hopke, P.K., Holsen, T.M., Crimmins, B.S. (2021). Nontargeted discovery of novel contaminants in the Great Lakes region: A comparison of fish fillets and fish consumers. \emph{Environmental Science & Technology}, 55(6), 3765-3774, \doi{10.1021/acs.est.0c08507}
}