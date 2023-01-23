\name{CSA_fragmentationPeakDetection}
\alias{CSA_fragmentationPeakDetection}
\title{
CSA peakList MSP generation
}
\description{
This function detects fragmentation peaks for the Composite Spectra Analysis (CSA) using IDSL.IPA peaklists.
}
\usage{
CSA_fragmentationPeakDetection(CSA_hrms_address, CSA_hrms_file,
tempAlignedTableSubsetsFolder = NULL, peaklist, selectedIPApeaks = NULL,
RTtolerance, massError, minSNRbaseline, smoothingWindowMS1, scanTolerance, nSpline,
topRatioPeakHeight, minIonRangeDifference, minNumCSApeaks, pearsonRHOthreshold,
outputCSAeic = NULL)
}
\arguments{
  \item{CSA_hrms_address}{
path to the HRMS file
}
  \item{CSA_hrms_file}{
CSA HRMS file
}
  \item{tempAlignedTableSubsetsFolder}{
tempAlignedTableSubsetsFolder
}
  \item{peaklist}{
IDSL.IPA peaklist
}
  \item{selectedIPApeaks}{
A vector of selected IDSL.IPA peaks only when a number of IDSL.IPA peaks from one peaklist is processed. When `NULL` is selected, the entire peaks in the peaklist are processed.
}
  \item{RTtolerance}{
retention time tolerance to detect common peaks
}
  \item{massError}{
Mass accuracy in Da
}
  \item{minSNRbaseline}{
A minimum baseline S/N threshold for IDSL.IPA pseudo-precursor m/z
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
  \item{minIonRangeDifference}{
Minimum distance (Da) between lowest and highest m/z to prevent clustering isotopic envelopes
}
  \item{minNumCSApeaks}{
Minumum number of ions in a CSA cluster
}
  \item{pearsonRHOthreshold}{
Minimum threshold for Pearson correlation rho values
}
  \item{outputCSAeic}{
When `NULL` CSA EICs are not plotted. `outputCSAeic` represents an address to save CSA EICs figures.
}
}
\value{
A dataframe peaklist of co-detected CSA analysis.
}
\references{
[1] Fakouri Baygi, S., Kumar, Y., Barupal, D.K. (2022). IDSL.IPA Characterizes the Organic Chemical Space in Untargeted LC/HRMS Data Sets. \emph{Journal of Proteome Research}, 21(6), 1485-1494, \doi{10.1021/acs.jproteome.2c00120}

[2] Fakouri Baygi, S., Fernando, S., Hopke, P.K., Holsen, T.M., Crimmins, B.S. (2021). Nontargeted discovery of novel contaminants in the Great Lakes region: A comparison of fish fillets and fish consumers. \emph{Environmental Science & Technology}, 55(6), 3765-3774, \doi{10.1021/acs.est.0c08507}
}