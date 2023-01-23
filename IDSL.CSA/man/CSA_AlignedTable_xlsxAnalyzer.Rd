\name{CSA_AlignedTable_xlsxAnalyzer}
\alias{CSA_AlignedTable_xlsxAnalyzer}
\title{
CSA Aligned Table xlsx Analyzer
}
\description{
This function processes the spreadsheet of the `AlignedTable` tab to ensure the parameter inputs are consistent with the requirements of the IDSL.CSA pipeline.
}
\usage{
CSA_AlignedTable_xlsxAnalyzer(spreadsheet)
}
\arguments{
  \item{spreadsheet}{
`AlignedTable` tab of the parameter spreadsheet
}
}
\value{
This function returns the aligned table parameters to feed the `aligned_fragmentation_spectra_annotator` function.
}