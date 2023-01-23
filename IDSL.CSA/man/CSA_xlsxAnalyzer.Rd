\name{CSA_xlsxAnalyzer}
\alias{CSA_xlsxAnalyzer}
\title{
CSA xlsx Analyzer
}
\description{
This function processes the spreadsheet of the CSA parameters to ensure the parameter inputs are consistent with the requirements of the IDSL.CSA pipeline.
}
\usage{
CSA_xlsxAnalyzer(spreadsheet)
}
\arguments{
  \item{spreadsheet}{
CSA tab of the parameter spreadsheet
}
}
\value{
This function returns the CSA parameters to feed the `CSA_workflow` function.
}
