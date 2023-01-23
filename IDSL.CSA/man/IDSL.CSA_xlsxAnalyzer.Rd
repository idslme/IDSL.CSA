\name{IDSL.CSA_xlsxAnalyzer}
\alias{IDSL.CSA_xlsxAnalyzer}
\title{
IDSL.CSA workflow xlsx Analyzer
}
\description{
This function processes the spreadsheet of the CSA parameters to ensure the parameter inputs are consistent with the requirements of the IDSL.CSA pipeline.
}
\usage{
IDSL.CSA_xlsxAnalyzer(spreadsheet)
}
\arguments{
  \item{spreadsheet}{
'Start' tab of the parameter spreadsheet
}
}
\value{
This function returns the CSA parameters to feed the `IDSL.CSA_workflow` function.
}
