\name{DDA_xlsxAnalyzer}
\alias{DDA_xlsxAnalyzer}
\title{
xlsx Analyzer for DDA analysis
}
\description{
This function processes the spreadsheet of the DDA spreadsheet tab to ensure the parameter inputs are in agreement with requirements of the Data-Dependent Acquisition (DDA) analysis.
}
\usage{
DDA_xlsxAnalyzer(spreadsheet)
}
\arguments{
  \item{spreadsheet}{
DDA spreadsheet tab
}
}
\value{
DDA parameters to feed the `DDA_workflow` function.
}