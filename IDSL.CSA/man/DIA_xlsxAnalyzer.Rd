\name{DIA_xlsxAnalyzer}
\alias{DIA_xlsxAnalyzer}
\title{
DIA xlsx Analyzer for DIA analysis
}
\description{
This function processes the spreadsheet of the DIA spreadsheet tab to ensure the parameter inputs are in agreement with requirements of the Data-Independent Acquisition (DIA) analysis.
}
\usage{
DIA_xlsxAnalyzer(spreadsheet)
}
\arguments{
  \item{spreadsheet}{
DIA spreadsheet tab
}
}
\value{
DIA parameters to feed the `DIA_workflow` function.
}