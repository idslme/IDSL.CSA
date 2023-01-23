\name{DDA_workflow}
\alias{DDA_workflow}
\title{
DDA Workflow
}
\description{
This function runs the Data-Dependent Acquisition (DDA) analysis.
}
\usage{
DDA_workflow(PARAM_DDA)
}
\arguments{
  \item{PARAM_DDA}{
DDA parameters
}
}
\value{
This module generates `.msp` files from DDA analysis.
}
\examples{
s_path <- system.file("extdata", package = "IDSL.CSA")
SSh1 <- paste0(s_path,"/CSA_parameters.xlsx")
## To see the results, use a known folder instead of the `tempdir()` command
temp_wd <- tempdir()
temp_wd_zip <- paste0(temp_wd, "/idsl_dda_test_files.zip")
spreadsheet <- readxl::read_xlsx(SSh1, sheet = "DDA")
PARAM_DDA <- cbind(spreadsheet[, 2], spreadsheet[, 4])
download.file(paste0("https://github.com/idslme/IDSL.CSA/blob/main/",
                     "CSA_educational_files/idsl_dda_test_files.zip?raw=true"),
              destfile = temp_wd_zip, mode = "wb")
unzip(temp_wd_zip, exdir = temp_wd)
PARAM_DDA[2, 2] <- "no"
PARAM_DDA[4, 2] <- temp_wd
PARAM_DDA[7, 2] <- temp_wd
PARAM_DDA[8, 2] <- "NA"
PARAM_DDA[11, 2] <- temp_wd
## To ensure `PARAM_DDA` is consistent with the `DDA_workflow`
PARAM_DDA <- DDA_xlsxAnalyzer(PARAM_DDA)
##
DDA_workflow(PARAM_DDA)
}