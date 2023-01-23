\name{CSA_workflow}
\alias{CSA_workflow}
\title{
CSA workflow
}
\description{
This function executes the CSA workflow.
}
\usage{
CSA_workflow(PARAM_CSA)
}
\arguments{
  \item{PARAM_CSA}{
PARAM_CSA
}
}
\value{
This module generates `.msp` files from DDA analysis.
}
\examples{
\donttest{
s_path <- system.file("extdata", package = "IDSL.CSA")
SSh1 <- paste0(s_path,"/CSA_parameters.xlsx")
## To see the results, use a known folder instead of the `tempdir()` command
temp_wd <- tempdir()
temp_wd_zip <- paste0(temp_wd, "/idsl_csa_test_files.zip")
spreadsheet <- readxl::read_xlsx(SSh1, sheet = "CSA")
PARAM_CSA <- cbind(spreadsheet[, 2], spreadsheet[, 4])
download.file(paste0("https://github.com/idslme/IDSL.CSA/blob/main/",
                     "CSA_educational_files/idsl_csa_test_files.zip?raw=true"),
              destfile = temp_wd_zip, mode = "wb")
unzip(temp_wd_zip, exdir = temp_wd)
PARAM_CSA[2, 2] <- "NO"
PARAM_CSA[3, 2] <- "NO"
PARAM_CSA[5, 2] <- temp_wd
PARAM_CSA[8, 2] <- temp_wd
PARAM_CSA[9, 2] <- "NA"
PARAM_CSA[11, 2] <- temp_wd
## To ensure `PARAM_CSA` is consistent with the `CSA_workflow`
PARAM_CSA <- CSA_xlsxAnalyzer(PARAM_CSA)
##
CSA_workflow(PARAM_CSA)
}
}