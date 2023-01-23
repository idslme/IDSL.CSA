\name{DDA2msp}
\alias{DDA2msp}
\title{
DDA to msp
}
\description{
DDA to msp
}
\usage{
DDA2msp(input_path_hrms, file_name_hrms = NULL, number_processing_threads = 1)
}
\arguments{
  \item{input_path_hrms}{
path to the HRMS file
}
  \item{file_name_hrms}{
file_name_hrms
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
This module generates `.msp` files from DDA analysis.
}
\examples{
## To see the results, use a known folder instead of the `tempdir()` command
temp_wd <- tempdir()
temp_wd_zip <- paste0(temp_wd, "/idsl_rawdda_test_files.zip")
download.file(paste0("https://github.com/idslme/IDSL.CSA/blob/main/",
                     "CSA_educational_files/idsl_rawdda_test_files.zip?raw=true"),
              destfile = temp_wd_zip, mode = "wb")
unzip(temp_wd_zip, exdir = temp_wd)
DDA2msp(input_path_hrms = temp_wd, file_name_hrms = NULL, number_processing_threads = 1)
}