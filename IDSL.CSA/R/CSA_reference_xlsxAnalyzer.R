CSA_reference_xlsxAnalyzer <- function(ref_xlsx_file, input_path_hrms = NULL, PARAM = NULL, PARAM_ID = "", checkpoint_parameter = TRUE) {
  ##
  ref_xlsx_file <- gsub("\\", "/", ref_xlsx_file, fixed = TRUE)
  ##
  if (file.exists(ref_xlsx_file)) {
    ##
    ############################################################################
    ##
    txtRef2csvREF <- function(txtFileAddress, input_path_hrms, PARAM_ID) {
      ##
      txtFile <- read.delim(txtFileAddress, header = TRUE, stringsAsFactors = FALSE)
      nROWtxtFile <- nrow(txtFile)
      ##
      file_name_hrms <- dir(path = input_path_hrms , pattern = ".mzML$|.mzXML$|.CDF$", ignore.case = TRUE)
      ##
      ref_csv_file <- do.call(rbind, lapply(file_name_hrms, function(f) {
        cbind(rep(f, nROWtxtFile), txtFile)
      }))
      ref_csv_file <- data.frame(ref_csv_file, stringsAsFactors = FALSE)
      colnames(ref_csv_file)[1] <- "Filename"
      ##
      csvFileAddress <- gsub(".txt$", ".csv", txtFileAddress, ignore.case = TRUE)
      write.csv(ref_csv_file, file = csvFileAddress, row.names = FALSE)
      ##
      FSA_message(paste0("The reference file represented by `", PARAM_ID, "` was converted into the `.csv` format!"), failedMessage = FALSE)
      ##
      return(csvFileAddress)
    }
    ##
    ############################################################################
    ##
    strSplitRefFileName <- strsplit(ref_xlsx_file, "[.]")[[1]]
    RefFileFormat <- tolower(strSplitRefFileName[length(strSplitRefFileName)])
    if (RefFileFormat == "txt") {
      ref_xlsx_file <- txtRef2csvREF(txtFileAddress = ref_xlsx_file, input_path_hrms, PARAM_ID)
      RefFileFormat <- "csv"
    }
    ##
    if (!is.null(PARAM)) {
      xPARAM_ID <- which(PARAM[, 1] == PARAM_ID)
      PARAM[xPARAM_ID, 2] <- ref_xlsx_file
    }
    ##
    if (RefFileFormat == "xlsx") {
      ref_table <- data.frame(readxl::read_xlsx(path = ref_xlsx_file))
    } else if (RefFileFormat == "csv") {
      ref_table <- data.frame(read.csv(file = ref_xlsx_file, header = TRUE))
    } else {
      stop("Incorrect reference file format")
    }
    ##
    col <- colnames(ref_table)
    x_fn <- which(col == 'Filename')
    x_name <- which(col == 'Name')
    x_mz <- which(col == 'PrecursorMZ')
    x_RT <- which(col == 'Precursor_RT')
    if (!(length(x_fn) > 0 & length(x_name) > 0 & length(x_mz) > 0 & length(x_RT) > 0)) {
      FSA_message(paste0("ERROR!!! Problem with `", PARAM_ID, "`! Incorrect column headers in the reference spreadsheet -> The following columns should be detected in the spreadsheet : 'Filename', 'Name', 'PrecursorMZ', 'Precursor_RT' - case sensitive"))
      checkpoint_parameter <- FALSE
    }
  } else {
    FSA_message(paste0("ERROR!!! Problem with `", PARAM_ID, "`! File was not detected!"))
    ##
    ref_table <- NULL
    checkpoint_parameter <- FALSE
  }
  ##
  listRefXlsxAnalyzer <- list(ref_table, PARAM, checkpoint_parameter)
  ##
  return(listRefXlsxAnalyzer)
}