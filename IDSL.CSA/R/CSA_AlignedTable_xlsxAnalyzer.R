CSA_AlignedTable_xlsxAnalyzer <- function(spreadsheet) {
  ##
  checkpoint_parameter <- FALSE
  #
  if (length(spreadsheet) >= 4) {
    if (typeof(spreadsheet) == "list") {
      PARAM_AT <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
    } else {
      FSA_message("The `AlignedTable` spreadsheet tab was not produced properly!")
    }
  } else if (length(spreadsheet) == 1) {
    if (typeof(spreadsheet) == "character") {
      if (file.exists(spreadsheet)) {
        PARAM_AT <- readxl::read_xlsx(spreadsheet, sheet = "AlignedTable")
        PARAM_AT <- cbind(PARAM_AT[, 2], PARAM_AT[, 4])
        checkpoint_parameter <- TRUE
      } else {
        FSA_message("The `AlignedTable` spreadsheet tab not found! It should be an Excel file with .xlsx extention!")
      }
    } else {
      FSA_message("The `AlignedTable` spreadsheet tab was not produced properly!")
    }
  } else {
    FSA_message("The `AlignedTable` spreadsheet tab was not produced properly!")
  }
  ##
  if (checkpoint_parameter) {
    ############################################################################
    listAlignmentFolderCheck <- IPA_peak_alignment_folder_xlsxAnalyzer(PARAM_AT, PARAM_ID = 'AT0001', checkpoint_parameter, correctedRTcheck = FALSE, CSAcheck = FALSE)
    PARAM_AT <- listAlignmentFolderCheck[[1]]
    checkpoint_parameter <- listAlignmentFolderCheck[[2]]
    listAlignmentFolderCheck <- NULL
    ##
    number_processing_threads <- as.numeric(PARAM_AT[which(PARAM_AT[, 1] == 'AT0002'), 2])
    if (length(number_processing_threads) == 0) {
      FSA_message("ERROR!!! Problem with AT0002! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (number_processing_threads >= 1) {
        if ((number_processing_threads %% 1) != 0) {
          FSA_message("ERROR!!! Problem with AT0002! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        FSA_message("ERROR!!! Problem with AT0002! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0003 <- as.numeric(PARAM_AT[which(PARAM_AT[, 1] == 'AT0003'), 2])
    if (length(x0003) == 0) {
      FSA_message("ERROR!!! Problem with AT0003! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (x0003 >= 1) {
        if ((x0003 %% 1) != 0) {
          FSA_message("ERROR!!! Problem with AT0003! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        FSA_message("ERROR!!! Problem with AT0003! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0004 <- as.numeric(PARAM_AT[which(PARAM_AT[, 1] == 'AT0004'), 2])
    if (length(x0004) == 0) {
      FSA_message("ERROR!!! Problem with AT0004! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (x0004 >= 1) {
        if ((x0004 %% 1) != 0) {
          FSA_message("ERROR!!! Problem with AT0004! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        FSA_message("ERROR!!! Problem with AT0004! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0005 <- which(PARAM_AT[, 1] == 'AT0005')
    adjustFreqRankCheck <- tolower(gsub(" ", "", PARAM_AT[x0005, 2]))
    if (adjustFreqRankCheck == "1" | adjustFreqRankCheck == "t" | adjustFreqRankCheck == "true") {
      adjustFreqRankCheck <- TRUE
    } else {
      adjustFreqRankCheck <- FALSE
    }
    PARAM_AT[x0005, 2] <- adjustFreqRankCheck
    ##
    RTtolerance <- as.numeric(PARAM_AT[which(PARAM_AT[, 1] == 'AT0006'), 2])
    if (length(RTtolerance) == 0) {
      FSA_message("ERROR!!! Problem with AT0006! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (RTtolerance < 0) {
        FSA_message("ERROR!!! Problem with AT0006! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
  }
  ##############################################################################
  if (!checkpoint_parameter) {
    PARAM_AT <- NULL
  }
  ##
  return(PARAM_AT)
}