IDSL.CSA_xlsxAnalyzer <- function(spreadsheet) {
  FSA_message("Initiated testing the IDSL.CSA workflow spreadsheet consistency!", failedMessage = FALSE)
  ##
  checkpoint_parameter <- FALSE
  ##
  if (length(spreadsheet) >= 4) {
    if (typeof(spreadsheet) == "list") {
      PARAM_Start <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
    } else {
      FSA_message("The IDSL.CSA workflow spreadsheet tab was not produced properly!")
      checkpoint_parameter <- FALSE
    }
  } else if (length(spreadsheet) == 1) {
    if (typeof(spreadsheet) == "character") {
      if (file.exists(spreadsheet)) {
        PARAM_Start <- readxl::read_xlsx(spreadsheet, sheet = "Start")
        PARAM_Start <- cbind(PARAM_Start[, 2], PARAM_Start[, 4])
        checkpoint_parameter <- TRUE
      } else {
        FSA_message("The IDSL.CSA workflow spreadsheet not found! It should be an Excel file with .xlsx extention!")
      }
    } else {
      FSA_message("The IDSL.CSA workflow spreadsheet was not produced properly!")
    }
  } else {
    FSA_message("The IDSL.CSA workflow spreadsheet was not produced properly!")
  }
  ##
  if (checkpoint_parameter) {
    ############################################################################
    x0001 <- which(PARAM_Start[, 1] == 'PARAM0001')
    if (length(x0001) == 0) {
      FSA_message("ERROR!!! Problem with PARAM0001!")
      checkpoint_parameter <- FALSE
    } else {
      PARAM0001 <- tolower(PARAM_Start[x0001, 2])
      if (PARAM0001 == "yes" | PARAM0001 == "no") {
        PARAM_Start[x0001, 2] <- PARAM0001
      } else {
        FSA_message("ERROR!!! Problem with PARAM0001!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0002 <- which(PARAM_Start[, 1] == 'PARAM0002')
    if (length(x0002) == 0) {
      FSA_message("ERROR!!! Problem with PARAM0002!")
      checkpoint_parameter <- FALSE
    } else {
      PARAM0002 <- tolower(PARAM_Start[x0002, 2])
      if (PARAM0002 == "yes" | PARAM0002 == "no") {
        PARAM_Start[x0002, 2] <- PARAM0002
      } else {
        FSA_message("ERROR!!! Problem with PARAM0002!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0003 <- which(PARAM_Start[, 1] == 'PARAM0003')
    if (length(x0003) == 0) {
      FSA_message("ERROR!!! Problem with PARAM0003!")
      checkpoint_parameter <- FALSE
    } else {
      PARAM0003 <- tolower(PARAM_Start[x0003, 2])
      if (PARAM0003 == "yes" | PARAM0003 == "no") {
        PARAM_Start[x0003, 2] <- PARAM0003
      } else {
        FSA_message("ERROR!!! Problem with PARAM0003!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    if (length(which(c(PARAM0001, PARAM0002, PARAM0003) == "yes")) > 1) {
      FSA_message("ERROR!!! PARAM0001, PARAM0002, PARAM0003 can not be 'YES' at the same time!")
      checkpoint_parameter <- FALSE
    }
    ##
    x0004 <- which(PARAM_Start[, 1] == 'PARAM0004')
    if (length(x0004) == 0) {
      FSA_message("ERROR!!! Problem with PARAM0004!")
      checkpoint_parameter <- FALSE
    } else {
      PARAM0004 <- tolower(PARAM_Start[x0004, 2])
      if (PARAM0004 == "yes" | PARAM0004 == "no") {
        PARAM_Start[x0004, 2] <- PARAM0004
      } else {
        FSA_message("ERROR!!! Problem with PARAM0004!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0005 <- which(PARAM_Start[, 1] == 'PARAM0005')
    if (length(x0005) == 0) {
      FSA_message("ERROR!!! Problem with PARAM0005!")
      checkpoint_parameter <- FALSE
    } else {
      PARAM0005 <- tolower(PARAM_Start[x0005, 2])
      if (PARAM0005 == "yes" | PARAM0005 == "no") {
        PARAM_Start[x0005, 2] <- PARAM0005
      } else {
        FSA_message("ERROR!!! Problem with PARAM0005!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0006 <- which(PARAM_Start[, 1] == 'PARAM0006')
    if (length(x0006) == 0) {
      FSA_message("ERROR!!! Problem with PARAM0006!")
      checkpoint_parameter <- FALSE
    } else {
      PARAM0006 <- tolower(PARAM_Start[x0006, 2])
      if (PARAM0006 == "yes" | PARAM0006 == "no") {
        PARAM_Start[x0006, 2] <- PARAM0006
      } else {
        FSA_message("ERROR!!! Problem with PARAM0006!")
        checkpoint_parameter <- FALSE
      }
    }
    ############################################################################
    if (!checkpoint_parameter) {
      FSA_message("ERROR!!! Problem with the `Start` spreadsheet tab!")
    }
  }
  ##############################################################################
  ##############################################################################
  ##############################################################################
  if (checkpoint_parameter) {
    ##
    if (PARAM0001 == "yes") {
      FSA_message("Initiated testing the `CSA` spreadsheet tab consistency!", failedMessage = FALSE)
      ##
      PARAM_CSA <- CSA_xlsxAnalyzer(spreadsheet)
      if (is.null(PARAM_CSA)) {
        FSA_message("ERROR!!! Problem with the `CSA` spreadsheet tab!")
        checkpoint_parameter <- FALSE
      }
    } else {
      PARAM_CSA <- NULL
    }
    ##
    if (PARAM0002 == "yes") {
      FSA_message("Initiated testing the `DDA` spreadsheet tab consistency!", failedMessage = FALSE)
      ##
      PARAM_DDA <- DDA_xlsxAnalyzer(spreadsheet)
      if (is.null(PARAM_DDA)) {
        FSA_message("ERROR!!! Problem with the `DDA` spreadsheet tab!")
        checkpoint_parameter <- FALSE
      }
    } else {
      PARAM_DDA <- NULL
    }
    ##
    if (PARAM0003 == "yes") {
      FSA_message("Initiated testing the `DIA` spreadsheet tab consistency!", failedMessage = FALSE)
      ##
      PARAM_DIA <- DIA_xlsxAnalyzer(spreadsheet)
      if (is.null(PARAM_DIA)) {
        FSA_message("ERROR!!! Problem with the `DIA` spreadsheet tab!")
        checkpoint_parameter <- FALSE
      }
    } else {
      PARAM_DIA <- NULL
    }
    ##
    if (PARAM0004 == "yes") {
      FSA_message("Initiated testing the `FSDB` spreadsheet tab consistency!", failedMessage = FALSE)
      ##
      PARAM_FSdb <- FSA_FSdb_xlsxAnalyzer(spreadsheet)
      if (is.null(PARAM_FSdb)) {
        FSA_message("ERROR!!! Problem with the `FSDB` spreadsheet tab!")
        checkpoint_parameter <- FALSE
      }
    } else {
      PARAM_FSdb <- NULL
    }
    ##
    if (PARAM0005 == "yes") {
      FSA_message("Initiated testing the `SpectraSimilarity` spreadsheet tab consistency!", failedMessage = FALSE)
      ##
      PARAM_SPEC <- FSA_SpectraSimilarity_xlsxAnalyzer(spreadsheet)
      if (is.null(PARAM_SPEC)) {
        FSA_message("ERROR!!! Problem with the `SpectraSimilarity` spreadsheet tab!")
        checkpoint_parameter <- FALSE
      }
    } else {
      PARAM_SPEC <- NULL
    }
    ##
    if (PARAM0006 == "yes") {
      FSA_message("Initiated testing the `AlignedTable` spreadsheet tab consistency!", failedMessage = FALSE)
      ##
      PARAM_AT <- CSA_AlignedTable_xlsxAnalyzer(spreadsheet)
      if (is.null(PARAM_AT)) {
        FSA_message("ERROR!!! Problem with the `AlignedTable` spreadsheet tab!")
        checkpoint_parameter <- FALSE
      }
    } else {
      PARAM_AT <- NULL
    }
    ############################################################################
    ######################### Sample import and export #########################
    ############################################################################
    x0007 <- which(PARAM_Start[, 1] == 'PARAM0007')
    if (length(x0007) == 0) {
      FSA_message("ERROR!!! Problem with PARAM0007!")
      checkpoint_parameter <- FALSE
    } else {
      ##
      if (PARAM0004 == "no") {
        if ((PARAM0005 == "yes")) {
          ##
          FSdb_file <- PARAM_Start[x0007, 2]
          FSdb_file <- gsub("\\", "/", FSdb_file, fixed = TRUE)
          PARAM_Start[x0007, 2] <- FSdb_file
          if (!file.exists(FSdb_file)) {
            FSA_message("ERROR!!! Problem with PARAM0007! Please ensure the full path is provided for the FSDB in .Rdata format OR select 'YES' for PARAM0005!")
            checkpoint_parameter <- FALSE
          }
        }
      }
    }
    ##
    if (PARAM0005 == "yes") {
      x0008 <- which(PARAM_Start[, 1] == 'PARAM0008')
      ##
      if (PARAM0001 == "no" & PARAM0002 == "no" & PARAM0003 == "no") {
        if (length(x0008) == 0) {
          FSA_message("ERROR!!! Problem with PARAM0008!")
          checkpoint_parameter <- FALSE
        } else {
          address_input_msp <- PARAM_Start[x0008, 2]
          address_input_msp <- gsub("\\", "/", address_input_msp, fixed = TRUE)
          PARAM_Start[x0008, 2] <- address_input_msp
          if (!dir.exists(address_input_msp)) {
            FSA_message("ERROR!!! Problem with PARAM0008! Please make sure the full path is provided!")
            checkpoint_parameter <- FALSE
          } else {
            file_name_sample_msp <- dir(path = address_input_msp)
            file_name_sample_msp <- file_name_sample_msp[grepl(".msp$", file_name_sample_msp, ignore.case = TRUE)]
            if (length(file_name_sample_msp) == 0) {
              FSA_message("ERROR!!! Problem with PARAM0008! No .msp file was detected in the designated folder!")
              checkpoint_parameter <- FALSE
            }
          }
        }
      } else {
        if (PARAM0001 == "yes") {
          if (!is.null(PARAM_CSA)) {
            ##
            refMSPcreationCheck <- file.exists(as.character(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0004'), 2]))
            if (refMSPcreationCheck) {
              annex_MSP <- "/CSA_REF_MSP"
            } else {
              annex_MSP <- "/CSA_MSP"
            }
            ##
            xCSA0011 <- which(PARAM_CSA[, 1] == 'CSA0011')
            address_input_msp <- PARAM_CSA[xCSA0011, 2]
            address_input_msp <- gsub("\\", "/", address_input_msp, fixed = TRUE)
            address_input_msp <- paste0(address_input_msp, annex_MSP)
            PARAM_Start[x0008, 2] <- address_input_msp
          }
          ##
        } else if (PARAM0002 == "yes") {
          if (!is.null(PARAM_DDA)) {
            ##
            refMSPcreationCheck <- file.exists(as.character(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0002'), 2]))
            if (refMSPcreationCheck) {
              annex_MSP <- "/DDA_REF_MSP"
            } else {
              annex_MSP <- "/DDA_MSP"
            }
            ##
            xDDA0011 <- which(PARAM_DDA[, 1] == 'DDA0011')
            address_input_msp <- PARAM_DDA[xDDA0011, 2]
            address_input_msp <- gsub("\\", "/", address_input_msp, fixed = TRUE)
            address_input_msp <- paste0(address_input_msp, annex_MSP)
            PARAM_Start[x0008, 2] <- address_input_msp
          }
          ##
        } else if (PARAM0003 == "yes") {
          if (!is.null(PARAM_DIA)) {
            ##
            refMSPcreationCheck <- file.exists(as.character(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0002'), 2]))
            if (refMSPcreationCheck) {
              annex_MSP <- "/DIA_REF_MSP"
            } else {
              annex_MSP <- "/DIA_MSP"
            }
            ##
            xDIA0011 <- which(PARAM_DIA[, 1] == 'DIA0011')
            address_input_msp <- PARAM_DIA[xDIA0011, 2]
            address_input_msp <- gsub("\\", "/", address_input_msp, fixed = TRUE)
            address_input_msp <- paste0(address_input_msp, annex_MSP)
            PARAM_Start[x0008, 2] <- address_input_msp
          }
        }
      }
    }
    ##
    if (PARAM0005 == "yes" | PARAM0006 == "yes") {
      x0009 <- which(PARAM_Start[, 1] == 'PARAM0009')
      if ((PARAM0001 == "no" & PARAM0002 == "no" & PARAM0003 == "no")) {
        if (length(x0009) == 0) {
          FSA_message("ERROR!!! Problem with PARAM0009!")
          checkpoint_parameter <- FALSE
        } else {
          output_sample <- PARAM_Start[x0009, 2]
          output_sample <- gsub("\\", "/", output_sample, fixed = TRUE)
          PARAM_Start[x0009, 2] <- output_sample
          if (!dir.exists(output_sample)) {
            tryCatch(dir.create(output_sample, recursive = TRUE), warning = function(w) {warning("Problem with PARAM0009! R cannot create the folder!")})
            if (!dir.exists(output_sample)) {
              checkpoint_parameter <- FALSE
            }
          }
        }
      } else {
        if (PARAM0001 == "yes") {
          output_path_sample <- PARAM_CSA[xCSA0011, 2]
          output_path_sample <- gsub("\\", "/", output_path_sample, fixed = TRUE)
          PARAM_Start[x0009, 2] <- output_path_sample
        } else if (PARAM0002 == "yes") {
          output_path_sample <- PARAM_DDA[xDDA0011, 2]
          output_path_sample <- gsub("\\", "/", output_path_sample, fixed = TRUE)
          PARAM_Start[x0009, 2] <- output_path_sample
        } else if (PARAM0003 == "yes") {
          output_path_sample <- PARAM_DIA[xDIA0011, 2]
          output_path_sample <- gsub("\\", "/", output_path_sample, fixed = TRUE)
          PARAM_Start[x0009, 2] <- output_path_sample
        }
      }
    }
    ############################################################################
    if ((PARAM0004 == "yes") & !is.null(PARAM_SPEC) & (PARAM0005 == "yes") & !is.null(PARAM_FSdb)) {
      ##
      FSdb0007 <- as.numeric(PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0007'), 2])
      SPEC0008 <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0008'), 2])
      if (FSdb0007 != SPEC0008) {
        FSA_message("Inconsistency between `FSdb0007` and `SPEC0008` parameters in the `FSDB` and `SpectraSimilarity` tabs!")
        checkpoint_parameter <- FALSE
      }
      ##
      FSdb0008 <- as.numeric(PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0008'), 2])
      SPEC0012 <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0012'), 2])
      if (FSdb0008 != SPEC0012) {
        FSA_message("Inconsistency between `FSdb0008` and `SPEC0012` parameters in the `FSDB` and `SpectraSimilarity` tabs!")
        checkpoint_parameter <- FALSE
      }
      ##
      FSdb0009 <- eval(parse(text = PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0009'), 2]))
      SPEC0013 <- eval(parse(text = PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0013'), 2]))
      if (FSdb0009 != SPEC0013) {
        FSA_message("Inconsistency between `FSdb0009` and `SPEC0013` parameters in the `FSDB` and `SpectraSimilarity` tabs!")
        checkpoint_parameter <- FALSE
      }
    }
  }
  ##############################################################################
  if (!checkpoint_parameter) {
    PARAM_total <- vector(mode = "list", 7)
    ##
  } else {
    ##
    PARAM_total <- list(PARAM_Start, PARAM_CSA, PARAM_DDA, PARAM_DIA, PARAM_FSdb, PARAM_SPEC, PARAM_AT)
    ##
    if (PARAM0001 == "yes") {
      FSA_message("The `CSA` spreadsheet tab is consistent with the IDSL.CSA workflow!", failedMessage = FALSE)
    }
    ##
    if (PARAM0002 == "yes") {
      FSA_message("The `DDA` spreadsheet tab is consistent with the IDSL.CSA workflow!", failedMessage = FALSE)
    }
    ##
    if (PARAM0003 == "yes") {
      FSA_message("The `DIA` spreadsheet tab is consistent with the IDSL.CSA workflow!", failedMessage = FALSE)
    }
    ##
    if (PARAM0004 == "yes") {
      FSA_message("The `FSDB` spreadsheet tab is consistent with the IDSL.CSA workflow!", failedMessage = FALSE)
    }
    ##
    if (PARAM0005 == "yes") {
      ##
      massErrorPrecursor <- PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0005'), 2]
      if (!is.na(massErrorPrecursor)) {
        FSA_message(paste0("NOTICE: Mass accuracy = '" , massErrorPrecursor, " Da' for precursor m/z (SPEC0005) was selected for spectra annotations! Empty annotation tables will be generated when .msp files do not contain 'PrecursorMZ' information!"), failedMessage = FALSE)
      } else {
        FSA_message("NOTICE: Mass accuracy for precursor m/z (SPEC0005) was not selected and precursor values will not be used for spectra annotations!", failedMessage = FALSE)
      }
      ##
      RTtolerance <- PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0006'), 2]
      if (!is.na(RTtolerance)) {
        FSA_message(paste0("NOTICE: Retention time tolerance = '" , RTtolerance, " min' (SPEC0006) was selected for spectra annotations! Empty annotation tables will be generated when .msp files do not contain 'Retention Time' information!"), failedMessage = FALSE)
      } else {
        FSA_message("NOTICE: Retention time tolerance (SPEC0006) was not selected and retention time values will not be used for spectra annotations!", failedMessage = FALSE)
      }
      ##
      FSA_message("The `SpectraSimilarity` spreadsheet tab is consistent with the IDSL.CSA workflow!", failedMessage = FALSE)
    }
    ##
    if (PARAM0006 == "yes") {
      FSA_message("The `AlignedTable` spreadsheet tab is consistent with the IDSL.CSA workflow!", failedMessage = FALSE)
    }
    ##
    FSA_message("The spreadsheet is consistent with the IDSL.CSA workflow!", failedMessage = FALSE)
  }
  ##
  names(PARAM_total) <- c("PARAM_Start", "PARAM_CSA", "PARAM_DDA", "PARAM_DIA", "PARAM_FSdb", "PARAM_SPEC", "PARAM_AT")
  ##############################################################################
  return(PARAM_total)
}