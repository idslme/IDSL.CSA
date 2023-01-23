DIA_xlsxAnalyzer <- function(spreadsheet) {
  ##
  checkpoint_parameter <- FALSE
  ##
  if (typeof(spreadsheet) == "list") {
    if (ncol(spreadsheet) >= 4) {
      PARAM_DIA <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
      ##
    } else if (ncol(spreadsheet) == 2) {
      PARAM_DIA <- spreadsheet
      checkpoint_parameter <- TRUE
      ##
    } else {
      FSA_message("The DIA spreadsheet tab was not produced properly!")
    }
  } else if (typeof(spreadsheet) == "character") {
    if (length(spreadsheet) == 1) {
      if (file.exists(spreadsheet)) {
        PARAM_DIA <- readxl::read_xlsx(spreadsheet, sheet = "DIA")
        PARAM_DIA <- cbind(PARAM_DIA[, 2], PARAM_DIA[, 4])
        checkpoint_parameter <- TRUE
      } else {
        FSA_message("The DIA spreadsheet tab not found! It should be an Excel file with .xlsx extention!")
      }
    } else {
      FSA_message("The DIA spreadsheet tab was not produced properly!")
    }
  } else {
    FSA_message("The DIA spreadsheet tab was not produced properly!")
  }
  ##############################################################################
  if (checkpoint_parameter == TRUE) {
    ##
    x0001 <- which(PARAM_DIA[, 1] == 'DIA0001')
    DIA0001 <- PARAM_DIA[x0001, 2]
    if (is.na(DIA0001)) {
      FSA_message("ERROR!!! Problem with DIA0001!")
      checkpoint_parameter <- FALSE
    } else {
      DIA0001 <- gsub(" ", "", tolower(DIA0001))
      if (DIA0001 == "yes" | DIA0001 == "no") {
        PARAM_DIA[x0001, 2] <- DIA0001
      } else {
        FSA_message("ERROR!!! Problem with DIA0001!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0002 <- which(PARAM_DIA[, 1] == 'DIA0002')
    DIA0002 <- PARAM_DIA[x0002, 2]
    if (is.na(DIA0002)) {
      FSA_message("ERROR!!! Problem with DIA0002!")
      checkpoint_parameter <- FALSE
    } else {
      DIA0002 <- gsub(" ", "", tolower(DIA0002))
      if (DIA0002 == "yes" | DIA0002 == "no") {
        PARAM_DIA[x0002, 2] <- DIA0002
      } else {
        FSA_message("ERROR!!! Problem with DIA0002!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    number_processing_threads <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0003'), 2])
    if (length(number_processing_threads) == 0) {
      FSA_message("ERROR!!! Problem with DIA0003! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (number_processing_threads >= 1) {
        if ((number_processing_threads %% 1) != 0) {
          FSA_message("ERROR!!! Problem with DIA0003! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        FSA_message("ERROR!!! Problem with DIA0003! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    if (DIA0001 == "yes") {
      x0004 <- which(PARAM_DIA[, 1] == 'DIA0004')
      if (length(x0004) == 0) {
        FSA_message("ERROR!!! Problem with DIA0004!")
        checkpoint_parameter <- FALSE
      } else {
        input_path_hrms <- PARAM_DIA[x0004, 2]
        input_path_hrms <- gsub("\\", "/", input_path_hrms, fixed = TRUE)
        PARAM_DIA[x0004, 2] <- input_path_hrms
        if (!dir.exists(input_path_hrms)) {
          FSA_message("ERROR!!! Problem with DIA0004! Please make sure the full path is provided!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    refMSPcreationCheck <- FALSE
    ref_xlsx_file <- as.character(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0005'), 2])
    if (tolower(ref_xlsx_file) != "na") {
      refMSPcreationCheck <- TRUE
      ##
      listRefXlsxAnalyzer <- CSA_reference_xlsxAnalyzer(ref_xlsx_file, input_path_hrms, PARAM = PARAM_DIA, PARAM_ID = 'DIA0005', checkpoint_parameter)
      ref_table <- listRefXlsxAnalyzer[[1]]
      PARAM_DIA <- listRefXlsxAnalyzer[[2]]
      checkpoint_parameter <- listRefXlsxAnalyzer[[3]]
      listRefXlsxAnalyzer <- NULL
      samples_string <- paste0(unique(ref_table$Filename), collapse = ";")
      ref_table <- NULL
    }
    ##
    if (DIA0001 == "yes") {
      LHRMS <- 0
      x0006 <- which(PARAM_DIA[, 1] == 'DIA0006')
      if (is.na(PARAM_DIA[x0006, 2])) {
        FSA_message("ERROR!!! Problem with DIA0006!")
        checkpoint_parameter <- FALSE
      } else {
        ##
        if (refMSPcreationCheck) {
          PARAM_DIA[x0006, 2] <- samples_string
        }
        ##
        if (tolower(PARAM_DIA[x0006, 2]) == "all") {
          file_name_hrms <- dir(path = input_path_hrms)
          file_name_hrms <- file_name_hrms[grep(pattern = ".mzML$|.mzXML$|.CDF$", file_name_hrms, ignore.case = TRUE)]
          LHRMS <- length(file_name_hrms)
          if (LHRMS == 0) {
            FSA_message("ERROR!!! Problem with DIA0006! No mzML/mzXML/CDF file was detected in the folder!")
          }
        } else {
          samples_string <- PARAM_DIA[x0006, 2]
          file_name_hrms <- strsplit(samples_string, ";")[[1]]
          LHRMS <- length(file_name_hrms)
          ndHRMS <- do.call(c, lapply(file_name_hrms, function(i) {
            if (!file.exists(paste0(input_path_hrms, "/", i))) {
              i
            }
          }))
          ##
          if (!is.null(ndHRMS)) {
            if (refMSPcreationCheck) {
              FSA_message("ERROR!!! The following file(s) can not be detected in the reference xlsx file (DIA0005) (case sensitive even for file extensions):")
            } else {
              FSA_message("ERROR!!! Problem with DIA0006! not detected the following file(s) (case sensitive even for file extensions):")
            }
            ##
            LHRMS <- LHRMS - length(ndHRMS)
            for (i in ndHRMS) {
              FSA_message(i)
            }
            checkpoint_parameter <- FALSE
          }
        }
      }
      ##
      x0007 <- which(PARAM_DIA[, 1] == 'DIA0007')
      if (length(x0007) == 0) {
        FSA_message("ERROR!!! Problem with DIA0007!")
        checkpoint_parameter <- FALSE
      } else {
        inputPathPeaklist <- PARAM_DIA[x0007, 2]
        inputPathPeaklist <- gsub("\\", "/", inputPathPeaklist, fixed = TRUE)
        PARAM_DIA[x0007, 2] <- inputPathPeaklist
        if (!dir.exists(inputPathPeaklist)) {
          FSA_message("ERROR!!! Problem with DIA0007! Please make sure the full path is provided!")
          checkpoint_parameter <- FALSE
        } else {
          ######################################################################
          ## To see if the entire peaklists were generated for all HRMS files ##
          ######################################################################
          if (LHRMS > 0) {
            peaklistFileNames <- dir(path = inputPathPeaklist, pattern = ".Rdata$")
            peaklistFileNames <- peaklistFileNames[grep("^peaklist_", peaklistFileNames)]
            L_PL <- length(peaklistFileNames)
            ##
            if (LHRMS > L_PL) {
              checkpoint_parameter <- FALSE
              peaklistHRMSfileNames <- paste0("peaklist_", file_name_hrms, ".Rdata")
              ndPeaklists <- setdiff(peaklistHRMSfileNames, peaklistFileNames)
              ndPeaklists <- gsub("^peaklist_|.Rdata$", "", ndPeaklists)
              FSA_message("Error!!! peaklist files are not available for the following HRMS file(s):")
              for (i in ndPeaklists) {
                FSA_message(i)
              }
            }
          }
        }
      }
      ##
      listAlignmentFolderCheck <- IPA_peak_alignment_folder_xlsxAnalyzer(PARAM_DIA, PARAM_ID = 'DIA0008', checkpoint_parameter, correctedRTcheck = FALSE, CSAcheck = FALSE, allowedVerbose = TRUE)
      PARAM_DIA <- listAlignmentFolderCheck[[1]]
      checkpoint_parameter <- listAlignmentFolderCheck[[2]]
      listAlignmentFolderCheck <- NULL
      ##
      if (LHRMS == 1) {
        DIA0009 <- PARAM_DIA[which(PARAM_DIA[, 1] == "DIA0009"), 2]
        if (is.na(DIA0009)) {
          checkpoint_parameter <- FALSE
          FSA_message("ERROR!!! Problem with DIA0009! This parameter should be 'All' or a vector of indices!")
        } else if (gsub(" ", "", tolower(DIA0009)) == "all") {
          FSA_message("The enitre 12C m/z values in the peaklist were placed in the processing row!", failedMessage = TRUE)
        } else {
          peaklist <- IDSL.IPA::loadRdata(paste0(inputPathPeaklist, "/peaklist_", file_name_hrms, ".Rdata"))
          n_peaks <- dim(peaklist)[1]
          ##
          selectedIPApeaks <- tryCatch(eval(parse(text = paste0("c(", DIA0009, ")"))), error = function(e){NULL})
          if (is.null(selectedIPApeaks) | (max(selectedIPApeaks) > n_peaks)) {
            checkpoint_parameter <- FALSE
            FSA_message("ERROR!!! Problem with DIA0009! The range of indices are out of the peaklist dimension!")
          } else {
            FSA_message("The following peak IDs were selected for processing: ")
            for (id in 1:length(selectedIPApeaks)) {
              FSA_message(paste0(selectedIPApeaks[id], " - ", peaklist[selectedIPApeaks[id], 3],  " - ", peaklist[selectedIPApeaks[id], 8]))
            }
          }
        }
      }
      ##
      x0010 <- which(PARAM_DIA[, 1] == 'DIA0010')
      if (length(x0010) == 0) {
        FSA_message("ERROR!!! Problem with DIA0010!")
        checkpoint_parameter <- FALSE
      } else {
        DIA0010 <- tolower(PARAM_DIA[x0010, 2])
        if (DIA0010 == "yes" | DIA0010 == "no") {
          PARAM_DIA[x0010, 2] <- DIA0010
        } else {
          FSA_message("ERROR!!! Problem with DIA0010!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    x0011 <- which(PARAM_DIA[, 1] == 'DIA0011')
    if (length(x0011) == 0) {
      FSA_message("ERROR!!! Problem with DIA0011!")
      checkpoint_parameter <- FALSE
    } else {
      output_path <- gsub("\\", "/", PARAM_DIA[x0011, 2], fixed = TRUE)
      PARAM_DIA[x0011, 2] <- output_path
      if (!dir.exists(output_path)) {
        tryCatch(dir.create(output_path, recursive = TRUE), warning = function(w){warning("Problem with DIA0011! R cannot create the folder!")})
        if (!dir.exists(output_path)) {
          checkpoint_parameter <- FALSE
        }
      }
    }
    ############################################################################
    ########################### DIA .msp generation ############################
    ############################################################################
    if (number_processing_threads > 1) {
      x0012 <- which(PARAM_DIA[, 1] == 'DIA0012')
      parallelizationMode <- PARAM_DIA[x0012, 2]
      if (is.na(parallelizationMode)) {
        FSA_message("ERROR!!! Problem with DIA0012!")
        checkpoint_parameter <- FALSE
      } else {
        parallelizationMode <- gsub(" ", "", tolower(parallelizationMode))
        if (parallelizationMode == "samplemode" | parallelizationMode == "peakmode") {
          PARAM_DIA[x0012, 2] <- parallelizationMode
        } else {
          FSA_message("ERROR!!! Problem with DIA0012!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##############################################################################
    PARAM_DIA <- rbind(PARAM_DIA, c('DIA0013', 2))
    ############################################################################
    msLevel <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0013'), 2])
    if (length(msLevel) == 0) {
      FSA_message("ERROR!!! Problem with DIA0013! This parameter should be 1 or 2 !")
      checkpoint_parameter <- FALSE
    } else {
      if (!((msLevel == 1) | (msLevel == 2))) {
        FSA_message("ERROR!!! Problem with DIA0013! This parameter should be 1 or 2 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    intensityThresholdFragment <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0014'), 2])
    if (length(intensityThresholdFragment) == 0) {
      FSA_message("ERROR!!! Problem with DIA0014! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (intensityThresholdFragment < 0) {
        FSA_message("ERROR!!! Problem with DIA0014! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    smoothingWindowMS1 <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0015'), 2])
    if (length(smoothingWindowMS1) == 0) {
      FSA_message("ERROR!!! Problem with DIA0015! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (smoothingWindowMS1 <= 0) {
        FSA_message("ERROR!!! Problem with DIA0015! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    if (msLevel == 2) {
      smoothingWindowMS2 <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0016'), 2])
      if (length(smoothingWindowMS2) == 0) {
        FSA_message("ERROR!!! Problem with DIA0016! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      } else {
        if (smoothingWindowMS2 <= 0) {
          FSA_message("ERROR!!! Problem with DIA0016! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    massError <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0017'), 2])
    if (length(massError) == 0) {
      FSA_message("ERROR!!! Problem with DIA0017! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (massError < 0) {
        FSA_message("ERROR!!! Problem with DIA0017! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    scanTolerance <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0018'), 2])
    if (length(scanTolerance) == 0) {
      FSA_message("ERROR!!! Problem with DIA0018! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (scanTolerance < 0) {
        FSA_message("ERROR!!! Problem with DIA0018! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    nSpline <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0019'), 2])
    if (length(nSpline) == 0) {
      FSA_message("ERROR!!! Problem with DIA0019! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (nSpline < 10) {
        FSA_message("ERROR!!! Problem with DIA0019! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    topPercentagePeakHeight <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0020'), 2])
    if (length(topPercentagePeakHeight) == 0) {
      FSA_message("ERROR!!! Problem with DIA0020! This parameter should be a positive numberbetween 50-100!")
      checkpoint_parameter <- FALSE
    } else {
      if (topPercentagePeakHeight < 50 | topPercentagePeakHeight > 100) {
        FSA_message("ERROR!!! Problem with DIA0020! This parameter should be a positive number between 50-100!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    pearsonRHOthreshold <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0021'), 2])
    if (length(pearsonRHOthreshold) == 0) {
      FSA_message("ERROR!!! Problem with DIA0021! This parameter should be a positive integer between 0.5-1!")
      checkpoint_parameter <- FALSE
    } else {
      if (pearsonRHOthreshold < 0.5 | pearsonRHOthreshold > 1) {
        FSA_message("ERROR!!! Problem with DIA0021! This parameter should be a positive integer between 0.5-1!")
        checkpoint_parameter <- FALSE
      }
    }
    ############################################################################
    #### Unique tag aggregation by spectra similarity across entire samples ####
    ############################################################################
    if (DIA0002 == "yes") {
      x0023 <- which(PARAM_DIA[, 1] == 'DIA0023')
      DIA0023 <- PARAM_DIA[x0023, 2]
      if (is.na(DIA0023)) {
        FSA_message("ERROR!!! Problem with DIA0023!")
        checkpoint_parameter <- FALSE
      } else {
        DIA0023 <- gsub(" ", "", tolower(DIA0023))
        if (DIA0023 == "yes" | DIA0023 == "no") {
          PARAM_DIA[x0023, 2] <- DIA0023
        } else {
          FSA_message("ERROR!!! Problem with DIA0023!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      if (refMSPcreationCheck) {
        massErrorRef <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0024'), 2])
        if (length(massErrorRef) == 0) {
          FSA_message("ERROR!!! Problem with DIA0024! This parameter should be a positive number!")
          checkpoint_parameter <- FALSE
        } else {
          if (massErrorRef < 0) {
            FSA_message("ERROR!!! Problem with DIA0024! This parameter should be a positive number!")
            checkpoint_parameter <- FALSE
          }
        }
      }
      ##
      RTtoleranceRef <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0025'), 2])
      if (length(RTtoleranceRef) == 0) {
        FSA_message("ERROR!!! Problem with DIA0025! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      } else {
        if (RTtoleranceRef < 0) {
          FSA_message("ERROR!!! Problem with DIA0025! This parameter should be a positive number!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      if (!refMSPcreationCheck) {
        minDIAdetectionFrequency <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0026'), 2])
        if (length(minDIAdetectionFrequency) == 0) {
          FSA_message("ERROR!!! Problem with DIA0026! This parameter should be a positive number between 0 - 100!")
          checkpoint_parameter <- FALSE
        } else {
          if (!((minDIAdetectionFrequency > 0) & (minDIAdetectionFrequency < 100))) {
            FSA_message("ERROR!!! Problem with DIA0026! This parameter should be a positive number between 0 - 100!")
            checkpoint_parameter <- FALSE
          }
        }
      }
      ##
      x0027 <- which(PARAM_DIA[, 1] == 'DIA0027')
      allowedWeightedSpectralEntropy <- tolower(gsub(" ", "", PARAM_DIA[x0027, 2]))
      if (allowedWeightedSpectralEntropy == "1" | allowedWeightedSpectralEntropy == "t" | allowedWeightedSpectralEntropy == "true") {
        allowedWeightedSpectralEntropy <- TRUE
      } else {
        allowedWeightedSpectralEntropy <- FALSE
      }
      PARAM_DIA[x0027, 2] <- allowedWeightedSpectralEntropy
      ##
      minEntropySimilarity <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0028'), 2])
      if (length(minEntropySimilarity) == 0) {
        FSA_message("ERROR!!! Problem with DIA0028! This parameter should be a positive number between 0 - 1!")
        checkpoint_parameter <- FALSE
      } else {
        if (!((minEntropySimilarity >= 0) & (minEntropySimilarity <= 1))) {
          FSA_message("ERROR!!! Problem with DIA0028! This parameter should be a positive number between 0 - 1!")
          checkpoint_parameter <- FALSE
        }
      }
    }
  }
  ##############################################################################  
  if (!checkpoint_parameter) {
    PARAM_DIA <- NULL
  }
  ##
  return(PARAM_DIA)
}