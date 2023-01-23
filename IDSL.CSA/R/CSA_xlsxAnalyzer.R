CSA_xlsxAnalyzer <- function(spreadsheet) {
  ##
  checkpoint_parameter <- FALSE
  ##
  if (typeof(spreadsheet) == "list") {
    if (ncol(spreadsheet) >= 4) {
      PARAM_CSA <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
      ##
    } else if (ncol(spreadsheet) == 2) {
      PARAM_CSA <- spreadsheet
      checkpoint_parameter <- TRUE
      ##
    } else {
      FSA_message("The `CSA` spreadsheet tab was not produced properly!")
    }
  } else if (typeof(spreadsheet) == "character") {
    if (length(spreadsheet) == 1) {
      if (file.exists(spreadsheet)) {
        PARAM_CSA <- readxl::read_xlsx(spreadsheet, sheet = "CSA")
        PARAM_CSA <- cbind(PARAM_CSA[, 2], PARAM_CSA[, 4])
        checkpoint_parameter <- TRUE
      } else {
        FSA_message("The `CSA` spreadsheet tab not found! It should be an Excel file with .xlsx extention!")
      }
    } else {
      FSA_message("The `CSA` spreadsheet tab was not produced properly!")
    }
  } else {
    FSA_message("The `CSA` spreadsheet tab was not produced properly!")
  }
  ##############################################################################
  if (checkpoint_parameter) {
    ##
    x0001 <- which(PARAM_CSA[, 1] == 'CSA0001')
    CSA0001 <- PARAM_CSA[x0001, 2]
    if (is.na(CSA0001)) {
      FSA_message("ERROR!!! Problem with CSA0001!")
      checkpoint_parameter <- FALSE
    } else {
      CSA0001 <- gsub(" ", "", tolower(CSA0001))
      if (CSA0001 == "yes" | CSA0001 == "no") {
        PARAM_CSA[x0001, 2] <- CSA0001
      } else {
        FSA_message("ERROR!!! Problem with CSA0001!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0002 <- which(PARAM_CSA[, 1] == 'CSA0002')
    CSA0002 <- PARAM_CSA[x0002, 2]
    if (is.na(CSA0002)) {
      FSA_message("ERROR!!! Problem with CSA0002!")
      checkpoint_parameter <- FALSE
    } else {
      CSA0002 <- gsub(" ", "", tolower(CSA0002))
      if (CSA0002 == "yes" | CSA0002 == "no") {
        PARAM_CSA[x0002, 2] <- CSA0002
      } else {
        FSA_message("ERROR!!! Problem with CSA0002!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0003 <- which(PARAM_CSA[, 1] == 'CSA0003')
    CSA0003 <- PARAM_CSA[x0003, 2]
    if (is.na(CSA0003)) {
      FSA_message("ERROR!!! Problem with CSA0003!")
      checkpoint_parameter <- FALSE
    } else {
      CSA0003 <- gsub(" ", "", tolower(CSA0003))
      if (CSA0003 == "yes" | CSA0003 == "no") {
        PARAM_CSA[x0003, 2] <- CSA0003
      } else {
        FSA_message("ERROR!!! Problem with CSA0003!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    number_processing_threads <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0004'), 2])
    if (length(number_processing_threads) == 0) {
      FSA_message("ERROR!!! Problem with CSA0004! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (number_processing_threads >= 1) {
        if ((number_processing_threads %% 1) != 0) {
          FSA_message("ERROR!!! Problem with CSA0004! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        FSA_message("ERROR!!! Problem with CSA0004! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    if (CSA0001 == "yes") {
      ##
      x0005 <- which(PARAM_CSA[, 1] == 'CSA0005')
      if (length(x0005) == 0) {
        FSA_message("ERROR!!! Problem with CSA0005!")
        checkpoint_parameter <- FALSE
      } else {
        input_path_hrms <- PARAM_CSA[x0005, 2]
        input_path_hrms <- gsub("\\", "/", input_path_hrms, fixed = TRUE)
        PARAM_CSA[x0005, 2] <- input_path_hrms
        if (!dir.exists(input_path_hrms)) {
          FSA_message("ERROR!!! Problem with CSA0005! Please make sure the full path is provided!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    ref_xlsx_file <- as.character(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0006'), 2])
    refMSPcreationCheck <- file.exists(ref_xlsx_file)
    ##
    if ((CSA0001 == "yes") & (tolower(ref_xlsx_file) != "na")) {
      refMSPcreationCheck <- TRUE
      ##
      listRefXlsxAnalyzer <- CSA_reference_xlsxAnalyzer(ref_xlsx_file, input_path_hrms, PARAM = PARAM_CSA, PARAM_ID = 'CSA0006', checkpoint_parameter)
      ref_table <- listRefXlsxAnalyzer[[1]]
      PARAM_CSA <- listRefXlsxAnalyzer[[2]]
      checkpoint_parameter <- listRefXlsxAnalyzer[[3]]
      listRefXlsxAnalyzer <- NULL
      samples_string <- paste0(unique(ref_table$Filename), collapse = ";")
      ref_table <- NULL
    }
    ##
    if (CSA0001 == "yes") {
      LHRMS <- 0
      x0007 <- which(PARAM_CSA[, 1] == 'CSA0007')
      if (is.na(PARAM_CSA[x0007, 2])) {
        FSA_message("ERROR!!! Problem with CSA0007!")
        checkpoint_parameter <- FALSE
      } else {
        ##
        if (refMSPcreationCheck) {
          PARAM_CSA[x0007, 2] <- samples_string
        }
        ##
        if (tolower(PARAM_CSA[x0007, 2]) == "all") {
          file_name_hrms <- dir(path = input_path_hrms)
          file_name_hrms <- file_name_hrms[grep(pattern = ".mzML$|.mzXML$|.CDF$", file_name_hrms, ignore.case = TRUE)]
          LHRMS <- length(file_name_hrms)
          if (LHRMS == 0) {
            FSA_message("ERROR!!! Problem with CSA0007! No mzML/mzXML/CDF file was detected in the folder!")
          }
        } else {
          samples_string <- PARAM_CSA[x0007, 2]
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
              FSA_message("ERROR!!! The following file(s) can not be detected in the reference xlsx file (CSA0006) (case sensitive even for file extensions):")
            } else {
              FSA_message("ERROR!!! Problem with CSA0007! not detected the following file(s) (case sensitive even for file extensions):")
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
      x0008 <- which(PARAM_CSA[, 1] == 'CSA0008')
      if (length(x0008) == 0) {
        FSA_message("ERROR!!! Problem with CSA0008!")
        checkpoint_parameter <- FALSE
      } else {
        inputPathPeaklist <- PARAM_CSA[x0008, 2]
        inputPathPeaklist <- gsub("\\", "/", inputPathPeaklist, fixed = TRUE)
        PARAM_CSA[x0008, 2] <- inputPathPeaklist
        if (!dir.exists(inputPathPeaklist)) {
          FSA_message("ERROR!!! Problem with CSA0008! Please make sure the full path is provided!")
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
      x0010 <- which(PARAM_CSA[, 1] == 'CSA0010')
      if (length(x0010) == 0) {
        FSA_message("ERROR!!! Problem with CSA0010!")
        checkpoint_parameter <- FALSE
      } else {
        CSA0010 <- tolower(PARAM_CSA[x0010, 2])
        if (CSA0010 == "yes" | CSA0010 == "no") {
          PARAM_CSA[x0010, 2] <- CSA0010
        } else {
          FSA_message("ERROR!!! Problem with CSA0010!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    x0011 <- which(PARAM_CSA[, 1] == 'CSA0011')
    if (length(x0011) == 0) {
      FSA_message("ERROR!!! Problem with CSA0011!")
      checkpoint_parameter <- FALSE
    } else {
      output_path <- gsub("\\", "/", PARAM_CSA[x0011, 2], fixed = TRUE)
      PARAM_CSA[x0011, 2] <- output_path
      if (!dir.exists(output_path)) {
        tryCatch(dir.create(output_path, recursive = TRUE), warning = function(w){warning("Problem with CSA0011! R cannot create the folder!")})
        if (!dir.exists(output_path)) {
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    ############################################################################
    ###################### Individual CSA .msp generation ######################
    ############################################################################
    if (CSA0001 == "yes") {
      x0012 <- which(PARAM_CSA[, 1] == 'CSA0012')
      if (length(x0012) == 0) {
        FSA_message("ERROR!!! Problem with CSA0012!")
        checkpoint_parameter <- FALSE
      } else {
        CSA0012 <- tolower(PARAM_CSA[x0012, 2])
        if (CSA0012 == "peaklist" | CSA0012 == "alignedtable") {
          PARAM_CSA[x0012, 2] <- CSA0012
        } else {
          FSA_message("ERROR!!! Problem with CSA0012!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      if (CSA0012 == "alignedtable") {
        listAlignmentFolderCheck <- IPA_peak_alignment_folder_xlsxAnalyzer(PARAM_CSA, PARAM_ID = 'CSA0009', checkpoint_parameter, correctedRTcheck = FALSE, CSAcheck = TRUE, allowedVerbose = TRUE)
        PARAM_CSA <- listAlignmentFolderCheck[[1]]
        checkpoint_parameter <- listAlignmentFolderCheck[[2]]
        listAlignmentFolderCheck <- NULL
      }
      ##
      RTtolerance <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0013'), 2])
      if (length(RTtolerance) == 0) {
        FSA_message("ERROR!!! Problem with CSA0013! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      } else {
        if (RTtolerance <= 0) {
          FSA_message("ERROR!!! Problem with CSA0013! This parameter should be a positive number!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      minSNRbaseline <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0014'), 2])
      if (length(minSNRbaseline) == 0) {
        FSA_message("ERROR!!! Problem with CSA0014! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      } else {
        if (minSNRbaseline <= 0) {
          FSA_message("ERROR!!! Problem with CSA0014! This parameter should be a positive number!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      smoothingWindowMS1 <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0015'), 2])
      if (length(smoothingWindowMS1) == 0) {
        FSA_message("ERROR!!! Problem with CSA0015! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      } else {
        if (smoothingWindowMS1 <= 0) {
          FSA_message("ERROR!!! Problem with CSA0015! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      massError <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0016'), 2])
      if (length(massError) == 0) {
        FSA_message("ERROR!!! Problem with CSA0016! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      } else {
        if (massError < 0) {
          FSA_message("ERROR!!! Problem with CSA0016! This parameter should be a positive number!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      scanTolerance <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0017'), 2])
      if (length(scanTolerance) == 0) {
        FSA_message("ERROR!!! Problem with CSA0017! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      } else {
        if (scanTolerance < 0) {
          FSA_message("ERROR!!! Problem with CSA0017! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      nSpline <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0018'), 2])
      if (length(nSpline) == 0) {
        FSA_message("ERROR!!! Problem with CSA0018! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      } else {
        if (nSpline < 10) {
          FSA_message("ERROR!!! Problem with CSA0018! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      topPercentagePeakHeight <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0019'), 2])
      if (length(topPercentagePeakHeight) == 0) {
        FSA_message("ERROR!!! Problem with CSA0019! This parameter should be a positive number between 50-100!")
        checkpoint_parameter <- FALSE
      } else {
        if (topPercentagePeakHeight < 50 | topPercentagePeakHeight > 100) {
          FSA_message("ERROR!!! Problem with CSA0019! This parameter should be a positive number between 50-100!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      minIonRangeDifference <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0020'), 2])
      if (length(minIonRangeDifference) == 0) {
        FSA_message("ERROR!!! Problem with CSA0020! This parameter should be a positive number greater than 0!")
        checkpoint_parameter <- FALSE
      } else {
        if (minIonRangeDifference < 0) {
          FSA_message("ERROR!!! Problem with CSA0020! This parameter should be a positive number greater than 0!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      minNumCSApeaks <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0021'), 2])
      if (length(minNumCSApeaks) == 0) {
        FSA_message("ERROR!!! Problem with CSA0021! This parameter should be a positive number >= 2!")
        checkpoint_parameter <- FALSE
      } else {
        minNumCSApeaks <- ceiling(minNumCSApeaks)
        if (minNumCSApeaks < 2) {
          FSA_message("ERROR!!! Problem with CSA0021! This parameter should be a positive number >= 2!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      pearsonRHOthreshold <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0022'), 2])
      if (length(pearsonRHOthreshold) == 0) {
        FSA_message("ERROR!!! Problem with CSA0022! This parameter should be a positive integer between 0.5-1!")
        checkpoint_parameter <- FALSE
      } else {
        if (pearsonRHOthreshold < 0.5 | pearsonRHOthreshold > 1) {
          FSA_message("ERROR!!! Problem with CSA0022! This parameter should be a positive integer between 0.5-1!")
          checkpoint_parameter <- FALSE
        }
      }
      ##########################################################################
      ### Unique tag aggregation by spectra similarity across entire samples ###
      ##########################################################################
      if (CSA0002 == "yes") {
        x0024 <- which(PARAM_CSA[, 1] == 'CSA0024')
        CSA0024 <- PARAM_CSA[x0024, 2]
        if (is.na(CSA0024)) {
          FSA_message("ERROR!!! Problem with CSA0024!")
          checkpoint_parameter <- FALSE
        } else {
          CSA0024 <- gsub(" ", "", tolower(CSA0024))
          if (CSA0024 == "yes" | CSA0024 == "no") {
            PARAM_CSA[x0024, 2] <- CSA0024
          } else {
            FSA_message("ERROR!!! Problem with CSA0024!")
            checkpoint_parameter <- FALSE
          }
        }
        ##
        if (refMSPcreationCheck) {
          massErrorRef <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0025'), 2])
          if (length(massErrorRef) == 0) {
            FSA_message("ERROR!!! Problem with CSA0025! This parameter should be a positive number!")
            checkpoint_parameter <- FALSE
          } else {
            if (massErrorRef < 0) {
              FSA_message("ERROR!!! Problem with CSA0025! This parameter should be a positive number!")
              checkpoint_parameter <- FALSE
            }
          }
        }
        ##
        RTtoleranceRef <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0026'), 2])
        if (length(RTtoleranceRef) == 0) {
          FSA_message("ERROR!!! Problem with CSA0026! This parameter should be a positive number!")
          checkpoint_parameter <- FALSE
        } else {
          if (RTtoleranceRef < 0) {
            FSA_message("ERROR!!! Problem with CSA0026! This parameter should be a positive number!")
            checkpoint_parameter <- FALSE
          }
        }
        ##
        if (!refMSPcreationCheck) {
          minCSAdetectionFrequency <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0027'), 2])
          if (length(minCSAdetectionFrequency) == 0) {
            FSA_message("ERROR!!! Problem with CSA0027! This parameter should be a positive number between 0 - 100!")
            checkpoint_parameter <- FALSE
          } else {
            if (!((minCSAdetectionFrequency > 0) & (minCSAdetectionFrequency < 100))) {
              FSA_message("ERROR!!! Problem with CSA0027! This parameter should be a positive number between 0 - 100!")
              checkpoint_parameter <- FALSE
            }
          }
        }
        ##
        x0028 <- which(PARAM_CSA[, 1] == 'CSA0028')
        allowedWeightedSpectralEntropy <- tolower(gsub(" ", "", PARAM_CSA[x0028, 2]))
        if (allowedWeightedSpectralEntropy == "1" | allowedWeightedSpectralEntropy == "t" | allowedWeightedSpectralEntropy == "true") {
          allowedWeightedSpectralEntropy <- TRUE
        } else {
          allowedWeightedSpectralEntropy <- FALSE
        }
        PARAM_CSA[x0028, 2] <- allowedWeightedSpectralEntropy
        ##
        minEntropySimilarity <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0029'), 2])
        if (length(minEntropySimilarity) == 0) {
          FSA_message("ERROR!!! Problem with CSA0029! This parameter should be a positive number between 0 - 1!")
          checkpoint_parameter <- FALSE
        } else {
          if (!((minEntropySimilarity >= 0) & (minEntropySimilarity <= 1))) {
            FSA_message("ERROR!!! Problem with CSA0029! This parameter should be a positive number between 0 - 1!")
            checkpoint_parameter <- FALSE
          }
        }
      }
    }
    ############################################################################
    ################### CSA aggregation on the aligned table ###################
    ############################################################################
    if (CSA0003 == "yes") {
      ##
      listAlignmentFolderCheck <- IPA_peak_alignment_folder_xlsxAnalyzer(PARAM_CSA, PARAM_ID = 'CSA0009', checkpoint_parameter, correctedRTcheck = FALSE, CSAcheck = TRUE, allowedVerbose = TRUE)
      PARAM_CSA <- listAlignmentFolderCheck[[1]]
      checkpoint_parameter <- listAlignmentFolderCheck[[2]]
      listAlignmentFolderCheck <- NULL
      ##
      RTtolerance_AT <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0030'), 2])
      if (length(RTtolerance_AT) == 0) {
        FSA_message("ERROR!!! Problem with CSA0030! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      } else {
        if (RTtolerance_AT <= 0) {
          FSA_message("ERROR!!! Problem with CSA0030! This parameter should be a positive number!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      minPercenetageDetection <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0031'), 2])
      if (length(minPercenetageDetection) == 0) {
        FSA_message("ERROR!!! Problem with CSA0031! This parameter should be a positive numberbetween 0-100!")
        checkpoint_parameter <- FALSE
      } else {
        if (minPercenetageDetection < 0 | minPercenetageDetection > 100) {
          FSA_message("ERROR!!! Problem with CSA0031! This parameter should be a positive number between 0-100!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      minNumberFragments <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0032'), 2])
      if (length(minNumberFragments) == 0) {
        FSA_message("ERROR!!! Problem with CSA0032! This parameter should be a positive integer >= 2 !")
        checkpoint_parameter <- FALSE
      } else {
        if (minNumberFragments >= 2) {
          if ((minNumberFragments %% 1) != 0) {
            FSA_message("ERROR!!! Problem with CSA0032! This parameter should be a positive integer >= 2 !")
            checkpoint_parameter <- FALSE
          }
        } else {
          FSA_message("ERROR!!! Problem with CSA0032! This parameter should be a positive integer >= 2 !")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      minTanimotoCoefficient1 <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0033'), 2])
      if (length(minTanimotoCoefficient1) == 0) {
        FSA_message("ERROR!!! Problem with CSA0033! This parameter should be a positive number between 0 - 1!")
        checkpoint_parameter <- FALSE
      } else {
        if (minTanimotoCoefficient1 < 0 | minTanimotoCoefficient1 > 1) {
          FSA_message("ERROR!!! Problem with CSA0033! This parameter should be a positive number between 0 - 1!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      minTanimotoCoefficient2 <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0034'), 2])
      if (length(minTanimotoCoefficient2) == 0) {
        FSA_message("ERROR!!! Problem with CSA0034! This parameter should be a positive number greater than `CSA0033` and less than 1!")
        checkpoint_parameter <- FALSE
      } else {
        if (minTanimotoCoefficient2 < minTanimotoCoefficient1 | minTanimotoCoefficient2 > 1) {
          FSA_message("ERROR!!! Problem with CSA0034! This parameter should be a positive number greater than `CSA0033` and less than 1!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0035 <- which(PARAM_CSA[, 1] == 'CSA0035')
      CSA0035 <- PARAM_CSA[x0035, 2]
      if (is.na(CSA0035)) {
        FSA_message("ERROR!!! Problem with CSA0035!")
        checkpoint_parameter <- FALSE
      } else {
        CSA0035 <- gsub(" ", "", tolower(CSA0035))
        if (CSA0035 == "yes" | CSA0035 == "no") {
          PARAM_CSA[x0035, 2] <- CSA0035
        } else {
          FSA_message("ERROR!!! Problem with CSA0035!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      massError <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0036'), 2])
      if (length(massError) == 0) {
        FSA_message("ERROR!!! Problem with CSA0036! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      } else {
        if (massError < 0) {
          FSA_message("ERROR!!! Problem with CSA0036! This parameter should be a positive number!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0037 <- which(PARAM_CSA[, 1] == 'CSA0037')
      allowedWeightedSpectralEntropy <- tolower(gsub(" ", "", PARAM_CSA[x0037, 2]))
      if (allowedWeightedSpectralEntropy == "1" | allowedWeightedSpectralEntropy == "t" | allowedWeightedSpectralEntropy == "true") {
        allowedWeightedSpectralEntropy <- TRUE
      } else {
        allowedWeightedSpectralEntropy <- FALSE
      }
      PARAM_CSA[x0037, 2] <- allowedWeightedSpectralEntropy
      ##
      minEntropySimilarity <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0038'), 2])
      if (length(minEntropySimilarity) == 0) {
        FSA_message("ERROR!!! Problem with CSA0038! This parameter should be a positive number between 0 - 1!")
        checkpoint_parameter <- FALSE
      } else {
        if (!((minEntropySimilarity >= 0) & (minEntropySimilarity <= 1))) {
          FSA_message("ERROR!!! Problem with CSA0038! This parameter should be a positive number between 0 - 1!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0039 <- which(PARAM_CSA[, 1] == 'CSA0039')
      if (length(x0039) == 0) {
        FSA_message("ERROR!!! Problem with CSA0039!")
        checkpoint_parameter <- FALSE
      } else {
        ##
        FSdb_file <- PARAM_CSA[x0039, 2]
        if (!(is.na(FSdb_file) | (tolower(FSdb_file) == "na"))) {
          ##
          FSdb_file <- gsub("\\", "/", FSdb_file, fixed = TRUE)
          PARAM_CSA[x0039, 2] <- FSdb_file
          if (!file.exists(FSdb_file)) {
            FSA_message("ERROR!!! Problem with CSA0039! Please ensure the full path is provided for the FSDB in .Rdata format OR select `NA`!")
            checkpoint_parameter <- FALSE
          }
        }
      }
    }
  }
  ##############################################################################
  if (!checkpoint_parameter) {
    PARAM_CSA <- NULL
  }
  ##
  return(PARAM_CSA)
}