DDA_workflow <- function(PARAM_DDA) {
  ##
  msLevelDDA <- 2
  ##
  x0000 <- which(PARAM_DDA[, 1] == 'DDA0000')
  if (length(x0000) > 0) {
    DDA_Analysis_Mode <- tolower(PARAM_DDA[x0000, 2])
  } else {
    DDA_Analysis_Mode <- "idsl.ipa"
  }
  ##
  DDA0001 <- tolower(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0001'), 2])
  DDA0002 <- tolower(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0002'), 2])
  NPT <- as.numeric(PARAM_DDA[which(PARAM_DDA[, 1] == "DDA0003"), 2])
  input_path_hrms <- PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0004'), 2]
  ##
  output_address <- PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0011'), 2]
  FSA_dir.create(output_address, allowedUnlink = FALSE)
  opendir(output_address)
  ##
  ##############################################################################
  ## To create log record for IDSL.CSA
  initiation_time <- Sys.time()
  timeZone <- tryCatch(Sys.timezone(), warning = function(w) {"UTC"}, error = function(e) {"UTC"})
  .logFSA <- NULL
  .logFSA <<- paste0(output_address, "/logDDA_performance.txt")
  FSA_logRecorder(paste0(rep("", 100), collapse = "="))
  FSA_logRecorder(paste0("mzML/mzXML/netCDF: ", input_path_hrms))
  FSA_logRecorder(paste0("OUTPUT: ", output_address))
  FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
  FSA_logRecorder("Initiated DDA workflow!")
  FSA_logRecorder(paste0(as.character(initiation_time), " ", timeZone))
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder(paste0(PARAM_DDA[, 1], "\t", PARAM_DDA[, 2]),  allowedPrinting = FALSE)
  FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  ##############################################################################
  ##
  ref_xlsx_file <- as.character(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0005'), 2])
  refMSPcreationCheck <- file.exists(ref_xlsx_file)
  ##
  ##############################################################################
  ##
  if (refMSPcreationCheck) {
    ##
    refDDAtable <- CSA_reference_xlsxAnalyzer(ref_xlsx_file, input_path_hrms)[[1]]
    ##
    refHRMSindexList <- FSA_R.aggregate(refDDAtable$`Filename`)
    file_name_hrms <- names(refHRMSindexList)
    mzRef <- as.numeric(refDDAtable$`PrecursorMZ`)
    RTref <- as.numeric(refDDAtable$`Precursor_RT`)
    refDDAtable$`PrecursorMZ` <- NULL
    refDDAtable$`Precursor_RT` <- NULL
    ##
    massErrorRef <- as.numeric(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0021'), 2])
    RTtoleranceRef <- as.numeric(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0022'), 2])
    refDDAtable <- data.frame(refDDAtable)
    ##
    output_DDA_MSP <- paste0(output_address, "/DDA_REF_MSP")
    tryCatch(FSA_dir.create(output_DDA_MSP, allowedUnlink = TRUE), error = function(e) {FSA_message("ERROR!!! The .msp files inside the `DDA_REF_MSP` folder may influence the final output!")})
    ##
    str_xlsx_file <- strsplit(ref_xlsx_file, "/")[[1]]
    mspFileName <- str_xlsx_file[length(str_xlsx_file)]
    mspFileName <- gsub("[.]xlsx$|[.]csv$", ".msp", mspFileName, ignore.case = TRUE)
    mspFileName <- paste0("DDA_REF_MSP_", mspFileName)
    ##
  } else {
    refHRMSindexList <- NULL
    refDDAtable <- NULL
    massErrorRef <- 0
    RTtoleranceRef <- 0
    ##
    samples_string <- PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0006'), 2]
    if (tolower(samples_string) == "all") {
      file_name_hrms <- dir(path = input_path_hrms)
      file_name_hrms <- file_name_hrms[grep(pattern = ".mzML$|.mzXML$|.CDF$", file_name_hrms, ignore.case = TRUE)]
    } else {
      file_name_hrms <- strsplit(samples_string, ";")[[1]]
    }
    ##
    output_DDA_MSP <- paste0(output_address, "/DDA_MSP")
    FSA_dir.create(output_DDA_MSP, allowedUnlink = FALSE)
  }
  ##
  ##############################################################################
  ##
  IDSL.IPA::opendir(output_DDA_MSP)
  ##
  LHRMS <- length(file_name_hrms)
  if (LHRMS == 0) {
    stop(FSA_logRecorder("EMPTY HRMS FOLDER!!!"))
  }
  ##
  ##############################################################################
  ##
  if (DDA0001 == "yes") {
    parallelizationMode <- tolower(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0012'), 2])
    ##
    if (DDA_Analysis_Mode == "rawddaspectra") {
      ##
      FSA_logRecorder("Initiated Data Dependent Acquisition (DDA) analysis using raw spectra!")
      ##
      if (refMSPcreationCheck) {
        rawDDAspectraVar <- vector(mode = "list", 4)
        names(rawDDAspectraVar) <- c("precursorMZvec", "precursorRTvec", "massError", "RTtolerance")
        rawDDAspectraVar[["massError"]] <- massErrorRef
        rawDDAspectraVar[["RTtolerance"]] <- RTtoleranceRef
      } else {
        rawDDAspectraVar <- NULL
        inputPathPeaklist <- NULL
        indexedIPApeaksCheck <- NULL
        selectedIPApeaks <- NULL
        plotDDAspectraCheck <- NULL
        outputDDAspectra <- NULL
        output_DDA_spectra_folder <- NULL
        massErrorPrecursor <- NULL
        DDAprocessingMode <- NULL
      }
      ##
      DDA_workflow_call <- function(inputPathPeaklist, iHRMSfilename, refMSPcreationCheck, refHRMSindexList, rawDDAspectraVar, massErrorRef, RTtoleranceRef,
                                    indexedIPApeaksCheck, selectedIPApeaks, plotDDAspectraCheck, outputDDAspectra, output_DDA_spectra_folder, input_path_hrms,
                                    massErrorPrecursor, DDAprocessingMode, NPT, refDDAtable, msLevelDDA, spectral_search_mode_option, output_DDA_MSP) {
        ##
        if (refMSPcreationCheck) {
          xRef <- refHRMSindexList[[iHRMSfilename]]
          ##
          rawDDAspectraVar[["precursorMZvec"]] <- mzRef[xRef]
          rawDDAspectraVar[["precursorRTvec"]] <- RTref[xRef]
        }
        ##
        DDA_peaklist <- DDA_rawSpectraDeconvolution(input_path_hrms, iHRMSfilename, rawDDAspectraVar, number_processing_threads = NPT)
        ##
        if (DDA_peaklist[1, 1] != 0) {
          if (refMSPcreationCheck) {
            ##
            precursorScanNumber_IDref <- do.call(rbind, lapply(xRef, function(j) {
              xprecursorScanNumber <- mzRTindexer(DDA_peaklist[, 2], DDA_peaklist[, 3], mzRef[j], RTref[j], massErrorRef, RTtoleranceRef)
              ##
              if (!is.null(xprecursorScanNumber)) {
                c(DDA_peaklist[xprecursorScanNumber, 1], j)
              }
            }))
            ##
            DDA_REF_MSP <- IDSL.CSA_referenceMSPgenerator(DDA_peaklist, refDDAtable, precursorScanNumber_IDref, msLevelDDA, spectral_search_mode = "dda", spectral_search_mode_option)
            write.table(DDA_REF_MSP, file = paste0(output_DDA_MSP, "/DDA_REF_MSP_", iHRMSfilename, ".msp"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
            ##
          } else {
            DDA_MSP <- IDSL.CSA_MSPgenerator(DDA_peaklist, msLevelDDA, spectral_search_mode = "dda", spectral_search_mode_option = "rawddaspectra", number_processing_threads = NPT)
            write.table(DDA_MSP, file = paste0(output_DDA_MSP, "/DDA_MSP_", iHRMSfilename, ".msp"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
          }
        } else {
          FSA_logRecorder(paste0("No peak was detected for `", iHRMSfilename, "`!"))
        }
        ##
        return()
      }
    } else {
      rawDDAspectraVar <- NULL
    }
    ############################################################################
    if (DDA_Analysis_Mode == "idsl.ipa") {
      spectral_search_mode_option <- NA
      ##
      inputPathPeaklist <- PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0007'), 2]
      peaklistFileNames <- dir(path = inputPathPeaklist, pattern = ".Rdata$")
      peaklistFileNames <- peaklistFileNames[grep("^peaklist_", peaklistFileNames)]
      L_PL <- length(peaklistFileNames)
      ##
      if (LHRMS > L_PL) {
        peaklistHRMSfileNames <- paste0("peaklist_", file_name_hrms, ".Rdata")
        ndPeaklists <- setdiff(peaklistHRMSfileNames, peaklistFileNames)
        ndPeaklists <- gsub("^peaklist_|.Rdata$", "", ndPeaklists)
        FSA_logRecorder("Error!!! peaklist files are not available for the following HRMS file(s):")
        for (i in ndPeaklists) {
          FSA_logRecorder(i)
        }
        stop()
      }
      ##
      indexedIPApeaksCheck <- FALSE
      if (LHRMS == 1 & !refMSPcreationCheck) {
        DDA0009 <- PARAM_DDA[which(PARAM_DDA[, 1] == "DDA0009"), 2]
        if (tolower(DDA0009) != "all") {
          indexedIPApeaksCheck <- TRUE
          selectedIPApeaks <- tryCatch(eval(parse(text = paste0("c(", DDA0009, ")"))), error = function(e){NULL})
          if (is.null(selectedIPApeaks)) {
            indexedIPApeaksCheck <- FALSE
          }
        } else {
          selectedIPApeaks <- NULL
        }
      } else {
        selectedIPApeaks <- NULL
      }
      ##
      plotDDAspectraCheck <- if (tolower(PARAM_DDA[which(PARAM_DDA[, 1] == "DDA0010"), 2]) == "yes") {TRUE} else {FALSE}
      if (plotDDAspectraCheck) {
        dev.offCheck <- TRUE
        while (dev.offCheck) {
          dev.offCheck <- tryCatch(dev.off(), error = function(e) {FALSE})
        }
        ##
        output_DDA_spectra_folder <- paste0(output_address, "/DDA_spectra")
        FSA_dir.create(output_DDA_spectra_folder, allowedUnlink = FALSE)
        opendir(output_DDA_spectra_folder)
        FSA_logRecorder("DDA spectra figures for MS2 ions are stored in the `DDA_spectra` folder!")
        ##
      } else {
        outputDDAspectra <- NULL
        output_DDA_spectra_folder <- NULL
      }
      ##
      massErrorPrecursor <- as.numeric(PARAM_DDA[which(PARAM_DDA[, 1] == "DDA0013"), 2])
      dda_processing_mode <- tolower(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0014'), 2])
      if (dda_processing_mode == "ddaspectraintegration") { # DDA spectra integration
        FSA_logRecorder("Initiated Data Dependent Acquisition (DDA) analysis on individual peaklists using the 'DDA spectra integration' method!")
        ##
        massErrorIntegration <- as.numeric(PARAM_DDA[which(PARAM_DDA[, 1] == "DDA0015"), 2])
        ##
        DDAprocessingMode <- c(dda_processing_mode, massErrorIntegration)
        ##
      } else if (dda_processing_mode == "ionfiltering") { # Ion Filtering
        FSA_logRecorder("Initiated Data Dependent Acquisition (DDA) analysis on individual peaklists using the 'ion filtering' method!")
        ##
        massErrorIonFiltering <- as.numeric(PARAM_DDA[which(PARAM_DDA[, 1] == "DDA0015"), 2])
        minPercentageDetectedScans <- as.numeric(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0016'), 2])
        rsdCutoff <- as.numeric(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0017'), 2])
        pearsonRHOthreshold <- as.numeric(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0018'), 2])
        ##
        DDAprocessingMode <- c(dda_processing_mode, massErrorIonFiltering, minPercentageDetectedScans, rsdCutoff, pearsonRHOthreshold)
        ##
      } else {# Most Intense DDA Spectra
        FSA_logRecorder("Initiated Data Dependent Acquisition (DDA) analysis on individual peaklists using the 'most intense DDA spectra' method!")
        ##
        DDAprocessingMode <- "mostintenseddaspectra"
      }
      ##
      if (refMSPcreationCheck) {
        FSA_logRecorder("Individual `.msp` files are stored in the `DDA_REF_MSP` folder!")
      } else {
        FSA_logRecorder("Individual `.msp` files are stored in the `DDA_MSP` folder!")
      }
      ##
      ##########################################################################
      ##
      DDA_workflow_call <- function(inputPathPeaklist, iHRMSfilename, refMSPcreationCheck, refHRMSindexList, rawDDAspectraVar, massErrorRef, RTtoleranceRef,
                                    indexedIPApeaksCheck, selectedIPApeaks, plotDDAspectraCheck, outputDDAspectra, output_DDA_spectra_folder, input_path_hrms,
                                    massErrorPrecursor, DDAprocessingMode, NPT, refDDAtable, msLevelDDA, spectral_search_mode_option, output_DDA_MSP) {
        ##
        peaklist <- loadRdata(paste0(inputPathPeaklist, "/peaklist_", iHRMSfilename, ".Rdata"))
        ##
        if (refMSPcreationCheck) {
          xRef <- refHRMSindexList[[iHRMSfilename]]
          selectedIPApeaks_IDref <- do.call(rbind, lapply(xRef, function(j) {
            xIPA <- mzRTindexer(peaklist[, 8], peaklist[, 3], mzRef[j], RTref[j], massErrorRef, RTtoleranceRef)
            ##
            if (!is.null(xIPA)) {
              c(xIPA, j)
            }
          }))
          ##
          if (length(selectedIPApeaks_IDref) > 0) {
            selectedIPApeaks <- unique(selectedIPApeaks_IDref[, 1])
          } else {
            selectedIPApeaks <- NULL
          }
          ##
        } else {
          if (!indexedIPApeaksCheck) {
            n_peaks <- dim(peaklist)[1]
            selectedIPApeaks <- 1:n_peaks
          }
        }
        ##
        if (plotDDAspectraCheck) {
          outputDDAspectra <- paste0(output_DDA_spectra_folder, "/DDA_spectra_", iHRMSfilename)
        }
        ##
        if (!is.null(selectedIPApeaks)) {
          DDA_peaklist <- DDA_fragmentationPeakDetection(input_path_hrms, iHRMSfilename, peaklist, selectedIPApeaks, massErrorPrecursor,
                                                         DDAprocessingMode, outputDDAspectra, number_processing_threads = NPT)
          ##
          if (DDA_peaklist[1, 1] != 0) {
            if (refMSPcreationCheck) {
              ##
              DDA_REF_MSP <- IDSL.CSA_referenceMSPgenerator(DDA_peaklist, refDDAtable, selectedIPApeaks_IDref, msLevelDDA, spectral_search_mode = "dda", spectral_search_mode_option)
              write.table(DDA_REF_MSP, file = paste0(output_DDA_MSP, "/DDA_REF_MSP_", iHRMSfilename, ".msp"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
              ##
            } else {
              DDA_MSP <- IDSL.CSA_MSPgenerator(DDA_peaklist, msLevelDDA, spectral_search_mode = "dda", spectral_search_mode_option, number_processing_threads = NPT)
              write.table(DDA_MSP, file = paste0(output_DDA_MSP, "/DDA_MSP_", iHRMSfilename, ".msp"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
            }
          } else {
            FSA_logRecorder(paste0("No peak was detected for `", iHRMSfilename, "`!"))
          }
        }
        ##
        return()
      }
    }
    ##
    ############################################################################
    ##
    if (NPT == 1 | parallelizationMode == "peakmode") {
      ##
      iCounter <- 0
      progressBARboundaries <- txtProgressBar(min = 0, max = LHRMS, initial = 0, style = 3)
      for (i in file_name_hrms) {
        ##
        null_variable <- tryCatch(DDA_workflow_call(inputPathPeaklist, iHRMSfilename = i, refMSPcreationCheck, refHRMSindexList, rawDDAspectraVar, massErrorRef, RTtoleranceRef,
                                                    indexedIPApeaksCheck, selectedIPApeaks, plotDDAspectraCheck, outputDDAspectra, output_DDA_spectra_folder, input_path_hrms,
                                                    massErrorPrecursor, DDAprocessingMode, NPT, refDDAtable, msLevelDDA, spectral_search_mode_option, output_DDA_MSP),
                                  error = function(e) {FSA_logRecorder(paste0("Problem with `", i,"`!"))})
        ##
        iCounter <- iCounter + 1
        setTxtProgressBar(progressBARboundaries, iCounter)
      }
      close(progressBARboundaries)
      ##
    } else if (parallelizationMode == "samplemode") {
      NPT0 <- NPT
      NPT <- 1
      ##
      osType <- Sys.info()[['sysname']]
      ##
      if (osType == "Linux") {
        ##
        null_variable <- mclapply(file_name_hrms, function(i) {
          ##
          tryCatch(DDA_workflow_call(inputPathPeaklist, iHRMSfilename = i, refMSPcreationCheck, refHRMSindexList, rawDDAspectraVar, massErrorRef, RTtoleranceRef,
                                     indexedIPApeaksCheck, selectedIPApeaks, plotDDAspectraCheck, outputDDAspectra, output_DDA_spectra_folder, input_path_hrms,
                                     massErrorPrecursor, DDAprocessingMode, NPT, refDDAtable, msLevelDDA, spectral_search_mode_option, output_DDA_MSP),
                   error = function(e) {FSA_logRecorder(paste0("Problem with `", i,"`!"))})
        }, mc.cores = NPT0)
        ##
        closeAllConnections()
        ##
      } else if (osType == "Windows") {
        ##
        clust <- makeCluster(NPT0)
        registerDoParallel(clust)
        ##
        null_variable <- foreach(i = file_name_hrms, .verbose = FALSE) %dopar% {
          ##
          tryCatch(DDA_workflow_call(inputPathPeaklist, iHRMSfilename = i, refMSPcreationCheck, refHRMSindexList, rawDDAspectraVar, massErrorRef, RTtoleranceRef,
                                     indexedIPApeaksCheck, selectedIPApeaks, plotDDAspectraCheck, outputDDAspectra, output_DDA_spectra_folder, input_path_hrms,
                                     massErrorPrecursor, DDAprocessingMode, NPT, refDDAtable, msLevelDDA, spectral_search_mode_option, output_DDA_MSP),
                   error = function(e) {FSA_logRecorder(paste0("Problem with `", i,"`!"))})
        }
        ##
        stopCluster(clust)
        ##
      }
      NPT <- NPT0
    }
    ##
    ############################################################################
    ##
    if (refMSPcreationCheck) {
      ref_msp_list <- dir(path = output_DDA_MSP, full.names = TRUE, pattern = ".msp$")
      if (length(ref_msp_list) > 0) {
        refMSP <- do.call(c, lapply(ref_msp_list, function(i) {
          readLines(i, warn = FALSE)
        }))
        ##
        write.table(refMSP, file = paste0(output_address, "/", mspFileName), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
        FSA_logRecorder(paste0("The reference `", mspFileName, "` file was stored in the output directory!!!"))
      }
    }
  }
  ##
  ##############################################################################
  ##### Unique tag aggregation by spectra similarity across entire samples #####
  ##############################################################################
  ##
  if (DDA0002 == "yes") {
    massError <- tryCatch(as.numeric(PARAM_DDA[which(PARAM_DDA[, 1] == "DDA0015"), 2]), warning = function(w) {as.numeric(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0021'), 2])})
    plotSpectra <- if (tolower(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0020'), 2]) == "yes") {TRUE} else {FALSE}
    allowedWeightedSpectralEntropy <- eval(parse(text = (PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0024'), 2])))
    minEntropySimilarity <- as.numeric(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0025'), 2])
    ##
    if (refMSPcreationCheck) {
      if (file.exists(paste0(output_address, "/", mspFileName))) {
        FSA_logRecorder("Initiated detecting unique DDA variants!")
        ##
        aggregateBy <- PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0019'), 2]
        xAggregateBy <- which(colnames(refDDAtable) == aggregateBy)
        if (length(xAggregateBy) == 0) {
          aggregateBy <- "Name"
        }
        FSA_logRecorder(paste0("The meta-variable for aggregation is `", aggregateBy, "`!"))
        ##
        listSimilarMSPvariants <- FSA_uniqueMSPblockTagger(path = output_address, MSPfile = mspFileName, aggregateBy, massError, RTtolerance = NA, minEntropySimilarity,
                                                           allowedNominalMass = FALSE, allowedWeightedSpectralEntropy, noiseRemovalRatio = 0, plotSpectra, number_processing_threads = NPT)
        FSA_logRecorder(paste0("Indices of similar MSP blocks for each compound are stored as `listSimilarMSPvariants.Rdata` in the `", output_address,"` folder!"))
        save(listSimilarMSPvariants, file = paste0(output_address, "/listSimilarMSPvariants.Rdata"))
        FSdb_address <- paste0(output_address, "/uniqueMSPtags_", gsub("[.]msp$|[.]Rdata$", ".Rdata", mspFileName, ignore.case = TRUE))
        FSA_logRecorder("Completed detecting unique DDA variants!")
      } else {
        FSdb_address <- ""
        FSA_logRecorder("No reference compound was detected!")
      }
      ##
    } else {
      ##
      massErrorPrecursor <- tryCatch(as.numeric(PARAM_DDA[which(PARAM_DDA[, 1] == "DDA0013"), 2]), error = function(e) {massError})
      RTtoleranceRef <- as.numeric(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0022'), 2])
      minDDAdetectionFrequency <- floor(as.numeric(PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0023'), 2])*LHRMS/100)
      MSPfile_vector <- dir(path = output_DDA_MSP, pattern = ".msp$", ignore.case = TRUE)
      ##
      FSA_logRecorder("Initiated detecting unique DDA variants!")
      FSA_uniqueMSPblockTaggerUntargeted(path = output_DDA_MSP, MSPfile_vector, minDDAdetectionFrequency, minEntropySimilarity, massError, massErrorPrecursor, RTtoleranceRef,
                                         noiseRemovalRatio = 0, allowedNominalMass = FALSE, allowedWeightedSpectralEntropy, plotSpectra, number_processing_threads = NPT)
      FSdb_address <- paste0(output_DDA_MSP, "/UNIQUETAGS/uniqueMSPtagsUntargeted.Rdata")
      FSA_logRecorder("Completed detecting unique DDA variants!")
    }
    ##
    ############################################################################
    ##
    IPA_PAxlsxCheck <- IPA_peak_alignment_folder_xlsxAnalyzer(PARAM = PARAM_DDA, PARAM_ID = 'DDA0008', checkpoint_parameter = TRUE, correctedRTcheck = FALSE, CSAcheck = TRUE, allowedVerbose = FALSE)
    PARAM_DDA <- IPA_PAxlsxCheck[[1]]
    checkpoint_parameter <- IPA_PAxlsxCheck[[2]]
    ##
    if (checkpoint_parameter & file.exists(FSdb_address)) {
      peak_alignment_folder <- PARAM_DDA[which(PARAM_DDA[, 1] == 'DDA0008'), 2]
      ##
      listXcolHeightAreaR13C <- FSdb2PeakXcolSubsetter(FSdb_address, peak_alignment_folder, metavariable = "idsl.ipa_peakid", number_processing_threads = NPT)
      if (!is.null(listXcolHeightAreaR13C)) {
        ##
        FSA_logRecorder("Initiated subsetting the peak alignment tables!")
        ##
        peak_alignment_subset <- paste0(output_address, "/peak_alignment_subset")
        FSA_dir.create(peak_alignment_subset, allowedUnlink = FALSE)
        ##
        peakXcol_subset <- listXcolHeightAreaR13C[["peakXcol"]]
        save(peakXcol_subset, file = paste0(peak_alignment_subset, "/peakXcol_subset.Rdata"))
        peakXcol_subset <- NULL
        ##
        peak_height_subset <- listXcolHeightAreaR13C[["peak_height"]]
        write.csv(peak_height_subset, file = paste0(peak_alignment_subset, "/peak_height_subset.csv"), row.names = TRUE)
        peak_height_subset <- NULL
        ##
        peak_area_subset <- listXcolHeightAreaR13C[["peak_area"]]
        write.csv(peak_area_subset, file = paste0(peak_alignment_subset, "/peak_area_subset.csv"), row.names = TRUE)
        peak_area_subset <- NULL
        ##
        peak_R13C_subset <- listXcolHeightAreaR13C[["peak_R13C"]]
        write.csv(peak_R13C_subset, file = paste0(peak_alignment_subset, "/peak_R13C_subset.csv"), row.names = TRUE)
        peak_R13C_subset <- NULL
        ##
        listXcolHeightAreaR13C <- NULL
        ##
        FSA_logRecorder("Completed subsetting the peak alignment tables!")
      }
    }
  }
  ##
  ##############################################################################
  ##
  completion_time <- Sys.time()
  FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
  required_time <- completion_time - initiation_time
  FSA_logRecorder(paste0("The required processing time was `", required_time, " ", attributes(required_time)$units, "`"))
  FSA_logRecorder(paste0(as.character(completion_time), " ", timeZone), allowedPrinting = FALSE)
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder("Completed the DDA analysis successfully!")
  FSA_logRecorder(paste0(rep("", 100), collapse = "="), allowedPrinting = FALSE)
  ##
  ##############################################################################
  ##
  return()
}