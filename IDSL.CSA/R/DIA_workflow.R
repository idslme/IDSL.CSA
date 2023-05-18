DIA_workflow <- function(PARAM_DIA) {
  ##
  if (length(PARAM_DIA) == 1) {
    if (typeof(PARAM_DIA) == "character") {
      stop("Please use `IDSL.CSA_workflow('spreadsheet')` to use the IDSL.CSA package!")
    }
  }
  ##
  DIA0001 <- tolower(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0001'), 2])
  DIA0002 <- tolower(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0002'), 2])
  NPT <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == "DIA0003"), 2])
  input_path_hrms <- PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0004'), 2]
  ##
  output_address <- PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0011'), 2]
  FSA_dir.create(output_address, allowedUnlink = FALSE)
  ##
  ##############################################################################
  ## To create log record for IDSL.CSA
  initiation_time <- Sys.time()
  timeZone <- tryCatch(Sys.timezone(), warning = function(w) {"UTC"}, error = function(e) {"UTC"})
  .logFSA <- NULL
  .logFSA <<- paste0(output_address, "/logDIA_performance.txt")
  FSA_logRecorder(paste0(rep("", 100), collapse = "="))
  FSA_logRecorder("Type <<< citation('IDSL.CSA') >>> for citing this R package in publications.")
  FSA_logRecorder(paste0("mzML/mzXML/netCDF:  ", input_path_hrms))
  FSA_logRecorder(paste0("OUTPUT:  ", output_address))
  FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
  FSA_logRecorder("Initiated DIA workflow!")
  FSA_logRecorder(paste0(as.character(initiation_time), " ", timeZone))
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder(paste0(PARAM_DIA[, 1], "\t", PARAM_DIA[, 2]),  allowedPrinting = FALSE)
  FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  ##############################################################################
  ##
  ref_xlsx_file <- as.character(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0005'), 2])
  refMSPcreationCheck <- file.exists(ref_xlsx_file)
  ##
  ##############################################################################
  ##
  if (refMSPcreationCheck) {
    ##
    refDIAtable <- CSA_reference_xlsxAnalyzer(ref_xlsx_file, input_path_hrms)[[1]]
    ##
    refHRMSindexList <- FSA_R.aggregate(refDIAtable$`Filename`)
    file_name_hrms <- names(refHRMSindexList)
    mzRef <- as.numeric(refDIAtable$`PrecursorMZ`)
    RTref <- as.numeric(refDIAtable$`Precursor_RT`)
    refDIAtable$`PrecursorMZ` <- NULL
    refDIAtable$`Precursor_RT` <- NULL
    ##
    massErrorRef <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0016'), 2])
    RTtoleranceRef <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0017'), 2])
    ##
    output_DIA_MSP <- paste0(output_address, "/DIA_REF_MSP")
    tryCatch(FSA_dir.create(output_DIA_MSP, allowedUnlink = TRUE), error = function(e) {FSA_message("ERROR!!! The .msp files inside the `DIA_REF_MSP` folder may influence the final output!")})
    ##
    str_xlsx_file <- strsplit(ref_xlsx_file, "/")[[1]]
    mspFileName <- str_xlsx_file[length(str_xlsx_file)]
    mspFileName <- gsub("[.]xlsx$|[.]csv$", ".msp", mspFileName, ignore.case = TRUE)
    mspFileName <- paste0("DIA_REF_MSP_", mspFileName)
    ##
  } else {
    refHRMSindexList <- NULL
    refDIAtable <- NULL
    massErrorRef <- 0
    RTtoleranceRef <- 0
    ##
    samples_string <- PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0006'), 2]
    if (tolower(samples_string) == "all") {
      file_name_hrms <- dir(path = input_path_hrms)
      file_name_hrms <- file_name_hrms[grep(pattern = ".mzML$|.mzXML$|.CDF$", file_name_hrms, ignore.case = TRUE)]
    } else {
      file_name_hrms <- strsplit(samples_string, ";")[[1]]
    }
    ##
    output_DIA_MSP <- paste0(output_address, "/DIA_MSP")
    FSA_dir.create(output_DIA_MSP, allowedUnlink = FALSE)
  }
  ##
  ##############################################################################
  ##
  LHRMS <- length(file_name_hrms)
  if (LHRMS == 0) {
    stop(FSA_logRecorder("EMPTY HRMS FOLDER!!!"))
  }
  ##
  ##############################################################################
  ##
  if (DIA0001 == "yes") {
    inputPathPeaklist <- PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0007'), 2]
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
      DIA0009 <- PARAM_DIA[which(PARAM_DIA[, 1] == "DIA0009"), 2]
      if (tolower(DIA0009) != "all") {
        indexedIPApeaksCheck <- TRUE
        selectedIPApeaks <- tryCatch(eval(parse(text = paste0("c(", DIA0009, ")"))), error = function(e){NULL})
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
    plotEICcheck <- if (tolower(PARAM_DIA[which(PARAM_DIA[, 1] == "DIA0010"), 2]) == "yes") {TRUE} else {FALSE}
    if (plotEICcheck) {
      dev.offCheck <- TRUE
      while (dev.offCheck) {
        dev.offCheck <- tryCatch(dev.off(), error = function(e) {FALSE})
      }
      ##
      output_DIA_EICs_folder <- paste0(output_address, "/DIA_EICs")
      FSA_dir.create(output_DIA_EICs_folder, allowedUnlink = FALSE)
      FSA_logRecorder("Aligned extracted ion chromatogram (EIC) figures for deconvoluted ions are stored in the `DIA_EICs` folder!")
      ##
    } else {
      output_DIA_EICs_folder <- NULL
      outputDIAeic <- NULL
    }
    ##
    parallelizationMode <- tolower(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0012'), 2]) 
    ##
    ############################################################################
    msLevelDIA <- 2#as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0013'), 2])
    ##
    FSA_logRecorder(paste0("Initiated Data Independent Acquisition (DIA) analysis at ms level = `", msLevelDIA,"` on individual IDSL.IPA peaklists using raw spectra!")) 
    ##
    intensityThresholdFragment <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0014'), 2])
    smoothingWindowMS1 <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0015'), 2])
    ##
    if (msLevelDIA == 2) {
      smoothingWindowMS2 <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0016'), 2])
    } else {
      smoothingWindowMS2 <- 0
    }
    ##
    massError <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0017'), 2])
    scanTolerance <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0018'), 2])
    nSpline <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0019'), 2])
    topRatioPeakHeight <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0020'), 2])/100 # ratio of top peak percentage of chromatographic peaks to measure peak similarities (%)
    pearsonRHOthreshold <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0021'), 2])
    ##
    if (refMSPcreationCheck) {
      FSA_logRecorder("Individual `.msp` files are stored in the `DIA_REF_MSP` folder!")
    } else {
      FSA_logRecorder("Individual `.msp` files are stored in the `DIA_MSP` folder!")
    }
    ##
    ############################################################################
    ##
    DIA_workflow_call <- function(iHRMSfilename) {
      ##
      peaklist <- loadRdata(paste0(inputPathPeaklist, "/peaklist_", iHRMSfilename, ".Rdata"))
      ##
      if (plotEICcheck) {
        outputDIAeic <- paste0(output_DIA_EICs_folder, "/DIA_EICs_", iHRMSfilename)
      }
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
      } else {
        if (!indexedIPApeaksCheck) {
          n_peaks <- dim(peaklist)[1]
          selectedIPApeaks <- 1:n_peaks
        }
      }
      ##
      if (!is.null(selectedIPApeaks)) {
        if (msLevelDIA == 1) {
          DIA_peaklist <- DIA_MS1_fragmentationPeakDetection(input_path_hrms, iHRMSfilename, peaklist, selectedIPApeaks, massError,
                                                             smoothingWindowMS1, scanTolerance, nSpline, topRatioPeakHeight,
                                                             intensityThresholdFragment, pearsonRHOthreshold, outputDIAeic, number_processing_threads = NPT)
          
        } else if (msLevelDIA == 2) {
          DIA_peaklist <- DIA_MS2_fragmentationPeakDetection(input_path_hrms, iHRMSfilename, peaklist, selectedIPApeaks, massError,
                                                             smoothingWindowMS1, smoothingWindowMS2, scanTolerance, nSpline, topRatioPeakHeight,
                                                             intensityThresholdFragment, pearsonRHOthreshold, outputDIAeic, number_processing_threads = NPT)
        }
      } else {
        DIA_peaklist <- matrix(rep(0, 12), nrow = 1)
      }
      ##
      if (DIA_peaklist[1, 1] != 0) {
        if (refMSPcreationCheck) {
          ##
          DIA_REF_MSP <- IDSL.CSA_referenceMSPgenerator(DIA_peaklist, refDIAtable, selectedIPApeaks_IDref, msLevelDIA, spectral_search_mode = "dia", spectral_search_mode_option = NA)
          write.table(DIA_REF_MSP, file = paste0(output_DIA_MSP, "/DIA_REF_MSP_", iHRMSfilename, ".msp"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
        } else {
          ##
          DIA_MSP <- IDSL.CSA_MSPgenerator(DIA_peaklist, msLevelDIA, spectral_search_mode = "dia", spectral_search_mode_option = NA, number_processing_threads = NPT)
          write.table(DIA_MSP, file = paste0(output_DIA_MSP, "/DIA_MSP_", iHRMSfilename, ".msp"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
        }
      } else {
        FSA_logRecorder(paste0("No peak was detected for `", iHRMSfilename, "`!"))
      }
      ##
      return()
    }
    ##
    ############################################################################
    ##
    if (NPT == 1 | parallelizationMode == "peakmode") {
      ##
      iCounter <- 0
      progressBARboundaries <- txtProgressBar(min = 0, max = LHRMS, initial = 0, style = 3)
      for (iHRMSfilename in file_name_hrms) {
        ##
        null_variable <- tryCatch(DIA_workflow_call(iHRMSfilename),
                                  error = function(e) {FSA_logRecorder(paste0("Problem with `", iHRMSfilename,"`!"))})
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
      if (osType == "Windows") {
        ##
        clust <- makeCluster(NPT0)
        clusterExport(clust, setdiff(ls(), c("clust", "file_name_hrms")), envir = environment())
        ##
        null_variable <- parLapply(clust, file_name_hrms, function(iHRMSfilename) {
          ##
          tryCatch(DIA_workflow_call(iHRMSfilename),
                   error = function(e) {FSA_logRecorder(paste0("Problem with `", iHRMSfilename,"`!"))})
        })
        ##
        stopCluster(clust)
        ##
      } else {
        ##
        null_variable <- mclapply(file_name_hrms, function(iHRMSfilename) {
          ##
          tryCatch(DIA_workflow_call(iHRMSfilename),
                   error = function(e) {FSA_logRecorder(paste0("Problem with `", iHRMSfilename,"`!"))})
        }, mc.cores = NPT0)
        ##
        closeAllConnections()
        ##
      }
      NPT <- NPT0
    }
    ##
    ############################################################################
    ##
    if (refMSPcreationCheck) {
      ref_msp_list <- dir(path = output_DIA_MSP, full.names = TRUE, pattern = ".msp$")
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
  if (DIA0002 == "yes") {
    massError <- tryCatch(as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0017'), 2]), warning = function(w) {as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0024'), 2])})
    plotSpectra <- if (tolower(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0023'), 2]) == "yes") {TRUE} else {FALSE}
    allowedWeightedSpectralEntropy <- eval(parse(text = (PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0027'), 2])))
    minEntropySimilarity <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0028'), 2])
    ##
    if (refMSPcreationCheck) {
      if (file.exists(paste0(output_address, "/", mspFileName))) {
        FSA_logRecorder("Initiated detecting unique DIA variants!")
        ##
        aggregateBy <- PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0022'), 2]
        xAggregateBy <- which(colnames(refDIAtable) == aggregateBy)
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
        FSA_logRecorder("Completed detecting unique DIA variants!")
      } else {
        FSdb_address <- ""
        FSA_logRecorder("No reference compound was detected!")
      }
      ##
    } else {
      ##
      massErrorPrecursor <- massError
      RTtoleranceRef <- as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0025'), 2])
      minDIAdetectionFrequency <- floor(as.numeric(PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0026'), 2])*LHRMS/100)
      MSPfile_vector <- dir(path = output_DIA_MSP, pattern = ".msp$", ignore.case = TRUE)
      ##
      FSA_logRecorder("Initiated detecting unique DIA variants!")
      FSA_uniqueMSPblockTaggerUntargeted(path = output_DIA_MSP, MSPfile_vector, minDIAdetectionFrequency, minEntropySimilarity, massError, massErrorPrecursor, RTtoleranceRef,
                                         noiseRemovalRatio = 0, allowedNominalMass = FALSE, allowedWeightedSpectralEntropy, plotSpectra, number_processing_threads = NPT)
      FSdb_address <- paste0(output_DIA_MSP, "/UNIQUETAGS/uniqueMSPtagsUntargeted.Rdata")
      FSA_logRecorder("Completed detecting unique DIA variants!")
    }
    ##
    ############################################################################
    ##
    IPA_PAxlsxCheck <- IPA_peak_alignment_folder_xlsxAnalyzer(PARAM = PARAM_DIA, PARAM_ID = 'DIA0008', checkpoint_parameter = TRUE, correctedRTcheck = FALSE, CSAcheck = TRUE, allowedVerbose = FALSE)
    PARAM_DIA <- IPA_PAxlsxCheck[[1]]
    checkpoint_parameter <- IPA_PAxlsxCheck[[2]]
    ##
    if (checkpoint_parameter & file.exists(FSdb_address)) {
      peak_alignment_folder <- PARAM_DIA[which(PARAM_DIA[, 1] == 'DIA0008'), 2]
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
  FSA_logRecorder("Completed the DIA analysis successfully!")
  FSA_logRecorder(paste0(rep("", 100), collapse = "="), allowedPrinting = FALSE)
  ##
  ##############################################################################
  ##
  return()
}