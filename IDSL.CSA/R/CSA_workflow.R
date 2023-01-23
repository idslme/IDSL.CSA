CSA_workflow <- function(PARAM_CSA) {
  ##
  CSA0001 <- tolower(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0001'), 2])
  CSA0002 <- tolower(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0002'), 2])
  CSA0003 <- tolower(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0003'), 2])
  NPT <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0004'), 2])
  input_path_hrms <- PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0005'), 2]
  ##
  output_address <- PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0011'), 2]
  FSA_dir.create(output_address, allowedUnlink = FALSE)
  opendir(output_address)
  ##
  ##############################################################################
  ## To create log record for IDSL.CSA
  initiation_time <- Sys.time()
  timeZone <- tryCatch(Sys.timezone(), warning = function(w) {"UTC"}, error = function(e) {"UTC"})
  .logFSA <- NULL
  .logFSA <<- paste0(output_address, "/logCSA_performance.txt")
  FSA_logRecorder(paste0(rep("", 100), collapse = "="))
  FSA_logRecorder(paste0("mzML/mzXML/netCDF: ", input_path_hrms))
  FSA_logRecorder(paste0("OUTPUT: ", output_address))
  FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
  FSA_logRecorder("Initiated CSA workflow!")
  FSA_logRecorder(paste0(as.character(initiation_time), " ", timeZone))
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder(paste0(PARAM_CSA[, 1], "\t", PARAM_CSA[, 2]), allowedPrinting = FALSE)
  FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  ##############################################################################
  ##
  ref_xlsx_file <- as.character(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0006'), 2])
  refMSPcreationCheck <- file.exists(ref_xlsx_file)
  ##
  ##############################################################################
  ##
  if (refMSPcreationCheck) {
    ##
    refCSAtable <- CSA_reference_xlsxAnalyzer(ref_xlsx_file, input_path_hrms)[[1]]
    ##
    refHRMSindexList <- FSA_R.aggregate(refCSAtable$`Filename`)
    file_name_hrms <- names(refHRMSindexList)
    mzRef <- as.numeric(refCSAtable$`PrecursorMZ`)
    RTref <- as.numeric(refCSAtable$`Precursor_RT`)
    refCSAtable$`PrecursorMZ` <- NULL
    refCSAtable$`Precursor_RT` <- NULL
    ##
    massErrorRef <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0025'), 2])
    RTtoleranceRef <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0026'), 2])
    refCSAtable <- data.frame(refCSAtable)
    ##
    output_CSA_MSP <- paste0(output_address, "/CSA_REF_MSP")
    tryCatch(FSA_dir.create(output_CSA_MSP, allowedUnlink = TRUE), error = function(e) {FSA_message("ERROR!!! The .msp files inside the `CSA_REF_MSP` folder may influence the final output!")})
    ##
    str_xlsx_file <- strsplit(ref_xlsx_file, "/")[[1]]
    mspFileName <- str_xlsx_file[length(str_xlsx_file)]
    mspFileName <- gsub("[.]xlsx$|[.]csv$", ".msp", mspFileName, ignore.case = TRUE)
    mspFileName <- paste0("CSA_REF_MSP_", mspFileName)
    ##
  } else {
    ##
    refHRMSindexList <- NULL
    refCSAtable <- NULL
    massErrorRef <- 0
    RTtoleranceRef <- 0
    ##
    output_CSA_MSP <- paste0(output_address, "/CSA_MSP")
    FSA_dir.create(output_CSA_MSP, allowedUnlink = FALSE)
    ##
    samples_string <- PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0007'), 2]
    if (tolower(samples_string) == "all") {
      file_name_hrms <- dir(path = input_path_hrms, pattern = ".mzML$|.mzXML$|.CDF$", ignore.case = TRUE)
    } else {
      file_name_hrms <- strsplit(samples_string, ";")[[1]]
    }
    ##
    output_CSA_adductAnnotator_folder <- paste0(output_address, "/CSA_adduct_annotation")
    FSA_dir.create(output_CSA_adductAnnotator_folder, allowedUnlink = FALSE)
    ##
  }
  ##
  ##############################################################################
  ##
  IDSL.IPA::opendir(output_CSA_MSP)
  ##
  LHRMS <- length(file_name_hrms)
  if (LHRMS == 0) {
    stop(FSA_logRecorder("EMPTY HRMS FOLDER!!!"))
  }
  ##
  ##############################################################################
  ##
  if (CSA0001 == "yes") {
    ##
    inputPathPeaklist <- PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0008'), 2]
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
    plotEICcheck <- if (tolower(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0010'), 2]) == "yes") {TRUE} else {FALSE}
    if (plotEICcheck) {
      dev.offCheck <- TRUE
      while (dev.offCheck) {
        dev.offCheck <- tryCatch(dev.off(), error = function(e) {FALSE})
      }
      ##
      output_CSA_EICs_folder <- paste0(output_address, "/CSA_EICs")
      FSA_dir.create(output_CSA_EICs_folder, allowedUnlink = FALSE)
      opendir(output_CSA_EICs_folder)
      FSA_logRecorder("Aligned extracted ion chromatogram (EIC) figures for deconvoluted ions are stored in the `CSA_EICs` folder!")
      ##
    } else {
      outputCSAeic <- NULL
      output_CSA_EICs_folder <- NULL
    }
    ############################################################################
    msLevelCSA <- 1  
    ##
    CSAanalysisMethod <- gsub(" ", "", tolower(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0012'), 2]))
    if (CSAanalysisMethod == "peaklist") {
      FSA_logRecorder("Initiated Composite Spectra Analysis (CSA) by grouping IDSL.IPA peaks on individual peaklists!")
      ##
      tempAlignedTableSubsetsFolder <- NULL
      ##
    } else if (CSAanalysisMethod == "alignedtable") {
      FSA_logRecorder("Initiated Composite Spectra Analysis (CSA) by grouping IDSL.IPA peaks on individual peaklists using aligned peak height table correlations!")
      ##
      peak_alignment_folder <- PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0009'), 2]
      peakXcol_FN <- paste0(peak_alignment_folder, "/peakXcol.Rdata")
      peakXcol <- IDSL.IPA::loadRdata(peakXcol_FN)
      namePeakXcol <- colnames(peakXcol)
      ##
      fileNameHRMScheck <- file_name_hrms[!(file_name_hrms %in% namePeakXcol)]
      if (length(fileNameHRMScheck) > 0) {
        FSA_logRecorder("peak alignement information are not available for the following HRMS file(s):")
        for (f in fileNameHRMScheck) {
          FSA_logRecorder(f)
        }
        stop()
      }
      ##
      tempAlignedTableSubsetsFolder <- paste0(output_address, "/tempAlignedTableSubsetsFolder")
      FSA_logRecorder(paste0("Initiated subsetting the `alignedPeakHeightTableCorrelationList.Rdata` for each sample! Temporary subsetted data are stored in the `", output_address, "` folder!"))
      tryCatch(FSA_dir.create(tempAlignedTableSubsetsFolder, allowedUnlink = TRUE), error = function(e) {stop(paste0("Temporary subsetted folder can't be created in the `", output_address, "` folder!"))})
      ##
      correlationList <- IDSL.IPA::loadRdata(paste0(peak_alignment_folder, "/alignedPeakHeightTableCorrelationList.Rdata"))
      ##
      progressBARboundaries <- txtProgressBar(min = 0, max = LHRMS, initial = 0, style = 3)
      for (i in 1:LHRMS) {
        xXcol <- which(namePeakXcol == file_name_hrms[i])
        peakXcolSubset <- peakXcol[, xXcol]
        xNon0Xcol <- which(peakXcolSubset != 0)
        correlationListSubset <- correlationList[xNon0Xcol]
        ##
        subsetAlignedFolder <- paste0(tempAlignedTableSubsetsFolder, "/", file_name_hrms[i])
        FSA_dir.create(subsetAlignedFolder, allowedUnlink = TRUE)
        ##
        save(xNon0Xcol, file = paste0(subsetAlignedFolder, "/xNon0Xcol.Rdata"))
        xNon0Xcol <- NULL
        ##
        save(peakXcolSubset, file = paste0(subsetAlignedFolder, "/peakXcolSubset.Rdata"))
        peakXcolSubset <- NULL
        ##
        save(correlationListSubset, file = paste0(subsetAlignedFolder, "/correlationListSubset.Rdata"))
        correlationListSubset <- NULL
        ##
        setTxtProgressBar(progressBARboundaries, i)
      }
      correlationList <- NULL
      peakXcol <- NULL
      ##
      close(progressBARboundaries)
      FSA_logRecorder("Completed subsetting the `alignedPeakHeightTableCorrelationList.Rdata`!")
      ##
    } else {
      stop(FSA_logRecorder("`CSA0016` was not detected"))
    }
    ##
    RTtolerance <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0013'), 2])
    minSNRbaseline <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0014'), 2])
    smoothingWindowMS1 <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0015'), 2])
    massError <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0016'), 2])
    scanTolerance <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0017'), 2])
    nSpline <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0018'), 2])
    topRatioPeakHeight <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0019'), 2])/100 # ratio of top peak percentage of chromatographic peaks to measure peak similarities (%)
    minIonRangeDifference <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0020'), 2])
    minNumCSApeaks <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0021'), 2])
    pearsonRHOthreshold <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0022'), 2])
    ##
    if (refMSPcreationCheck) {
      FSA_logRecorder("Individual `.msp` files are stored in the `CSA_REF_MSP` folder!")
    } else {
      FSA_logRecorder("Individual `.msp` files are stored in the `CSA_MSP` folder!")
      FSA_logRecorder("Individual adduct annotated IDSL.IPA peaklists are stored in `.Rdata` and `.csv` formats in the `CSA_adduct_annotation` folder!")
    }
    ##
    ############################################################################
    ##
    call_CSA_workflow <- function(inputPathPeaklist, iHRMSfilename, refMSPcreationCheck, refCSAtable, refHRMSindexList, massErrorRef,
                                  RTtoleranceRef, plotEICcheck, outputCSAeic, output_CSA_EICs_folder, input_path_hrms, tempAlignedTableSubsetsFolder,
                                  RTtolerance, massError, minSNRbaseline, smoothingWindowMS1, scanTolerance, nSpline, topRatioPeakHeight,
                                  minIonRangeDifference, minNumCSApeaks, pearsonRHOthreshold, msLevelCSA, CSAanalysisMethod, output_CSA_MSP) {
      ##
      peaklist <- loadRdata(paste0(inputPathPeaklist, "/peaklist_", iHRMSfilename, ".Rdata"))
      ##
      if (plotEICcheck) {
        outputCSAeic <- paste0(output_CSA_EICs_folder, "/CSA_EICs_", iHRMSfilename)
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
        selectedIPApeaks <- NULL
      }
      ##
      CSA_peaklist <- CSA_fragmentationPeakDetection(input_path_hrms, iHRMSfilename, tempAlignedTableSubsetsFolder, peaklist, selectedIPApeaks,
                                                     RTtolerance, massError, minSNRbaseline, smoothingWindowMS1, scanTolerance, nSpline,
                                                     topRatioPeakHeight, minIonRangeDifference, minNumCSApeaks, pearsonRHOthreshold, outputCSAeic)
      ##
      if (CSA_peaklist[1, 1] != 0) {
        if (refMSPcreationCheck) {
          ## To update CSA peaklist of CSA_peakListMSPgeneration to generate reference .msp files
          IDSL.IPA_Collective_PeakIDs <- CSA_peaklist$IDSL.IPA_PeakID
          CSA_peaklist <- cbind(CSA_peaklist, IDSL.IPA_Collective_PeakIDs)
          xMatchedIPApeakIDs <- which(CSA_peaklist$`IDSL.IPA_PeakID` %in% selectedIPApeaks)
          ##
          LxMatchedIPApeakIDs <- length(xMatchedIPApeakIDs)
          if (LxMatchedIPApeakIDs > 0) {
            ##
            MatchedIPApeakIDs <- CSA_peaklist$`IDSL.IPA_PeakID`[xMatchedIPApeakIDs]
            xCSApeakGroupingCSA <- CSA_peaklist$`CSApeakGrouping_ID`[xMatchedIPApeakIDs]
            ##
            CSA_peaklist <- do.call(rbind, lapply(1:LxMatchedIPApeakIDs, function(j) {
              CSArefList <- CSA_peaklist[CSA_peaklist$`CSApeakGrouping_ID` %in% xCSApeakGroupingCSA[j], ]
              CSArefList$`IDSL.IPA_PeakID` <- MatchedIPApeakIDs[j]
              CSArefList$`Precursor_INT` <- CSA_peaklist$`CSA_int_fragment`[xMatchedIPApeakIDs[j]]
              CSArefList$`PrecursorMZ` <- CSA_peaklist$`CSA_mz_fragment`[xMatchedIPApeakIDs[j]]
              ##
              CSArefList
            }))
            ##
            CSA_REF_MSP <- IDSL.CSA_referenceMSPgenerator(CSA_peaklist, refCSAtable, selectedIPApeaks_IDref, msLevelCSA, spectral_search_mode = "csa", spectral_search_mode_option = CSAanalysisMethod)
            write.table(CSA_REF_MSP, file = paste0(output_CSA_MSP, "/CSA_REF_MSP_", iHRMSfilename, ".msp"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
          }
          ##
        } else {
          peaklist <- CSA_adductAnnotator(peaklist, CSA_peaklist, massError)
          ##
          save(peaklist, file = paste0(output_CSA_adductAnnotator_folder, "/peaklist_CSA_adduct_annotation_", iHRMSfilename, ".Rdata"))
          write.csv(peaklist, file = paste0(output_CSA_adductAnnotator_folder, "/peaklist_CSA_adduct_annotation_", iHRMSfilename, ".csv"), row.names = TRUE)
          ##
          CSA_MSP <- IDSL.CSA_MSPgenerator(CSA_peaklist, msLevelCSA, spectral_search_mode = "csa", spectral_search_mode_option = CSAanalysisMethod, number_processing_threads = NPT)
          write.table(CSA_MSP, file = paste0(output_CSA_MSP, "/CSA_MSP_", iHRMSfilename, ".msp"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
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
    if (NPT == 1) {
      ##
      iCounter <- 0
      progressBARboundaries <- txtProgressBar(min = 0, max = LHRMS, initial = 0, style = 3)
      for (i in file_name_hrms) {
        ##
        null_variable <- tryCatch(call_CSA_workflow(inputPathPeaklist, iHRMSfilename = i, refMSPcreationCheck, refCSAtable, refHRMSindexList, massErrorRef,
                                                    RTtoleranceRef, plotEICcheck, outputCSAeic, output_CSA_EICs_folder, input_path_hrms, tempAlignedTableSubsetsFolder,
                                                    RTtolerance, massError, minSNRbaseline, smoothingWindowMS1, scanTolerance, nSpline, topRatioPeakHeight,
                                                    minIonRangeDifference, minNumCSApeaks, pearsonRHOthreshold, msLevelCSA, CSAanalysisMethod, output_CSA_MSP),
                                  error = function(e) {FSA_logRecorder(paste0("Problem with `", i,"`!"))})
        ##
        iCounter <- iCounter + 1
        setTxtProgressBar(progressBARboundaries, iCounter)
      }
      close(progressBARboundaries)
      ##
    } else {
      osType <- Sys.info()[['sysname']]
      ##
      if (osType == "Linux") {
        ##
        null_variable <- mclapply(file_name_hrms, function(i) {
          ##
          tryCatch(call_CSA_workflow(inputPathPeaklist, iHRMSfilename = i, refMSPcreationCheck, refCSAtable, refHRMSindexList, massErrorRef,
                                     RTtoleranceRef, plotEICcheck, outputCSAeic, output_CSA_EICs_folder, input_path_hrms, tempAlignedTableSubsetsFolder,
                                     RTtolerance, massError, minSNRbaseline, smoothingWindowMS1, scanTolerance, nSpline, topRatioPeakHeight,
                                     minIonRangeDifference, minNumCSApeaks, pearsonRHOthreshold, msLevelCSA, CSAanalysisMethod, output_CSA_MSP),
                   error = function(e) {FSA_logRecorder(paste0("Problem with `", i,"`!"))})
        }, mc.cores = NPT)
        ##
        closeAllConnections()
        ##
      } else if (osType == "Windows") {
        ##
        clust <- makeCluster(NPT)
        registerDoParallel(clust)
        ##
        null_variable <- foreach(i = file_name_hrms, .verbose = FALSE) %dopar% {
          ##
          tryCatch(call_CSA_workflow(inputPathPeaklist, iHRMSfilename = i, refMSPcreationCheck, refCSAtable, refHRMSindexList, massErrorRef,
                                     RTtoleranceRef, plotEICcheck, outputCSAeic, output_CSA_EICs_folder, input_path_hrms, tempAlignedTableSubsetsFolder,
                                     RTtolerance, massError, minSNRbaseline, smoothingWindowMS1, scanTolerance, nSpline, topRatioPeakHeight,
                                     minIonRangeDifference, minNumCSApeaks, pearsonRHOthreshold, msLevelCSA, CSAanalysisMethod, output_CSA_MSP),
                   error = function(e) {FSA_logRecorder(paste0("Problem with `", i,"`!"))})
        }
        ##
        stopCluster(clust)
        ##
      }
    }
    ##
    ############################################################################
    ##
    if (!is.null(tempAlignedTableSubsetsFolder)) {
      tryCatch(unlink(tempAlignedTableSubsetsFolder, recursive = TRUE), error = function(e) {FSA_logRecorder(paste0("Can't delete `", tempAlignedTableSubsetsFolder, "`!"))})
    }
    ##
    ############################################################################
    ##
    if (refMSPcreationCheck) {
      ref_msp_list <- dir(path = output_CSA_MSP, full.names = TRUE, pattern = ".msp$")
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
  if (CSA0002 == "yes") {
    massError <- tryCatch(as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0016'), 2]), warning = function(w) {as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0025'), 2])})
    plotSpectra <- if (tolower(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0024'), 2]) == "yes") {TRUE} else {FALSE}
    allowedWeightedSpectralEntropy <- eval(parse(text = (PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0028'), 2])))
    minEntropySimilarity <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0029'), 2])
    ##
    if (refMSPcreationCheck) {
      if (file.exists(paste0(output_address, "/", mspFileName))) {
        FSA_logRecorder("Initiated detecting unique CSA variants!")
        ##
        aggregateBy <- PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0023'), 2]
        xAggregateBy <- which(colnames(refCSAtable) == aggregateBy)
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
        FSA_logRecorder("Completed detecting unique CSA variants!")
      } else {
        FSdb_address <- ""
        FSA_logRecorder("No reference compound was detected!")
      }
      ##
    } else {
      ##
      RTtoleranceRef <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0026'), 2])
      minCSAdetectionFrequency <- floor(as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0027'), 2])*LHRMS/100)
      MSPfile_vector <- dir(path = output_CSA_MSP, pattern = ".msp$", ignore.case = TRUE)
      ##
      FSA_logRecorder("Initiated detecting unique CSA variants!")
      FSA_uniqueMSPblockTaggerUntargeted(path = output_CSA_MSP, MSPfile_vector, minCSAdetectionFrequency, minEntropySimilarity, massError, massErrorPrecursor = NA, RTtoleranceRef,
                                         noiseRemovalRatio = 0, allowedNominalMass = FALSE, allowedWeightedSpectralEntropy, plotSpectra, number_processing_threads = NPT)
      FSdb_address <- paste0(output_CSA_MSP, "/UNIQUETAGS/uniqueMSPtagsUntargeted.Rdata")
      FSA_logRecorder("Completed detecting unique CSA variants!")
    }
    ##
    ############################################################################
    ##
    IPA_PAxlsxCheck <- IPA_peak_alignment_folder_xlsxAnalyzer(PARAM = PARAM_CSA, PARAM_ID = 'CSA0009', checkpoint_parameter = TRUE, correctedRTcheck = FALSE, CSAcheck = TRUE, allowedVerbose = FALSE)
    PARAM_CSA <- IPA_PAxlsxCheck[[1]]
    checkpoint_parameter <- IPA_PAxlsxCheck[[2]]
    ##
    if (checkpoint_parameter & file.exists(FSdb_address)) {
      ##
      peak_alignment_folder <- PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0009'), 2]
      ##
      listXcolHeightAreaR13C <- FSdb2PeakXcolSubsetter(FSdb_address, peak_alignment_folder, metavariable = "idsl.ipa_collective_peakids", number_processing_threads = NPT)
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
  #################### CSA aggregation on the aligned table ####################
  ##############################################################################
  ##
  if (CSA0003 == "yes") {
    if (length(dir(path = output_CSA_MSP, pattern = ".msp$", ignore.case = TRUE)) > 0) {
      ##
      peak_alignment_folder <- PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0009'), 2]
      peakXcol_FN <- paste0(peak_alignment_folder, "/peakXcol.Rdata")
      peakXcol <- IDSL.IPA::loadRdata(peakXcol_FN)
      ##
      RTtolerance_AT <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0030'), 2])
      minPercenetageDetection <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0031'), 2])
      minNumberFragments <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0032'), 2])
      minTanimotoCoefficient1 <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0033'), 2])
      ##
      codetectedIDTC <- CSA_alignedPeaksTanimotoCoefficientCalculator(output_CSA_MSP, peakXcol, minPercenetageDetection, minNumberFragments,
                                                                      minTanimotoCoefficient1, RTtolerance_AT, number_processing_threads = NPT)
      ##
      peak_height <- IDSL.IPA::loadRdata(paste0(peak_alignment_folder, "/peak_height.Rdata"))
      peak_R13C <- IDSL.IPA::loadRdata(paste0(peak_alignment_folder, "/peak_R13C.Rdata"))
      ##
      CSA_aligned_table <- cbind(peak_height[, 1:4], peak_R13C[, 4:5], codetectedIDTC)
      CSA_aligned_table <- data.frame(CSA_aligned_table)
      colnames(CSA_aligned_table) <- c("m/z", "RT", "IPAdetectionFrequency", "medianPeakHeight", "medianR13C", "Flag", "coDetectedGroupingID", "TanimotoCoefficient", "CSAdetectionFrequency")
      rownames(CSA_aligned_table) <- NULL
      ##
      output_path_aligned_table <- paste0(output_address, "/aligned_spectra_table")
      FSA_dir.create(output_path_aligned_table, allowedUnlink = FALSE)
      opendir(output_path_aligned_table)
      ##    
      save(CSA_aligned_table, file = paste0(output_path_aligned_table, "/CSA_aligned_table.Rdata"))
      write.csv(CSA_aligned_table, file = paste0(output_path_aligned_table, "/CSA_aligned_table.csv"), row.names = TRUE)
      ##
      FSA_logRecorder("Co-detected aligned peaks were stored as `CSA_aligned_table` in `.Rdata` and `.csv` formats in the `aligned_spectra_table` folder!")
      ##
      ##########################################################################
      ##
      minTanimotoCoefficient2 <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0034'), 2])
      ##
      listCSAaverageAlignedSpectra <- CSA_alignedMetaSpectraCataloger(output_CSA_MSP, peakXcol, peak_height, CSA_aligned_table, codetectedIDTC,
                                                                      minTanimotoCoefficient2, number_processing_threads = NPT)
      peakXcol <- NULL
      peak_height <- NULL
      CSA_aligned_table <- NULL
      ##
      output_path_aligned_table_integrated <- paste0(output_path_aligned_table, "/MSP_integrated_aligned_spectra")
      ##
      FSA_dir.create(output_path_aligned_table_integrated, allowedUnlink = FALSE)
      write.table(listCSAaverageAlignedSpectra[[1]], file = paste0(output_path_aligned_table_integrated, "/MSP_integrated_aligned_spectra.msp"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
      ##
      FSdb_integrated_aligned_spectra <- msp2FSdb(path = output_path_aligned_table_integrated, MSPfile_vector = "MSP_integrated_aligned_spectra.msp", massIntegrationWindow = massError,
                                                  allowedNominalMass = FALSE, allowedWeightedSpectralEntropy, noiseRemovalRatio = 0.01, number_processing_threads = NPT)
      save(FSdb_integrated_aligned_spectra, file = paste0(output_path_aligned_table_integrated, "/FSdb_integrated_aligned_spectra.Rdata"))
      
      ##
      FSA_logRecorder("Stored integrated aligned MSP units in the `/aligned_spectra_table/MSP_integrated_aligned_spectra` folder!")
      ##
      output_path_aligned_table_abundant <- paste0(output_path_aligned_table, "/MSP_most_abundant_aligned_spectra")
      ##
      FSA_dir.create(output_path_aligned_table_abundant, allowedUnlink = FALSE)
      write.table(listCSAaverageAlignedSpectra[[2]], file = paste0(output_path_aligned_table_abundant, "/MSP_most_abundant_aligned_spectra.msp"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
      FSA_logRecorder("Stored abundant aligned MSP units in the `/aligned_spectra_table/MSP_most_abundant_aligned_spectra` folder!")
      ##
      FSdb_most_abundant_aligned_spectra <- msp2FSdb(path = output_path_aligned_table_abundant, MSPfile_vector = "MSP_most_abundant_aligned_spectra.msp", massIntegrationWindow = massError,
                                                     allowedNominalMass = FALSE, allowedWeightedSpectralEntropy, noiseRemovalRatio = 0.01, number_processing_threads = NPT)
      save(FSdb_most_abundant_aligned_spectra, file = paste0(output_path_aligned_table_abundant, "/FSdb_most_abundant_aligned_spectra.Rdata"))
      ##
      ##########################################################################
      #################### Plot aligned CSA variant spectra ####################
      ##########################################################################
      ##
      CSA0035 <- tolower(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0035'), 2])
      if (CSA0035 == "yes") {
        alignedCSAvariantFolder <- paste0(output_path_aligned_table_integrated, "/spectra_folder/")
        ##
        FSA_logRecorder(paste0("Tanimato integrated aligned CSA spectra figures are stored in the `", output_path_aligned_table_integrated,"` folder!"))
        ##
        FSA_plotFSdb2Spectra(path = alignedCSAvariantFolder, allowedUnlink = TRUE, annexName = "Aligned_IAS_CSA",
                             FSdb = FSdb_integrated_aligned_spectra, selectedFSdbIDs = NULL, number_processing_threads = NPT, allowedVerbose = TRUE)
        ##
        ########################################################################
        ##
        alignedCSAvariantFolder <- paste0(output_path_aligned_table_abundant, "/spectra_folder/")
        ##
        FSA_logRecorder(paste0("Tanimato abundant aligned CSA spectra figures are stored in the `", output_path_aligned_table_integrated,"` folder!"))
        ##
        FSA_plotFSdb2Spectra(path = alignedCSAvariantFolder, allowedUnlink = TRUE, annexName = "Aligned_MAAS_CSA",
                             FSdb = FSdb_most_abundant_aligned_spectra, selectedFSdbIDs = NULL, number_processing_threads = NPT, allowedVerbose = TRUE)
      }
      ##
      ##########################################################################
      ########################### Cytoscape network ############################
      ##########################################################################
      ##
      massError <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0036'), 2])
      allowedWeightedSpectralEntropy <- eval(parse(text = PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0037'), 2]))
      minEntropySimilarity_AT <- as.numeric(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0038'), 2])
      ##
      ##########################################################################
      ##
      listCytoscape_integrated <- FSA_msp2Cytoscape(path = output_path_aligned_table_integrated, MSPfile = "FSdb_integrated_aligned_spectra.Rdata", mspVariableVector = c("Retention_time", "Tanimoto_IDSL.IPA_PeakHeight"),
                                                    mspNodeID = "Tanimoto_CSA_aligned_cluster_ID", massError, RTtolerance = NA, minEntropySimilarity_AT, noiseRemovalRatio = 0.01, allowedNominalMass = FALSE,
                                                    allowedWeightedSpectralEntropy, number_processing_threads = NPT)
      ##
      write.table(listCytoscape_integrated[["node_attributes_dataFrame"]], paste0(output_path_aligned_table_integrated, "/node_attributes_dataFrame_unAnnotated.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
      ##
      write.table(listCytoscape_integrated[["correlation_network"]], paste0(output_path_aligned_table_integrated, "/correlation_network.sif"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
      ##
      write.table(listCytoscape_integrated[["edge_dataFrame"]], paste0(output_path_aligned_table_integrated, "/edge_dataFrame.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
      ##
      ##########################################################################
      ##########################################################################
      ##
      listCytoscape_abundant <- FSA_msp2Cytoscape(path = output_path_aligned_table_abundant, MSPfile = "FSdb_most_abundant_aligned_spectra.Rdata", mspVariableVector = c("Retention_time", "Tanimoto_IDSL.IPA_PeakHeight", "Name"),
                                                  mspNodeID = "Tanimoto_CSA_aligned_cluster_ID", massError, RTtolerance = NA, minEntropySimilarity_AT, noiseRemovalRatio = 0.01, allowedNominalMass = FALSE,
                                                  allowedWeightedSpectralEntropy, number_processing_threads = NPT)
      ##
      ##########################################################################
      ##########################################################################
      ## Remove redundant MSP blocks
      node_attributes <- listCytoscape_abundant[["node_attributes_dataFrame"]]
      if (!is.null(node_attributes)) {
        IPAheight <- as.numeric(node_attributes$`Tanimoto_IDSL.IPA_PeakHeight`)
        retentionTimeNodeAttributes <- as.numeric(node_attributes$`Retention_time`)
        lNode <- length(IPAheight)
        ##
        if (lNode > 0) {
          nodeIntRT <- cbind(seq(1, lNode, 1), IPAheight, retentionTimeNodeAttributes)
          orderNodeIntRT <- order(nodeIntRT[, 2], decreasing = TRUE)
          ##
          tCounter <- 0
          matchedXtable <- rep(0, lNode)
          for (t in orderNodeIntRT) {
            if (nodeIntRT[t, 1] != 0) {
              x_t <- which(abs(nodeIntRT[t, 3] - nodeIntRT[, 3]) <= RTtolerance_AT)
              xMax <- which.max(nodeIntRT[x_t, 2])
              ##
              tCounter <- tCounter + 1
              matchedXtable[tCounter] <- nodeIntRT[x_t[xMax[1]], 1]
              ##
              nodeIntRT[x_t, ] <- 0
            }
          }
          ##
          matchedXtable <- matchedXtable[1:tCounter]
          nodeid <- node_attributes$nodeid[matchedXtable]
          ##
          correlation_network <- listCytoscape_abundant[["correlation_network"]]
          ##
          xSpecNet <- which((correlation_network[, 1] %in% nodeid) | (correlation_network[, 3] %in% nodeid))
          ##
          edge_dataFrame <- listCytoscape_abundant[["edge_dataFrame"]]
          ##
          listCytoscape_abundant[["node_attributes_dataFrame"]] <- node_attributes[matchedXtable, ]
          listCytoscape_abundant[["edge_dataFrame"]] <- edge_dataFrame[xSpecNet, ]
          listCytoscape_abundant[["correlation_network"]] <- correlation_network[xSpecNet, ]
          ##
          ######################################################################
          ######################################################################
          ##
          write.table(listCytoscape_abundant[["node_attributes_dataFrame"]], paste0(output_path_aligned_table_abundant, "/node_attributes_dataFrame_unAnnotated.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
          ##
          write.table(listCytoscape_abundant[["correlation_network"]], paste0(output_path_aligned_table_abundant, "/correlation_network.sif"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
          ##
          write.table(listCytoscape_abundant[["edge_dataFrame"]], paste0(output_path_aligned_table_abundant, "/edge_dataFrame.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
          ##
          ######################################################################
          ######################################################################
          ##
          FSA_logRecorder("Stored `node_attributes_dataFrame_unAnnotated.txt`, `correlation_network.sif`, and `edge_dataFrame.txt` files for network analysis by Cytoscape software!")
        }
      }
      ##
      FSdb_file <- as.character(PARAM_CSA[which(PARAM_CSA[, 1] == 'CSA0039'), 2])
      FSdb_file <- gsub("\\", "/", FSdb_file, fixed = TRUE)
      ##
      if (file.exists(FSdb_file)) {
        nAnnotation <- 30 # Number of annotations
        ##
        alignedTableIntegratedMSPcheck <- (length(dir(path = output_path_aligned_table_integrated, pattern = ".msp$", ignore.case = TRUE)) > 0)
        alignedTableMSPabundantCheck <- (length(dir(path = output_path_aligned_table_abundant, pattern = ".msp$", ignore.case = TRUE)) > 0)
        ##
        if (alignedTableIntegratedMSPcheck | alignedTableMSPabundantCheck) {
          FSA_logRecorder("Initiated annotating `.msp` network of spectra using FSDB in the CSA0039 using default values!")
          ######################################################################
          PARAM_SPEC <- IDSL.IPA::loadRdata(paste0(system.file("data", package = "IDSL.CSA"), "/CSA_PARAM_SPEC.rda"))
          ##
          PARAM_SPEC[which(PARAM_SPEC[, 1] == "SPEC0001"), 2] <- 1 # Number of parallel processing threads
          PARAM_SPEC[which(PARAM_SPEC[, 1] == "SPEC0012"), 2] <- massError
          PARAM_SPEC[which(PARAM_SPEC[, 1] == "SPEC0013"), 2] <- allowedWeightedSpectralEntropy
          PARAM_SPEC[which(PARAM_SPEC[, 1] == "SPEC0016"), 2] <- minEntropySimilarity_AT
          ######################################################################
          ##
          libFSdb <- IDSL.IPA::loadRdata(FSdb_file)
          ##
          ######################################################################
          ##
          get_annotated_node_attributes_dataFrame <- function(annotatationTable, node_attributes_dataFrame, nAnnotation) {
            ##
            repN <- rep("", nAnnotation)
            ##
            listIDannotatationTable <- FSA_R.aggregate(annotatationTable$analyte_tanimoto_csa_aligned_cluster_id)
            ##
            annotationNodeAttributes <- do.call(rbind, lapply(node_attributes_dataFrame[, 1], function(i) {
              ##
              rowIndex <- listIDannotatationTable[[i]]
              nRowIndex <- length(rowIndex)
              if (nRowIndex > 0) {
                ##
                minK <- min(c(nAnnotation, nRowIndex))
                maxK <- max(c(0, (nAnnotation - nRowIndex)))
                ##
                c(annotatationTable$FSDB_name[rowIndex[c(1:minK)]], rep("", maxK))
              } else {
                repN
              }
            }))
            ##
            annotationNodeAttributes <- data.frame(annotationNodeAttributes)
            colnames(annotationNodeAttributes) <- paste0("Annotation_", seq(1, nAnnotation, 1))
            node_attributes_dataFrame <- cbind(node_attributes_dataFrame, annotationNodeAttributes)
            ##
            return(node_attributes_dataFrame)
          }
          ##
          ######################################################################
          ##
          if (alignedTableIntegratedMSPcheck) {
            FSA_msp_annotator(PARAM_SPEC, libFSdb, output_path_aligned_table_integrated, output_path_aligned_table_integrated, allowedVerbose = FALSE)
            ##
            annotatationTable <- IDSL.IPA::loadRdata(paste0(output_path_aligned_table_integrated, "/annotated_spectra_tables/SpectraAnnotationTable_MSP_integrated_aligned_spectra.msp.Rdata"))
            ##
            node_attributes_dataFrame <- get_annotated_node_attributes_dataFrame(annotatationTable, listCytoscape_integrated[["node_attributes_dataFrame"]], nAnnotation)
            ##
            write.table(node_attributes_dataFrame, paste0(output_path_aligned_table_integrated, "/node_attributes_dataFrame.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
          }
          ##
          ######################################################################
          ##
          if (alignedTableMSPabundantCheck) {
            FSA_msp_annotator(PARAM_SPEC, libFSdb, output_path_aligned_table_abundant, output_path_aligned_table_abundant, allowedVerbose = FALSE)
            ##
            annotatationTable <- IDSL.IPA::loadRdata(paste0(output_path_aligned_table_abundant, "/annotated_spectra_tables/SpectraAnnotationTable_MSP_most_abundant_aligned_spectra.msp.Rdata"))
            ##
            node_attributes_dataFrame <- get_annotated_node_attributes_dataFrame(annotatationTable, listCytoscape_abundant[["node_attributes_dataFrame"]], nAnnotation)
            ##
            write.table(node_attributes_dataFrame, paste0(output_path_aligned_table_abundant, "/node_attributes_dataFrame.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
          }
          ##
          ######################################################################
          ##
          if (exists('annotatationTable')) {
            FSA_logRecorder("Completed annotating `.msp` network of spectra using FSDB in the CSA0039 using default values!")
            FSA_logRecorder("Stored annotated `node_attributes_dataFrame.txt` files for network analysis by Cytoscape software!")
          }
        }
      }
    } else {
      FSA_logRecorder("No `.msp` file was detected!")
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
  FSA_logRecorder("Completed the CSA analysis successfully!")
  FSA_logRecorder(paste0(rep("", 100), collapse = "="), allowedPrinting = FALSE)
  ##
  ##############################################################################
  ##
  return()
}