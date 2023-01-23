aligned_fragmentation_spectra_annotator <- function(PARAM_AT, output_path) {
  ##
  ##############################################################################
  ## To create log record for IDSL.CSA
  initiation_time <- Sys.time()
  timeZone <- tryCatch(Sys.timezone(), warning = function(w) {"UTC"}, error = function(e) {"UTC"})
  .logFSA <- NULL
  .logFSA <<- paste0(output_path, "/logCSA_AlignedTable.txt")
  FSA_logRecorder(paste0(rep("", 100), collapse = "="))
  FSA_logRecorder(paste0("OUTPUT:  ", output_path))
  FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
  FSA_logRecorder("Initiated generating the aligned spectra annotated table!")
  FSA_logRecorder(paste0(as.character(initiation_time), " ", timeZone))
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder(paste0(PARAM_AT[, 1], "\t", PARAM_AT[, 2]),  allowedPrinting = FALSE)
  FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  ##############################################################################
  ##
  peak_alignment_folder <- PARAM_AT[which(PARAM_AT[, 1] == 'AT0001'), 2]
  number_processing_threads <- as.numeric(PARAM_AT[which(PARAM_AT[, 1] == "AT0002"), 2])
  ##
  output_path_annotated_spectra_tables <- paste0(output_path, "/annotated_spectra_tables")
  spectra_table_list <- dir(path = output_path_annotated_spectra_tables, pattern = ".Rdata")
  ##
  peakXcol <- IDSL.IPA::loadRdata(paste0(peak_alignment_folder, "/peakXcol.Rdata"))
  ColPL <- colnames(peakXcol)
  L_peaks <- dim(peakXcol)[1]
  Lsamples3 <- dim(peakXcol)[2]
  Lsamples <- Lsamples3 - 3
  ColPL <- ColPL[4:Lsamples3]
  ##
  ##############################################################################
  ##
  seqXcolSample <- do.call(rbind, lapply(1:length(ColPL), function(i) {
    patternSampleName <- paste0("_MSP_", ColPL[i], ".msp.Rdata")
    ##
    xPatternCheck <- grep(patternSampleName, spectra_table_list)
    if (length(xPatternCheck) > 0) {
      c(i, xPatternCheck)
    } else {
      c(i, 0)
    }
  }))
  ##
  MissedPL <- which(seqXcolSample[, 2] == 0)
  if (length(MissedPL) > 0) {
    FSA_logRecorder("WARNING!!! SpectraAnnotationTables are not avialable for the following HRMS files:")
    for (i in MissedPL) {
      FSA_logRecorder(ColPL[i])
    }
    ##
    seqXcolSample <- matrix(seqXcolSample[-MissedPL, ], ncol = 2)
  }
  ##
  seqXcolSample <- matrix(seqXcolSample[order(seqXcolSample[, 2], decreasing = FALSE), ], ncol = 2)
  orderSeqSample <- seqXcolSample[, 1]
  seqSample <- seqXcolSample[, 2]
  peakXcol <- peakXcol[, c(seq(1, 3, 1), (orderSeqSample + 3))]
  ##
  ##############################################################################
  ##
  maxRankSample <- as.numeric(PARAM_AT[which(PARAM_AT[, 1] == 'AT0003'), 2])
  Ncandidate <- as.numeric(PARAM_AT[which(PARAM_AT[, 1] == 'AT0004'), 2])
  adjustFreqRankCheck <- eval(parse(text = (PARAM_AT[which(PARAM_AT[, 1] == 'AT0005'), 2])))
  ##
  ##############################################################################
  ##
  RTtolerance <- as.numeric(PARAM_AT[which(PARAM_AT[, 1] == 'AT0006'), 2])
  ##
  call_Rank_Zcol <- "call_Rank_Zcol <- function(i) {
    ##
    peak_table_id <- peakXcol[, (i + 3)]
    SpectraAnnotationTable <- IDSL.IPA::loadRdata(paste0(output_path_annotated_spectra_tables, '/', spectra_table_list[i]))
    matched_peak_ids <- as.numeric(SpectraAnnotationTable$analyte_idsl.ipa_peakid)
    x_peak_ids <- which(peak_table_id %in% unique(matched_peak_ids))
    ##
    if (length(x_peak_ids) > 0) {
      ##
      matchedMetaVariable1 <- SpectraAnnotationTable$MetaVariable1str
      if (is.null(matchedMetaVariable1)) {
        matchedMetaVariable1 <- rep('', dim(SpectraAnnotationTable)[1])
      }
      ##
      matchedMetaVariable2 <- SpectraAnnotationTable$MetaVariable2str
      if (is.null(matchedMetaVariable2)) {
        matchedMetaVariable2 <- rep('', dim(SpectraAnnotationTable)[1])
      }
      ##
      matchedMetaVariable3 <- SpectraAnnotationTable$MetaVariable3str
      if (is.null(matchedMetaVariable3)) {
        matchedMetaVariable3 <- rep('', dim(SpectraAnnotationTable)[1])
      }
      ##
      do.call(rbind, lapply(x_peak_ids, function(j) {
        ##
        x_j <- which(matched_peak_ids == peak_table_id[j])
        max_j <- min(c(maxRankSample, length(x_j)))
        ##
        jRow <- cbind(matchedMetaVariable1[x_j[1:max_j]], matchedMetaVariable2[x_j[1:max_j]], matchedMetaVariable3[x_j[1:max_j]])
        jRow <- data.frame(jRow)
        colnames(jRow) <- c('MetaVariable1', 'MetaVariable2', 'MetaVariable3')
        ##
        uMetaVariable1 <- unique(jRow$MetaVariable1)
        ##
        do.call(rbind, lapply(1:length(uMetaVariable1), function(k) {
          kRow <- subset(jRow, MetaVariable1 == uMetaVariable1[k])
          ##
          cbind(j, k, kRow[1, ])
        }))
      }))
    }
  }"
  ##
  iRandomNumber <- ceiling(runif(1, min(seqSample), max(seqSample)))
  FSA_logRecorder(paste0("The `", spectra_table_list[iRandomNumber], "` file was randomly selected for the sanity check to ensure AT0007-AT0009 meta-variables are available in the annotated tables!"))
  SAnT <- IDSL.IPA::loadRdata(paste0(output_path_annotated_spectra_tables, '/', spectra_table_list[iRandomNumber]))
  colSAnT <- colnames(SAnT)
  ##
  Meta1Variable <- PARAM_AT[which(PARAM_AT[, 1] == 'AT0007'), 2]
  Meta1Variable <- gsub("^FSDB_", "", Meta1Variable, ignore.case = TRUE)
  if (tolower(Meta1Variable) == "id") {
    MetaVariable1str <- "IDSL.FSA_FSDBreferenceID"
    Meta1Variable <- "IDSL.FSA_FSDBreferenceID"
  } else {
    MetaVariable1str <- paste0("FSDB_", tolower(Meta1Variable))
    if (length(which(colSAnT == MetaVariable1str)) == 0) {
      FSA_logRecorder("WARNING!!! The AT0007 metadata variable is not available in in the selected annotated tables, and the default metadata variable ('IDSL.FSA_FSDBreferenceID') will be used as the 1st chemical metadata variable!")
      MetaVariable1str <- "IDSL.FSA_FSDBreferenceID"
      Meta1Variable <- "IDSL.FSA_FSDBreferenceID"
    }
  }
  call_Rank_Zcol <- gsub("MetaVariable1str", MetaVariable1str, call_Rank_Zcol)
  ##
  Meta2Variable <- PARAM_AT[which(PARAM_AT[, 1] == 'AT0008'), 2]
  Meta2Variable <- gsub("^FSDB_", "", Meta2Variable, ignore.case = TRUE)
  MetaVariable2str <- paste0("FSDB_", tolower(Meta2Variable))
  if (length(which(colSAnT == MetaVariable2str)) == 0) {
    FSA_logRecorder("WARNING!!! The AT0008 metadata variable is not available in in the selected annotated tables, and the default metadata variable ('name') will be used as the 2nd chemical metadata variable!")
    MetaVariable1str <- "FSDB_name"
    Meta2Variable <- "name"
  }
  call_Rank_Zcol <- gsub("MetaVariable2str", MetaVariable2str, call_Rank_Zcol)
  ##
  Meta3Variable <- PARAM_AT[which(PARAM_AT[, 1] == 'AT0009'), 2]
  Meta3Variable <- gsub("^FSDB_", "", Meta3Variable, ignore.case = TRUE)
  MetaVariable3str <- paste0("FSDB_", tolower(Meta3Variable))
  if (length(which(colSAnT == MetaVariable3str)) == 0) {
    FSA_logRecorder("WARNING!!! The AT0009 metadata variable is not available in in the selected annotated tables, and the default metadata variable ('MSPfilename') will be used as the 3rd chemical metadata variable!")
    MetaVariable1str <- "FSDB_MSPfilename"
    Meta3Variable <- "MSPfilename"
  }
  ##
  call_Rank_Zcol <- gsub("MetaVariable3str", MetaVariable3str, call_Rank_Zcol)
  eval(parse(text = call_Rank_Zcol))
  ##
  ##############################################################################
  ##
  rep3NA00Ncandidate <- rep(c(NA, NA, NA, 0, 0), Ncandidate)
  ##
  call_calculating_median_ranks <- function(i) {
    metaVariableFreqRank <- rep3NA00Ncandidate
    ##
    if (xZcol[i, 1] != 0) {
      rankMetaVariable123 <- Rank_Zcol[xZcol[i, 1]:xZcol[i, 2], 2:5]
      rankMetaVariable123 <- rankMetaVariable123[order(rankMetaVariable123$Rank, decreasing = FALSE), ]
      ##
      uMetaVariable1 <- unique(rankMetaVariable123$MetaVariable1)
      ##
      metaVariableFreqRank5 <- do.call(rbind, lapply(uMetaVariable1, function(j) {
        jRow <- subset(rankMetaVariable123, MetaVariable1 == j)
        jFreq <- dim(jRow)[1]
        jMed <- median(jRow$Rank)
        ##
        j3T <- table(jRow$MetaVariable2)
        if (length(j3T) == 1) {
          j4T <- table(jRow$MetaVariable3)
          xj4T <- which.max(j4T)
          jX <- which(names(j4T[xj4T[1]]) == jRow$MetaVariable3)
          jX <- jX[1]
        } else {
          xj3T <- which.max(j3T)
          jX <- which(names(j3T[xj3T[1]]) == jRow$MetaVariable2)
          jX <- jX[1]
        }
        ##
        cbind(jRow[jX, 2:4], jFreq, jMed)
      }))
      ##
      if (adjustFreqRankCheck) {    # To adjust ranking and frequencies
        oderAdjustFreqRank <- order(sqrt(as.numeric(metaVariableFreqRank5[, 4]))/as.numeric(metaVariableFreqRank5[, 5]), decreasing = TRUE)
        metaVariableFreqRank5 <- metaVariableFreqRank5[oderAdjustFreqRank, ]
      } else {
        metaVariableFreqRank5 <- metaVariableFreqRank5[order(as.numeric(metaVariableFreqRank5[, 5]), decreasing = FALSE), ]
        metaVariableFreqRank5 <- metaVariableFreqRank5[order(as.numeric(metaVariableFreqRank5[, 4]), decreasing = TRUE), ]
      }
      ##
      minNcandidate <- min(Ncandidate, dim(metaVariableFreqRank5)[1])
      for (k in 1:minNcandidate) {
        metaVariableFreqRank[5*k - 4] <- metaVariableFreqRank5[k, 1]
        metaVariableFreqRank[5*k - 3] <- metaVariableFreqRank5[k, 2]
        metaVariableFreqRank[5*k - 2] <- metaVariableFreqRank5[k, 3]
        metaVariableFreqRank[5*k - 1] <- metaVariableFreqRank5[k, 4]
        metaVariableFreqRank[5*k] <- metaVariableFreqRank5[k, 5]
      }
    }
    return(metaVariableFreqRank)
  }
  ##
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    ##
    FSA_logRecorder("Initiated matching peak IDs!")
    progressBARboundaries <- txtProgressBar(min = 0, max = Lsamples, initial = 0, style = 3)
    #
    Rank_Zcol <- do.call(rbind, lapply(seqSample, function(i) {
      setTxtProgressBar(progressBARboundaries, i)
      #
      call_Rank_Zcol(i)
    }))
    close(progressBARboundaries)
    #
    Rank_Zcol <- data.frame(Rank_Zcol)
    colnames(Rank_Zcol) <- c("XID", "Rank", "MetaVariable1", "MetaVariable2", "MetaVariable3")
    Rank_Zcol$XID <- as.numeric(Rank_Zcol$XID)
    Rank_Zcol$Rank <- as.numeric(Rank_Zcol$Rank)
    Rank_Zcol <- Rank_Zcol[order(Rank_Zcol[, 1], decreasing = FALSE), ]
    rownames(Rank_Zcol) <- NULL
    xDiff <- which(diff(Rank_Zcol[, 1]) > 0)
    #
    xZcol <- matrix(rep(0, 2*L_peaks), ncol = 2)
    #
    u_peakid <- unique(Rank_Zcol[, 1])
    xZcol[u_peakid, 1] <- c(1, (xDiff + 1))
    xZcol[u_peakid, 2] <- c(xDiff, dim(Rank_Zcol)[1])
    #
    FSA_logRecorder("Completed matching peak IDs!")
    ##
    FSA_logRecorder("Initiated calculating median ranks!")
    progressBARboundaries <- txtProgressBar(min = 0, max = L_peaks, initial = 0, style = 3)
    #
    aligned_spectra <- do.call(rbind, lapply(1:L_peaks, function(i) {
      setTxtProgressBar(progressBARboundaries, i)
      #
      call_calculating_median_ranks(i)
    }))
    close(progressBARboundaries)
    Rank_Zcol <- NULL
    FSA_logRecorder("Completed calculating median ranks!")
    ##
    title_mat <- do.call(c, lapply(1:Ncandidate, function(i) {
      c(paste0(Meta1Variable, "_", i), paste0(Meta2Variable, "_", i), paste0(Meta3Variable, "_", i), paste0("Frequency_", i), paste0("MedianRank_", i))
    }))
    ##
  } else {
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      FSA_logRecorder("Initiated matching peak IDs!")
      Rank_Zcol <- do.call(rbind, mclapply(seqSample, function(i) {
        call_Rank_Zcol(i)
      }, mc.cores = number_processing_threads))
      #
      Rank_Zcol <- data.frame(Rank_Zcol)
      colnames(Rank_Zcol) <- c("XID", "Rank", "MetaVariable1", "MetaVariable2", "MetaVariable3")
      Rank_Zcol$XID <- as.numeric(Rank_Zcol$XID)
      Rank_Zcol$Rank <- as.numeric(Rank_Zcol$Rank)
      Rank_Zcol <- Rank_Zcol[order(Rank_Zcol[, 1], decreasing = FALSE), ]
      rownames(Rank_Zcol) <- NULL
      xDiff <- which(diff(Rank_Zcol[, 1]) > 0)
      #
      xZcol <- matrix(rep(0, 2*L_peaks), ncol = 2)
      #
      u_peakid <- unique(Rank_Zcol[, 1])
      xZcol[u_peakid, 1] <- c(1, (xDiff + 1))
      xZcol[u_peakid, 2] <- c(xDiff, dim(Rank_Zcol)[1])
      #
      FSA_logRecorder("Completed matching peak IDs!")
      ##
      FSA_logRecorder("Initiated calculating median ranks!")
      aligned_spectra <- do.call(rbind, mclapply(1:L_peaks, function(i) {
        call_calculating_median_ranks(i)
      }, mc.cores = number_processing_threads))
      Rank_Zcol <- NULL
      FSA_logRecorder("Completed calculating median ranks!")
      ##
      title_mat <- do.call(c, mclapply(1:Ncandidate, function(i) {
        c(paste0(Meta1Variable, "_", i), paste0(Meta2Variable, "_", i), paste0(Meta3Variable, "_", i), paste0("Frequency_", i), paste0("MedianRank_", i))
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      ##
      FSA_logRecorder("Initiated matching peak IDs!")
      Rank_Zcol <- foreach(i = seqSample, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_Rank_Zcol(i)
      }
      #
      Rank_Zcol <- data.frame(Rank_Zcol)
      colnames(Rank_Zcol) <- c("XID", "Rank", "MetaVariable1", "MetaVariable2", "MetaVariable3")
      Rank_Zcol$XID <- as.numeric(Rank_Zcol$XID)
      Rank_Zcol$Rank <- as.numeric(Rank_Zcol$Rank)
      Rank_Zcol <- Rank_Zcol[order(Rank_Zcol[, 1], decreasing = FALSE), ]
      rownames(Rank_Zcol) <- NULL
      xDiff <- which(diff(Rank_Zcol[, 1]) > 0)
      #
      xZcol <- matrix(rep(0, 2*L_peaks), ncol = 2)
      #
      u_peakid <- unique(Rank_Zcol[, 1])
      xZcol[u_peakid, 1] <- c(1, (xDiff + 1))
      xZcol[u_peakid, 2] <- c(xDiff, dim(Rank_Zcol)[1])
      #
      FSA_logRecorder("Completed matching peak IDs!")
      ##
      FSA_logRecorder("Initiated calculating median ranks!")
      aligned_spectra <- foreach(i = 1:L_peaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_calculating_median_ranks(i)
      }
      Rank_Zcol <- NULL
      FSA_logRecorder("Completed calculating median ranks!")
      ##
      FSA_logRecorder("Initiated processing the peak property table!")
      ##
      title_mat <- foreach(i = 1:Ncandidate, .combine = 'c', .verbose = FALSE) %dopar% {
        c(paste0(Meta1Variable, "_", i), paste0(Meta2Variable, "_", i), paste0(Meta3Variable, "_", i), paste0("Frequency_", i), paste0("MedianRank_", i))
      }
      ##
      stopCluster(clust)
    }
  }
  ##############################################################################
  FSA_logRecorder(paste0("Initiated grouping compounds with similar `", Meta1Variable, "` on the first hit of the aligned MS/MS table annotation!"))
  ##
  medianPeakHeight <- IDSL.IPA::loadRdata(paste0(peak_alignment_folder, "/peak_height.Rdata"))[, 4:5]
  ##
  MetaVariable1RT <- data.frame(cbind(peakXcol[, 2], aligned_spectra[, 1], medianPeakHeight[, 1]))
  MetaVariable1RT[, 1] <- as.numeric(MetaVariable1RT[, 1])
  MetaVariable1RT[, 3] <- as.numeric(MetaVariable1RT[, 3])
  OrderMetaVariable1RT <- order(MetaVariable1RT[, 3], decreasing = TRUE)
  #
  groupMetaVariable1RT <- cbind(rep(NA, L_peaks), OrderMetaVariable1RT)
  #
  MetaVariable1RT <- MetaVariable1RT[OrderMetaVariable1RT, ]
  #
  uMetaVariable1RT <- unique(MetaVariable1RT[, 2])
  uMetaVariable1RT <- setdiff(uMetaVariable1RT, NA)
  #
  LuMetaVariable1RT <- length(uMetaVariable1RT)
  if (LuMetaVariable1RT > 0) {
    progressBARboundaries <- txtProgressBar(min = 0, max = LuMetaVariable1RT, initial = 0, style = 3)
    #
    counter <- 0
    for (i in 1:LuMetaVariable1RT) {
      setTxtProgressBar(progressBARboundaries, i)
      #
      x_inchikey <- which(MetaVariable1RT[, 2] == uMetaVariable1RT[i])
      #
      for (j in x_inchikey) {
        if (!is.na(MetaVariable1RT[j, 1])) {
          x_rt <- which(abs(MetaVariable1RT[x_inchikey, 1] - MetaVariable1RT[j, 1]) <= RTtolerance)
          #
          counter <- counter + 1
          groupMetaVariable1RT[x_inchikey[x_rt], 1] <- counter
          MetaVariable1RT[x_inchikey[x_rt], 1] <- NA
        }
      }
    }
    groupMetaVariable1RT <- matrix(groupMetaVariable1RT[order(groupMetaVariable1RT[, 2]), ], ncol = 2)
    close(progressBARboundaries)
  }
  #
  groupMetaVariable1RT <- matrix(groupMetaVariable1RT[, 1], ncol = 1)
  namesGroupMetaVariable1RT <- paste0("Group_", Meta1Variable, "_number")
  FSA_logRecorder(paste0("Completed grouping compounds with similar `", Meta1Variable, "` on the first hit of the aligned MS/MS table annotation!"))
  ##
  ##############################################################################
  ##
  medianPeakArea <- IDSL.IPA::loadRdata(paste0(peak_alignment_folder, "/peak_area.Rdata"))[, 4]
  medianR13C <- IDSL.IPA::loadRdata(paste0(peak_alignment_folder, "/peak_R13C.Rdata"))[, 4]
  ##
  aligned_spectra <- cbind(peakXcol[, 1:3], medianPeakHeight[, 2], medianPeakHeight[, 1], medianPeakArea, medianR13C, groupMetaVariable1RT, aligned_spectra)
  rownames(aligned_spectra) <- NULL
  colnames(aligned_spectra) <- c("mz", "RT", "frequencyPeakXcol", "Flag", "medianPeakHeight", "medianPeakArea", "medianR13C", namesGroupMetaVariable1RT, title_mat)
  ##
  output_path_aligned_table <- paste0(output_path, "/aligned_spectra_table")
  FSA_dir.create(output_path_aligned_table, allowedUnlink = FALSE)
  ##
  FSA_logRecorder("Initiated saving the aligned spectra table!")
  save(aligned_spectra, file = paste0(output_path_aligned_table, "/aligned_spectra.Rdata"))
  write.csv(aligned_spectra, file = paste0(output_path_aligned_table, "/aligned_spectra.csv"), row.names = TRUE)
  FSA_logRecorder("Stored annotated aligned table as `aligned_spectra_table` in the `.Rdata` and `.csv` formats in the `annotated_spectra_tables` folder!")
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
  FSA_logRecorder("Completed generating the aligned spectra annotated table!")
  FSA_logRecorder(paste0(rep("", 100), collapse = "="), allowedPrinting = FALSE)
  ##
  ##############################################################################
  ##
  return()
}