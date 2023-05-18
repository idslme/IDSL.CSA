CSA_alignedMetaSpectraCataloger <- function(address_input_msp, peakXcol, peak_height, CSA_aligned_property_table, groupedID, minTanimotoCoefficient = 0.5, number_processing_threads = 1) {
  ##
  FSA_logRecorder("Initiated generating integrated and the most abundant aligned MSP units!")
  ##
  file_name_sample_msp <- dir(path = address_input_msp)
  file_name_sample_msp <- file_name_sample_msp[grepl(".msp$", file_name_sample_msp, ignore.case = TRUE)]
  ##
  Lsamples3 <- dim(peakXcol)[2]
  Lsamples5 <- dim(peak_height)[2]
  ##############################################################################
  ##############################################################################
  ColPL <- colnames(peakXcol)[4:Lsamples3]
  ##
  seqXcolMSP <- do.call(rbind, lapply(1:length(ColPL), function(i) {
    patternSampleName <- paste0("_MSP_", ColPL[i], ".msp")
    ##
    xPatternCheck <- grep(patternSampleName, file_name_sample_msp)
    if (length(xPatternCheck) > 0) {
      c(i, xPatternCheck)
    } else {
      c(i, 0)
    }
  }))
  ##
  MissedPL <- which(seqXcolMSP[, 2] == 0)
  nMissedPL <- length(MissedPL)
  if (nMissedPL > 0) {
    FSA_logRecorder("WARNING!!! The following MSP files are not avialable:")
    for (i in MissedPL) {
      FSA_logRecorder(ColPL[i])
    }
    seqXcolMSP <- matrix(seqXcolMSP[-MissedPL, ], ncol = 2)
    Lsamples5 <- Lsamples5 - nMissedPL
  }
  ##
  seqXcolMSP <- matrix(seqXcolMSP[order(seqXcolMSP[, 2], decreasing = FALSE), ], ncol = 2)
  orderSeqMSP <- seqXcolMSP[, 1]
  peakXcol <- peakXcol[, c(seq(1, 3, 1), (orderSeqMSP + 3))]
  peak_height <- peak_height[, c(seq(1, 5, 1), (orderSeqMSP + 5))]
  ##############################################################################
  ##############################################################################
  groupedID <- cbind(seq(1, nrow(groupedID), 1), groupedID)
  groupedID <- groupedID[(groupedID[, 3] > 0), ]
  groupedID <- groupedID[order(groupedID[, 3], decreasing = TRUE), ]
  groupedID <- groupedID[order(groupedID[, 2], decreasing = FALSE), ]
  ##
  xDiff <- which(diff(groupedID[, 2]) > 0)
  xsDiff1 <- c(1, (xDiff + 1))
  xsDiff2 <- c(xDiff, nrow(groupedID))
  ##
  mzRTalignedTable <- peakXcol[, 1:2]
  mzRTalignedTable[, 2] <- round(mzRTalignedTable[, 2], digits = 3)
  ##
  ##############################################################################
  ##
  call_MSPvectorAverage <- function(i) {
    ##
    if (!is.na(groupedID[xsDiff1[i], 3])) {
      if (groupedID[xsDiff1[i], 3] >= minTanimotoCoefficient) {
        peakIDs <- groupedID[(xsDiff1[i]:xsDiff2[i]), 1]
        ##
        median_peakIDs_height <- CSA_aligned_property_table[peakIDs, 4]
        ##
        median_peakIDs_R13C <- CSA_aligned_property_table[peakIDs, 5]
        ##
        rep0000 <- rep(0, length(peakIDs))
        ##
        peaks <- rbind(cbind(peakIDs, mzRTalignedTable[peakIDs, 1], median_peakIDs_height, groupedID[(xsDiff1[i]:xsDiff2[i]), 3]), cbind(rep0000, (mzRTalignedTable[peakIDs, 1] + 1.003354835336), median_peakIDs_R13C*median_peakIDs_height/100, rep0000))
        peaks[, 2] <- round(peaks[, 2], digits = 5)
        peaks[, 3] <- round(peaks[, 3], digits = 0)
        ##
        peaks <- peaks[order(peaks[, 3], decreasing = TRUE), ]
        ##########################################################################
        MSPid <- paste0("Name: CSA_aligned_integrated_spectra_", groupedID[xsDiff1[i], 2], "_RT_", mzRTalignedTable[peakIDs[1], 2], "\n")
        MSPid <- paste0(MSPid, "Tanimoto_CSA_aligned_cluster_ID: ", groupedID[xsDiff1[i], 2], "\n")
        MSPid <- paste0(MSPid, "Retention_time: ", mzRTalignedTable[peakIDs[1], 2], "\n")
        MSPid <- paste0(MSPid, "Tanimoto_IDSL.IPA_AlignedTable_PeakIDs: ", paste0(peaks[, 1], collapse = ","), "\n")
        MSPid <- paste0(MSPid, "Tanimoto_IDSL.IPA_PeakHeight: ", max(median_peakIDs_height), "\n")
        MSPid <- paste0(MSPid, "Tanimoto_coefficients: ", paste0(peaks[, 4], collapse = ","), "\n")
        MSPid <- paste0(MSPid, "Num Peaks: ", (xsDiff2[i] - xsDiff1[i] + 1)*2, "\n")
        MSPid <- paste0(MSPid, paste0(peaks[, 2], " ", peaks[, 3], "\n", collapse = ""), "\n")
        return(MSPid)
      }
    }
  }
  ##
  ##############################################################################
  ##
  CSAmspIPApeakExtraction <- function(path, mspFileName, peakID) {
    subset.msp <- NULL
    ##
    msp <- readLines(paste0(path, "/", mspFileName), warn = FALSE)
    msp <- c("", msp, "")
    ##
    loc_collectivePeakIDs <- grep("IDSL.IPA_Collective_PeakIDs: ", msp, ignore.case = TRUE)
    ##
    strPeakID <- paste0(" ", peakID, ",")
    ##
    xMSPb <- grep(strPeakID, msp[loc_collectivePeakIDs], ignore.case = TRUE)
    ##
    if (length(xMSPb) == 0) {
      ##
      strPeakID <- paste0(",", peakID)
      ##
      xMSPbn <- grep(strPeakID, msp[loc_collectivePeakIDs], ignore.case = TRUE)
      ##
      xMSPb <- do.call(c, lapply(xMSPbn, function(j) {
        xM <- FSA_locate_regex(msp[loc_collectivePeakIDs[j]], strPeakID)[, 2]
        if (!is.na(xM[1])) {
          xW <- which(xM == nchar(msp[loc_collectivePeakIDs[j]]))
          if (length(xW) > 0) {
            j
          }
        }
      }))
    }
    if (length(xMSPb) != 0) {
      ##
      subset.msp2 <- lapply(xMSPb, function(k) {
        xMSPblock <- loc_collectivePeakIDs[k]
        ##
        xBlank <- which(msp == "")
        ##
        x1 <- which(xBlank < xMSPblock)
        x1 <- xBlank[x1[length(x1)]]
        ##
        x2 <- which(xBlank > xMSPblock)
        x2 <- xBlank[x2[1]]
        ##
        subset.msp <- msp[x1:x2]
      })
      ##
      subset.msp <- subset.msp2[[1]]
    }
    ##
    return(subset.msp)
  }
  ##
  call_MSPvectorMAX <- function(i) {
    ##
    if (groupedID[xsDiff1[i], 3] >= minTanimotoCoefficient) {
      alignedPeakID <- groupedID[xsDiff1[i], 1]
      ##
      orderH <- order(peak_height[alignedPeakID, 6:Lsamples5], decreasing = TRUE)
      lnon0 <- as.numeric(peak_height[alignedPeakID, 3])
      ##
      k <- 0
      while (k < lnon0) {
        k <- k + 1
        xHmax <- orderH[k]
        individualPeakID <- peakXcol[alignedPeakID, (xHmax + 3)]
        subset.msp <- CSAmspIPApeakExtraction(path = address_input_msp, mspFileName = file_name_sample_msp[xHmax], peakID = individualPeakID)
        if (!is.null(subset.msp)) {
          k <- lnon0
        }
      }
      ##
      if (!is.null(subset.msp)) {
        loc_Name <- grep("Name: ", subset.msp, ignore.case = TRUE)
        ##
        if (length(loc_Name) > 0) {
          loc_Name <- loc_Name[1]
          strName <- strsplit(subset.msp[loc_Name], "_")[[1]]
          subset.msp[loc_Name] <- paste0("Name: CSA_intense_spectra_", groupedID[xsDiff1[i], 2], "_", file_name_sample_msp[xHmax], "_", strName[3], "_RT_", strName[5])
          ##
          subset.msp2 <- c(subset.msp[1:2],
                           paste0("Tanimoto_CSA_aligned_cluster_ID: ", groupedID[xsDiff1[i], 2]),
                           paste0("Tanimoto_IDSL.IPA_AlignedTable_PeakIDs: ", paste0(groupedID[(xsDiff1[i]:xsDiff2[i]), 1], collapse = ",")),
                           paste0("Tanimoto_IDSL.IPA_PeakHeight: ", peak_height[alignedPeakID, (xHmax + 2)]),
                           paste0("Tanimoto_coefficients: ", paste0(groupedID[(xsDiff1[i]:xsDiff2[i]), 3], collapse = ",")),
                           subset.msp[3:length(subset.msp)])
          ##
          return(subset.msp2)
        }
      }
    }
  }
  ##############################################################################
  ##############################################################################
  ##############################################################################
  if (number_processing_threads == 1) {
    ##
    ############################################################################
    ##
    MSPvectorAverage <- do.call(c, lapply(1:length(xsDiff1), function(i) {
      call_MSPvectorAverage(i)
    }))
    CSA_aligned_property_table <- NULL
    ##
    ############################################################################
    ##
    MSPvectorMAX <- do.call(c, lapply(1:length(xsDiff1), function(i) {
      call_MSPvectorMAX(i)
    }))
    ##
    ############################################################################
    ##
  } else {
    ##
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Windows") {
      ##
      ############################################################################
      ##
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, setdiff(ls(), c("clust")), envir = environment())
      ##
      MSPvectorAverage <- do.call(c, parLapply(clust, 1:length(xsDiff1), function(i) {
        call_MSPvectorAverage(i)
      }))
      ##
      stopCluster(clust)
      ##
      CSA_aligned_property_table <- NULL
      ##
      ############################################################################
      ##
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, setdiff(ls(), c("clust")), envir = environment())
      ##
      MSPvectorMAX <- do.call(c, parLapply(clust, 1:length(xsDiff1), function(i) {
        call_MSPvectorMAX(i)
      }))
      ##
      stopCluster(clust)
      ############################################################################
      ##
    } else {
      ##
      ##########################################################################
      ##
      MSPvectorAverage <- do.call(c, mclapply(1:length(xsDiff1), function(i) {
        call_MSPvectorAverage(i)
      }, mc.cores = number_processing_threads))
      CSA_aligned_property_table <- NULL
      ##
      ##########################################################################
      ##
      MSPvectorMAX <- do.call(c, mclapply(1:length(xsDiff1), function(i) {
        call_MSPvectorMAX(i)
      }, mc.cores = number_processing_threads))
      ##
      ##########################################################################
      ##
      closeAllConnections()
      ##
    }
  }
  ##############################################################################
  ##
  MSPvectorMAX <- MSPvectorMAX[-1]
  ##
  ##############################################################################
  ##
  listCSAaverageAlignedSpectra <- list(MSPvectorAverage, MSPvectorMAX)
  ##
  FSA_logRecorder("Completed generating integrated and the most abundant aligned MSP units!")
  ##
  return(listCSAaverageAlignedSpectra)
}