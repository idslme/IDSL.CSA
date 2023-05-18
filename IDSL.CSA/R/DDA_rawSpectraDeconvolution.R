DDA_rawSpectraDeconvolution <- function(DDA_hrms_address, DDA_hrms_file, rawDDAspectraVar = NULL, number_processing_threads = 1) {
  ##
  DDA_peaklist <- data.frame(matrix(rep(0, 10), ncol = 10))
  ##############################################################################
  if (is.null(rawDDAspectraVar)) {
    rawDDAFilteringTrue <- FALSE
  } else {
    rawDDAFilteringTrue <- TRUE
    precursorMZvec <- rawDDAspectraVar[["precursorMZvec"]]
    precursorRTvec <- rawDDAspectraVar[["precursorRTvec"]]
    massError <- rawDDAspectraVar[["massError"]]
    RTtolerance <- rawDDAspectraVar[["RTtolerance"]]
    lengthVec <- length(precursorMZvec)
  }
  ##############################################################################
  p2l <- IDSL.MXP::peak2list(DDA_hrms_address, DDA_hrms_file)
  scanTable <- p2l[["scanTable"]] # this gets table of details for each spectra
  spectraList <- p2l[["spectraList"]] # this gets the spectra values
  p2l <- NULL
  #
  x_MS <- which(scanTable$peaksCount > 0) # peaks from soft ionization channel ## some files may not have data in the column re-calibration period.
  spectraList <- spectraList[x_MS]
  scanTable <- scanTable[x_MS, ]
  #
  x_MS2 <- which(scanTable$msLevel == 2)
  ##
  if (length(x_MS2) > 0) {
    ##
    precursorCE <- as.matrix(scanTable$collisionEnergy) # Collision energy
    MS_polarity <- as.numeric(scanTable$polarity)
    precursorIntensity <- as.matrix(scanTable$precursorIntensity)
    precursorMZ <- as.matrix(scanTable$precursorMZ)
    RetentionTime <- scanTable$retentionTime
    ##
    DDA_peaklist_call <- function(i) {
      if (!is.na(precursorMZ[i])) {
        MS2Peaks <- spectraList[[x_MS2[i]]]
        MS2Peaks <- MS2Peaks[order(MS2Peaks[, 2], decreasing = TRUE), ]
        L_MS2Peaks <- dim(MS2Peaks)[1]
        ##
        spectralEntropy <- round(spectral_entropy_calculator(MS2Peaks, allowedWeightedSpectralEntropy = TRUE, noiseRemovalRatio = 1e-16)[[1]], 5)
        ##
        cbind(rep(i , L_MS2Peaks), rep(precursorMZ[i], L_MS2Peaks), rep(RetentionTime[i], L_MS2Peaks), rep(precursorIntensity[i], L_MS2Peaks),
              MS2Peaks, rep(spectralEntropy, L_MS2Peaks), rep(MS_polarity[i], L_MS2Peaks), rep(precursorCE[i], L_MS2Peaks), rep(1, L_MS2Peaks))
      }
    }
    ############################################################################
    if (number_processing_threads == 1) {
      ##
      if (rawDDAFilteringTrue) {
        xScanTable <- do.call(c, lapply(1:lengthVec, function(i) {
          mzRTindexer(precursorMZ, RetentionTime, precursorMZvec[i], precursorRTvec[i], massError, RTtolerance)
        }))
      } else {
        xScanTable <- seq(1, dim(scanTable)[1], 1)
      }
      ##
      if (!is.null(xScanTable)) {
        DDA_peaklist <- do.call(rbind, lapply(xScanTable, function(i) {
          DDA_peaklist_call(i)
        }))
      }
      ##
    } else {
      osType <- Sys.info()[['sysname']]
      if (osType == "Windows") {
        ##
        if (rawDDAFilteringTrue) {
          ##
          clust <- makeCluster(number_processing_threads)
          clusterExport(clust, setdiff(ls(), c("clust", "lengthVec")), envir = environment())
          ##
          xScanTable <- do.call(c, parLapply(clust, 1:lengthVec, function(i) {
            mzRTindexer(precursorMZ, RetentionTime, precursorMZvec[i], precursorRTvec[i], massError, RTtolerance)
          }))
          ##
          stopCluster(clust)
          ##
        } else {
          xScanTable <- seq(1, dim(scanTable)[1], 1)
        }
        ##
        if (!is.null(xScanTable)) {
          ##
          clust <- makeCluster(number_processing_threads)
          clusterExport(clust, setdiff(ls(), c("clust", "xScanTable")), envir = environment())
          ##
          DDA_peaklist <- do.call(rbind, parLapply(clust, xScanTable, function(i) {
            DDA_peaklist_call(i)
          }))
          ##
          stopCluster(clust)
          ##
        }
        ##
      } else {
        ##
        if (rawDDAFilteringTrue) {
          xScanTable <- do.call(c, mclapply(1:lengthVec, function(i) {
            mzRTindexer(precursorMZ, RetentionTime, precursorMZvec[i], precursorRTvec[i], massError, RTtolerance)
          }, mc.cores = number_processing_threads))
        } else {
          xScanTable <- seq(1, dim(scanTable)[1], 1)
        }
        ##
        if (!is.null(xScanTable)) {
          DDA_peaklist <- do.call(rbind, mclapply(xScanTable, function(i) {
            DDA_peaklist_call(i)
          }, mc.cores = number_processing_threads))
        }
        ##
        closeAllConnections()
      }
    }
    ##
    if (!is.null(DDA_peaklist)) {
      DDA_peaklist <- data.frame(matrix(DDA_peaklist, ncol = 10))
      DDA_peaklist[, 1] <- as.numeric(DDA_peaklist[, 1])
      DDA_peaklist[, 2] <- round(as.numeric(DDA_peaklist[, 2]), 5)
      DDA_peaklist[, 3] <- round(as.numeric(DDA_peaklist[, 3]), 3)
      DDA_peaklist[, 4] <- round(as.numeric(DDA_peaklist[, 4]), 0)
      DDA_peaklist[, 5] <- round(as.numeric(DDA_peaklist[, 5]), 5)
      DDA_peaklist[, 6] <- round(as.numeric(DDA_peaklist[, 6]), 0)
      DDA_peaklist[, 7] <- round(as.numeric(DDA_peaklist[, 7]), 5)
    } else {
      DDA_peaklist <- data.frame(matrix(rep(0, 10), ncol = 10))
    }
  }
  ##
  rownames(DDA_peaklist) <- NULL
  colnames(DDA_peaklist) <- c("PrecursorScanNumber", "PrecursorMZ", "Precursor_RT", "Precursor_Intensity", "CSA_mz_fragment", "CSA_int_fragment", "Weighted_spectral_entropy_0noiseRemoval", "Ion_mode", "Collision_energy", "Count_DDA_scan")
  ##
  return(DDA_peaklist)
}