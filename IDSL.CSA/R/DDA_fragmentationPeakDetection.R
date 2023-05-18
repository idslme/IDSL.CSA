DDA_fragmentationPeakDetection <- function(DDA_hrms_address, DDA_hrms_file, peaklist, selectedIPApeaks,
                                           massErrorPrecursor, DDAprocessingMode = 'MostIntenseDDAspectra',
                                           outputDDAspectra = NULL, number_processing_threads = 1) {
  ##
  DDA_peaklist <- matrix(rep(0, 10), ncol = 10)
  ##############################################################################
  DDAspectraIntegrationCheck <- FALSE
  DDAionFilteringCheck <- FALSE
  ##
  if (tolower(DDAprocessingMode[1]) == 'ddaspectraintegration') {
    DDAspectraIntegrationCheck <- TRUE
    massErrorIntegration <- as.numeric(DDAprocessingMode[2])
    ##
  } else if (tolower(DDAprocessingMode[1]) == 'ionfiltering') {
    DDAionFilteringCheck <- TRUE
    massErrorIonFiltering <- as.numeric(DDAprocessingMode[2])
    minPercentageDetectedScans <- as.numeric(DDAprocessingMode[3])
    rsdCutoff <- as.numeric(DDAprocessingMode[4])
    pearsonRHOthreshold <- as.numeric(DDAprocessingMode[5])
  }
  ##
  if (is.null(outputDDAspectra)) {
    plotDDAspectraCheck <- FALSE
  } else {
    plotDDAspectraCheck <- TRUE
    ##
    FSA_dir.create(outputDDAspectra, allowedUnlink = FALSE)
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
  x_MS1 <- which(scanTable$msLevel == 1)
  x_MS2 <- which(scanTable$msLevel == 2)
  if (length(x_MS2) > 0) {
    ##
    precursorCE <- as.matrix(scanTable$collisionEnergy) # Collision energy
    MS_polarity <- as.matrix(scanTable$polarity)
    ##
    precursorIntensity <- as.matrix(scanTable$precursorIntensity)
    ##
    precursorMZ <- as.matrix(scanTable$precursorMZ)
    ############################################################################
    scanNumberStartPL <- x_MS1[peaklist[, 1]]
    scanNumberEndPL <- x_MS1[peaklist[, 2]]
    mz_MS1 <- peaklist[, 8]
    RT_MS1 <- peaklist[, 3]
    Int_MS1 <- peaklist[, 4]
    ############################################################################
    call_DDA_peaklist <- function(i) {
      DDA_list <- NULL
      ##
      x_pl_MS2 <- x_MS2[x_MS2 %in% ((scanNumberStartPL[i] + 1):(scanNumberEndPL[i] - 1))]
      ##
      x_precursor_MS2 <- which(abs(precursorMZ[x_pl_MS2] - mz_MS1[i]) <= massErrorPrecursor)
      L_x_precursor_MS2 <- length(x_precursor_MS2)
      if (L_x_precursor_MS2 > 0) {
        x_precursor <- x_pl_MS2[x_precursor_MS2]
        ##
        if (L_x_precursor_MS2 > 1) {
          ##
          if (DDAspectraIntegrationCheck) {
            ## Integrated DDA spectra
            stackedDDAspectra <- do.call(rbind, lapply(x_precursor, function(j) {
              spectraList[[j]]
            }))
            ##
            DDA_spectra <- spectra_integrator(stackedDDAspectra, massErrorIntegration)
            ##
            DDA_CE <- precursorCE[x_precursor[1]]
            DDA_polarity <- MS_polarity[x_precursor[1]]
            ##
          } else if (DDAionFilteringCheck & L_x_precursor_MS2 >= 3) {
            ## ion filtering
            DDA_spectra <- spectra_ion_filter(spectraList, x_precursor, massErrorIonFiltering, minPercentageDetectedScans, rsdCutoff, pearsonRHOthreshold)
            ##
            DDA_CE <- precursorCE[x_precursor[1]]
            DDA_polarity <- MS_polarity[x_precursor[1]]
            ##
          } else {
            ## The most intense DDA spectra
            x_intensity <- which.max(precursorIntensity[x_precursor])
            scan_max_intensity <- x_precursor[x_intensity[1]]
            DDA_spectra <- spectraList[[scan_max_intensity]]
            ##
            DDA_CE <- precursorCE[scan_max_intensity]
            DDA_polarity <- MS_polarity[scan_max_intensity]
            L_x_precursor_MS2 <- 1
            ####
          }
        } else {
          DDA_spectra <- spectraList[[x_precursor]]
          DDA_CE <- precursorCE[x_precursor]
          DDA_polarity <- MS_polarity[x_precursor]
        }
        ##
        L_DDA_spectra <- dim(DDA_spectra)[1]
        if (L_DDA_spectra > 0) {
          DDA_spectra <- matrix(DDA_spectra[order(DDA_spectra[, 2], decreasing = TRUE), ], ncol = 2)
          DDA_spectra[, 1] <- round(DDA_spectra[, 1], 5)
          ##
          spectralEntropy <- round(spectral_entropy_calculator(DDA_spectra, allowedWeightedSpectralEntropy = TRUE, noiseRemovalRatio = 1e-16)[[1]], 5)
          ##
          if (plotDDAspectraCheck) {
            ##
            pngfilename <- paste0(outputDDAspectra, "/", i, "_peakID", i,"_DDA_spectra_", mz_MS1[i], ".png")
            png(pngfilename, width = 16, height = 8, units = "in", res = 100)
            ##
            plotlabelsMZ <- FSA_annotation_text_repel(FSAspectra = DDA_spectra, nGridX = 8, nGridY = 4)
            ##
            max_int <- max(DDA_spectra[, 2])
            nudge_y = max_int/80
            middle_position <- (DDA_spectra[1, 1] + DDA_spectra[L_DDA_spectra, 1])/2
            ##
            plot(DDA_spectra[, 1], DDA_spectra[, 2], type = "h", ylim = c(0, 1.1*max_int), yaxs = "i", lwd = 2, cex = 4, xlab = "m/z", ylab = "Intensity")
            text(x = DDA_spectra[, 1], y = (DDA_spectra[, 2] + nudge_y), label = plotlabelsMZ, col = "red")
            mtext(DDA_hrms_file, side = 3, adj = 0, line = 0.25, cex = 1.0)
            mtext(text = paste0("S = ", spectralEntropy), side = 3, adj = 1, line = 0.25, cex = 1.0)
            ##
            dev.off()
            ##
          }
          ##
          DDA_list <- cbind(rep(i, L_DDA_spectra), rep(mz_MS1[i], L_DDA_spectra), rep(RT_MS1[i], L_DDA_spectra), rep(Int_MS1[i], L_DDA_spectra), DDA_spectra, rep(spectralEntropy, L_DDA_spectra), rep(DDA_polarity, L_DDA_spectra), rep(DDA_CE, L_DDA_spectra), rep(L_x_precursor_MS2, L_DDA_spectra))
        }
      }
      return(DDA_list)
    }
    ############################################################################
    if (number_processing_threads == 1) {
      DDA_peaklist <- do.call(rbind, lapply(selectedIPApeaks, function(i) {
        call_DDA_peaklist(i)
      }))
    } else {
      osType <- Sys.info()[['sysname']]
      if (osType == "Windows") {
        ##
        clust <- makeCluster(number_processing_threads)
        clusterExport(clust, setdiff(ls(), c("clust", "selectedIPApeaks")), envir = environment())
        ##
        DDA_peaklist <- do.call(rbind, parLapply(clust, selectedIPApeaks, function(i) {
          call_DDA_peaklist(i)
        }))
        ##
        stopCluster(clust)
        ##
      } else {
        ##
        DDA_peaklist <- do.call(rbind, mclapply(selectedIPApeaks, function(i) {
          call_DDA_peaklist(i)
        }, mc.cores = number_processing_threads))
        ##
        closeAllConnections()
      }
    }
    ##
    if (!is.null(DDA_peaklist)) {
      DDA_peaklist <- data.frame(DDA_peaklist)
      DDA_peaklist[, 1] <- as.numeric(DDA_peaklist[, 1])
      DDA_peaklist[, 6] <- round(as.numeric(DDA_peaklist[, 6]), 0)
    } else {
      DDA_peaklist <- data.frame(matrix(rep(0, 10), ncol = 10))
    }
  }
  ##
  rownames(DDA_peaklist) <- NULL
  colnames(DDA_peaklist) <- c("IDSL.IPA_PeakID", "PrecursorMZ", "Precursor_RT", "Precursor_Intensity", "CSA_mz_fragment", "CSA_int_fragment", "Weighted_spectral_entropy_0noiseRemoval", "Ion_mode", "Collision_energy", "Count_DDA_scan")
  ##
  return(DDA_peaklist)
}