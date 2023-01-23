DIA_MS1_fragmentationPeakDetection <- function(DIA_hrms_address, DIA_hrms_file, peaklist, selectedIPApeaks, massError,
                                               smoothingWindowMS1, scanTolerance, nSpline, topRatioPeakHeight,
                                               intensityThresholdFragment, pearsonRHOthreshold, outputDIAeic = NULL,
                                               number_processing_threads = 1) {
  ##
  DIA_peaklist <- data.frame(matrix(rep(0, 10), ncol = 10))
  ##############################################################################
  if (is.null(outputDIAeic)) {
    plotEICcheck <- FALSE
  } else {
    ##
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    ##
    plotEICcheck <- TRUE
    ##
    FSA_dir.create(outputDIAeic, allowedUnlink = FALSE)
  }
  ############################ MS level = 1 ####################################
  p2l <- IDSL.MXP::peak2list(DIA_hrms_address, DIA_hrms_file)
  scanTable <- p2l[["scanTable"]] # this gets table of details for each spectra
  spectraList <- p2l[["spectraList"]] # this gets the spectra values
  p2l <- NULL
  #
  x_MS <- which(scanTable$peaksCount > 0 & scanTable$msLevel == 1) # peaks from soft ionization channel ## some files may not have data in the re-calibration period.
  spectraList <- spectraList[x_MS]
  scanTable <- scanTable[x_MS, ]
  aggregatedSpectraList <- IPA_spectraListAggregator(spectraList)
  #
  RetentionTime <- scanTable$retentionTime
  ##
  precursorCE <- as.matrix(scanTable$collisionEnergy) # Collision energy
  MS1polarity <- as.matrix(scanTable$polarity)
  ##
  n_RT <- length(RetentionTime)     # n_RT is the maximum number of scan number
  ##############################################################################
  scanNumberStartPL <- peaklist[, 1]
  scanNumberEndPL <- peaklist[, 2]
  mz12CIPA <- peaklist[, 8]
  RTIPA <- peaklist[, 3]
  IntIPA <- peaklist[, 4]
  ##############################################################################
  call_DIA_peaklist <- function(i) {
    DIA_list <- NULL
    ##
    scanNumberApex <- which.min(abs(RetentionTime - RTIPA[i]))
    scanNumberStart <- scanNumberStartPL[i] - (scanTolerance + 1)
    if (scanNumberStart < 1) {
      scanNumberStart <- 1
    }
    scanNumberEnd <- scanNumberEndPL[i] + (scanTolerance + 1)
    if (scanNumberEnd > n_RT) {
      scanNumberEnd <- n_RT
    }
    chromatogramMatrixPrecursor <- XIC(aggregatedSpectraList, scanNumberStart, scanNumberEnd, mz12CIPA[i], massError)
    if (!is.null(chromatogramMatrixPrecursor)) {
      chromatogramMatrixPrecursor <- cbind(chromatogramMatrixPrecursor[, 1], chromatogramMatrixPrecursor[, 3], chromatogramMatrixPrecursor[, 3])
      SZC <- nrow(chromatogramMatrixPrecursor)
      ##
      chromatogramMatrixPrecursor <- data.frame(chromatogramMatrixPrecursor)
      colnames(chromatogramMatrixPrecursor) <- c("scanNumber", "smoothChromatogram", "rawChromatogram")
      loess_SZC <- loess(smoothChromatogram ~ scanNumber, data = chromatogramMatrixPrecursor, span = smoothingWindowMS1/SZC, control = loess.control(surface = "direct"))
      chromatogramMatrixPrecursor[, 2] <- predict(loess_SZC)
      x_neg <- which(chromatogramMatrixPrecursor[, 2] < 0)
      chromatogramMatrixPrecursor[x_neg, 2] <- 0
      ##
      if (chromatogramMatrixPrecursor[1, 1] == 1) {
        chromatogramMatrixPrecursor[1, 2] <- 0
      }
      if (chromatogramMatrixPrecursor[SZC, 1] == n_RT) {
        chromatogramMatrixPrecursor[SZC, 2] <- 0
      }
      ##
      Segment <- chromatographicPeakDetector(chromatogramMatrixPrecursor[, 2])
      if (!is.null(Segment)) {
        Segment1 <- Segment + chromatogramMatrixPrecursor[1, 1] - 1
        x_seg_apex <- which(Segment1[, 1] <= scanNumberApex & Segment1[, 2] >= scanNumberApex)
        L_x_seg_apex <- length(x_seg_apex)
        if (L_x_seg_apex > 0) {
          if (L_x_seg_apex > 1) {
            s_x <- do.call(rbind, lapply(x_seg_apex, function(s) {
              x_sh <- which.max(chromatogramMatrixPrecursor[Segment[s, 1]:Segment[s, 2], 3])
              x_sh[1] + Segment[s, 1] - 1
            }))
            x_max <- which.max(chromatogramMatrixPrecursor[s_x, 3])
            x_seg_apex <- x_seg_apex[x_max[1]]
          }
          ##
          chromatogramMatrixPrecursor <- chromatogramMatrixPrecursor[Segment[x_seg_apex, 1]:Segment[x_seg_apex, 2], ]
          ##
          RT_chrom_precursor <- RetentionTime[chromatogramMatrixPrecursor[, 1]]
          Int_chrom_precursor <- chromatogramMatrixPrecursor[, 2]
          ##
          W_precursor <- spline(RT_chrom_precursor, Int_chrom_precursor , n = nSpline, method = "fmm",
                                xmin = RT_chrom_precursor[1], xmax = RT_chrom_precursor[length(RT_chrom_precursor)], ties = mean) # To smooth the curve for derivative calculations
          RT_spline_precursor <- W_precursor[[1]]
          Int_spline_precursor <- W_precursor[[2]]
          #
          x_topRatioPeakHeight <- which(Int_spline_precursor/max(Int_spline_precursor) >= (1 - topRatioPeakHeight))
          RT_spline_precursor <- RT_spline_precursor[x_topRatioPeakHeight]
          Int_spline_precursor <- Int_spline_precursor[x_topRatioPeakHeight]
          ######################################################################
          peaks <- spectraList[[scanNumberApex]]
          # To remove peaks below intensity threshold
          x_mz_fragment <- which((peaks[, 2] >= intensityThresholdFragment) & (abs(peaks[, 1] - mz12CIPA[i]) > massError))
          L_mz_fragment <- length(x_mz_fragment)
          if (L_mz_fragment > 1) {
            mz_fragment <- peaks[x_mz_fragment, 1]
            ##
            DIA_EICs <- lapply(1:L_mz_fragment, function(k) {
              ##
              chromatogramMatrixFragment <- XIC(aggregatedSpectraList, scanNumberStart, scanNumberEnd, mz_fragment[k], massError)
              if (!is.null(chromatogramMatrixFragment)) {
                chromatogramMatrixFragment <- cbind(chromatogramMatrixFragment[, 1], chromatogramMatrixFragment[, 3], chromatogramMatrixFragment[, 3])
                ##
                Top_ScN <- (scanNumberStart - smoothingWindowMS1 - 1):(scanNumberStart - 1)
                x_Top <- which(Top_ScN > 0)
                L_Top <- length(x_Top)
                if (L_Top > 0) {
                  Top_chrom_builder <- cbind(Top_ScN[x_Top], rep(0, L_Top), rep(0, L_Top))
                } else {
                  Top_chrom_builder <- NULL
                }
                Bottom_ScN <- (scanNumberEnd + 1):(scanNumberEnd + smoothingWindowMS1 + 1)
                x_Bottom <- which(Bottom_ScN <= n_RT)
                L_Bottom <- length(x_Bottom)
                if (L_Bottom > 0) {
                  Bottom_chrom_builder <- cbind(Bottom_ScN[x_Bottom], rep(0, L_Bottom), rep(0, L_Bottom))
                } else {
                  Bottom_chrom_builder <- NULL
                }
                chromatogramMatrixFragment <- rbind(Top_chrom_builder, chromatogramMatrixFragment, Bottom_chrom_builder)
                SZC <- nrow(chromatogramMatrixFragment)
                ## Smoothing the chromatogram trace over a smoothing window
                chromatogramMatrixFragment <- data.frame(chromatogramMatrixFragment)
                colnames(chromatogramMatrixFragment) <- c("scan_number", "smooth_chrom", "raw_chrom")
                loess_SZC <- loess(smooth_chrom ~ scan_number, data = chromatogramMatrixFragment, span = smoothingWindowMS1/SZC, control = loess.control(surface = "direct"))
                chromatogramMatrixFragment[, 2] <- predict(loess_SZC)
                x_neg <- which(chromatogramMatrixFragment[, 2] < 0)
                chromatogramMatrixFragment[x_neg, 2] <- 0
                ##
                if (chromatogramMatrixFragment[1, 1] == 1) {
                  chromatogramMatrixFragment[1, 2] <- 0
                }
                if (chromatogramMatrixFragment[SZC, 1] == n_RT) {
                  chromatogramMatrixFragment[SZC, 2] <- 0
                }
                ## Peak detection module
                Segment <- chromatographicPeakDetector(chromatogramMatrixFragment[, 2])
                if (!is.null(Segment)) {
                  Segment2 <- Segment + chromatogramMatrixFragment[1, 1] - 1
                  x_seg_apex <- which(Segment2[, 1] <= scanNumberApex & Segment2[, 2] >= scanNumberApex)
                  L_x_seg_apex <- length(x_seg_apex)
                  if (L_x_seg_apex > 0) {
                    if (L_x_seg_apex > 1) {
                      s_x <- do.call(c, lapply(x_seg_apex, function(s) {
                        x_sh <- which.max(chromatogramMatrixFragment[Segment[s, 1]:Segment[s, 2], 3])
                        x_sh[1] + Segment[s, 1] - 1
                      }))
                      x_max <- which.max(chromatogramMatrixFragment[s_x, 3])
                      x_seg_apex <- x_seg_apex[x_max[1]]
                    }
                    ##
                    chromatogramMatrixFragment <- chromatogramMatrixFragment[Segment[x_seg_apex, 1]:Segment[x_seg_apex, 2], ]
                    ##
                    height_fragment <- max(chromatogramMatrixFragment[, 3]) # raw intensity
                    if (height_fragment > 0) {
                      ##
                      RT_chrom_fragment <- RetentionTime[chromatogramMatrixFragment[, 1]]
                      Int_chrom_fragment <- chromatogramMatrixFragment[, 2] # smooth chromatogram
                      #
                      Int_spline_fragment <- approx(RT_chrom_fragment, Int_chrom_fragment, RT_spline_precursor, method = "linear", 0, 0, rule = 2, f = 0, ties = mean)[[2]]
                      #
                      pearsonRHO <- suppressWarnings(cor(Int_spline_precursor, Int_spline_fragment, method = "pearson"))
                      ##
                      if (!is.na(pearsonRHO)) {
                        if (pearsonRHO >= pearsonRHOthreshold) {
                          if (plotEICcheck) {
                            L_chrom_fragment <- length(Int_chrom_fragment)
                            DIAEICdata <- cbind(rep(k, L_chrom_fragment), RT_chrom_fragment, Int_chrom_fragment)
                          } else {
                            DIAEICdata <- NULL
                          }
                          ##
                          DIA_fragments <- c(mz_fragment[k], height_fragment, pearsonRHO)
                          list(DIAEICdata, DIA_fragments)
                        }
                      }
                    }
                  }
                }
              }
            })
            ####################################################################
            x_fragment_dia <- do.call(c, lapply(1:L_mz_fragment, function(j) {
              if (!is.null(DIA_EICs[[j]])) {
                j
              }
            }))
            ##
            if (length(x_fragment_dia) > 1) {
              ##
              DIA_fragments <- do.call(rbind, lapply(x_fragment_dia, function(j) {
                DIA_EICs[[j]][[2]]
              }))
              ##
              IPA12Cmz <- c(mz12CIPA[i], IntIPA[i], 1)
              DIAfragmentationList <- rbind(IPA12Cmz, DIA_fragments)
              DIAfragmentationList <- DIAfragmentationList[order(DIAfragmentationList[, 2], decreasing = TRUE), ]
              ##
              spectralEntropy <- round(spectral_entropy_calculator(DIAfragmentationList[, 1:2], allowedWeightedSpectralEntropy = TRUE, noiseRemovalRatio = 1e-16)[[1]], 5)
              ##
              L_fragments <- nrow(DIAfragmentationList)
              DIA_list <- cbind(rep(i, L_fragments), rep(mz12CIPA[i], L_fragments), rep(RTIPA[i], L_fragments), rep(IntIPA[i], L_fragments), DIAfragmentationList[, 1:2],
                                rep(spectralEntropy, L_fragments), rep(MS1polarity[scanNumberApex], L_fragments), rep(precursorCE[scanNumberApex], L_fragments), DIAfragmentationList[, 3])
              ##################################################################
              if (plotEICcheck) {
                DIAEICdata <- do.call(rbind, lapply(x_fragment_dia, function(j) {
                  DIA_EICs[[j]][[1]]
                }))
                ##
                yMaxLimPlot <- max(c(Int_chrom_precursor, DIAEICdata[, 3]))
                ##
                xLinesDiff <- c(0, which(abs(diff(DIAEICdata[, 1])) > 0), nrow(DIAEICdata))
                nLines <- length(xLinesDiff)
                nLines1 <- nLines - 1
                ##
                legText <- rep("", nLines)
                legText[1] <- paste0("* m/z = ", round(mz12CIPA[i], 5))
                colors <- c("black", rainbow(nLines1, alpha = 1))
                legCol <- rep("", nLines)
                legCol[1] <- colors[1]
                ##
                orderMSP <- order(DIA_fragments[, 2], decreasing = TRUE)
                orderLines <- do.call(rbind, lapply(2:nLines, function(p) {
                  c((xLinesDiff[p - 1] + 1), xLinesDiff[p])
                }))
                orderLines <- matrix(orderLines[orderMSP, ], ncol = 2)
                ##
                alignedEICfilename <- paste0(outputDIAeic, "/", i, "_MS1_", mz12CIPA[i], "_RT_", RTIPA[i], ".png")
                png(alignedEICfilename, width = 16, height = 8, units = "in", res = 100)
                ##
                par(mar = c(5.1, 4.1, 4.1, 13.8))
                plot(RT_chrom_precursor, Int_chrom_precursor, type = "l", ylim = c(0, yMaxLimPlot*1.01), lwd = 4, col = colors[1], cex = 4, xlab = "", ylab = "")
                ##
                pCounter <- 1
                for (p in 1:nLines1) {
                  pCounter <- pCounter + 1
                  ##
                  xLines <- seq(orderLines[p, 1], orderLines[p, 2], 1)
                  ##
                  lines(DIAEICdata[xLines, 2], DIAEICdata[xLines, 3], lwd = 2, col = colors[pCounter], cex = 4)
                  ##
                  legText[pCounter] <- paste0("m/z = ", round(mz_fragment[DIAEICdata[xLines[1], 1]], 5))
                  legCol[pCounter] <- colors[pCounter]
                }
                ##
                mtext(text = paste0("S = ", spectralEntropy), side = 3, adj = 1, line = 0.25, cex = 1.0)
                mtext("Retention time (min)", side = 1, adj = 0.5, line = 2, cex = 1.35)
                mtext("Intensity", side = 2, adj = 0.5, line = 2, cex = 1.35)
                mtext(DIA_hrms_file, side = 3, adj = 0, line = 0.25, cex = 1.4)
                legend(x = "topright", inset = c(-0.22, 0), legend = legText, lwd = c(4, rep(2, nLines1)), cex = 1.125, bty = "n",
                       col = legCol, seg.len = 1, x.intersp = 0.5, y.intersp = 0.9, xpd = TRUE)
                ##
                dev.off()
                ################################################################
              }
            }
          }
        }
      }
    }
    return(DIA_list)
  }
  ##############################################################################
  if (number_processing_threads == 1) {
    DIA_peaklist <- do.call(rbind, lapply(selectedIPApeaks, function(i) {
      call_DIA_peaklist(i)
    }))
  } else {
    ##
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      DIA_peaklist <- do.call(rbind, mclapply(selectedIPApeaks, function(i) {
        call_DIA_peaklist(i)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      DIA_peaklist <- foreach(i = selectedIPApeaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_DIA_peaklist(i)
      }
      ##
      stopCluster(clust)
      ##
    }
  }
  ##############################################################################
  if (!is.null(DIA_peaklist)) {
    DIA_peaklist <- data.frame(DIA_peaklist)
    DIA_peaklist[, 1] <- as.numeric(DIA_peaklist[, 1])
    DIA_peaklist[, 5] <- round(as.numeric(DIA_peaklist[, 5]), 5)
    DIA_peaklist[, 6] <- round(as.numeric(DIA_peaklist[, 6]), 0)
    DIA_peaklist[, 10] <- round(as.numeric(DIA_peaklist[, 10]), 2)
  }
  ##
  rownames(DIA_peaklist) <- NULL
  colnames(DIA_peaklist) <- c("IDSL.IPA_PeakID", "PrecursorMZ", "Precursor_RT", "Precursor_Intensity", "CSA_mz_fragment", "CSA_int_fragment", "Weighted_spectral_entropy_0noiseRemoval", "Ion_mode", "Collision_energy", "Pearson_rho")
  ##
  return(DIA_peaklist)
}