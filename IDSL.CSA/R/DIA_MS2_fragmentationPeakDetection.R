DIA_MS2_fragmentationPeakDetection <- function(DIA_hrms_address, DIA_hrms_file, peaklist, selectedIPApeaks, massError,
                                               smoothingWindowMS1, smoothingWindowMS2, scanTolerance, nSpline, topRatioPeakHeight,
                                               intensityThresholdFragment, pearsonRHOthreshold, outputDIAeic = NULL, number_processing_threads = 1) {
  ##
  DIA_peaklist <- data.frame(matrix(rep(0, 10), ncol = 10))
  ##############################################################################
  if (is.null(outputDIAeic)) {
    plotEICcheck <- FALSE
  } else {
    ##
    oldpar <- par(no.readonly = TRUE)
    on.exit(suppressWarnings(par(oldpar)))
    ##
    plotEICcheck <- TRUE
    ##
    FSA_dir.create(outputDIAeic, allowedUnlink = FALSE)
  }
  ##############################################################################
  p2l <- IDSL.MXP::peak2list(DIA_hrms_address, DIA_hrms_file)
  scanTable <- p2l[["scanTable"]] # this gets table of details for each spectra
  spectraList <- p2l[["spectraList"]] # this gets the spectra values
  p2l <- NULL
  ##
  x_MS <- which(scanTable$peaksCount > 0) # peaks from soft ionization channel ## some files may not have data in the column re-calibration period.
  spectraList <- spectraList[x_MS]
  scanTable <- scanTable[x_MS, ]  
  ##
  RetentionTime <- scanTable$retentionTime
  ##
  x_MS1 <- which(scanTable$msLevel == 1)
  x_MS2 <- which(scanTable$msLevel == 2)
  n_RT_MS2 <- length(x_MS2)   # n_RT_MS2 is the maximum number of scan number
  if (n_RT_MS2 > 0) {
    precursorCE <- as.matrix(scanTable$collisionEnergy) # Collision energy
    ############################ MS level = 1 ##################################
    spectraList_MS1 <- spectraList[x_MS1]
    aggregatedSpectraListMS1 <- IPA_spectraListAggregator(spectraList_MS1)
    spectraList_MS1 <- NULL
    RetentionTime_MS1 <- RetentionTime[x_MS1]
    n_RT_MS1 <- length(x_MS1)   # n_RT_MS1 is the maximum number of scan number
    ############################ MS level = 2 ##################################
    spectraList_MS2 <- spectraList[x_MS2]
    aggregatedSpectraListMS2 <- IPA_spectraListAggregator(spectraList_MS2)
    spectraList <- NULL
    RetentionTime_MS2 <- RetentionTime[x_MS2]
    precursorCE_MS2 <- precursorCE[x_MS2]
    MS2polarity <- as.numeric(scanTable$polarity[x_MS2])
    ##
    isolationWindowTargetMZ <- as.numeric(scanTable$isolationWindowTargetMZ)
    isolationWindowLowerOffset <- as.numeric(scanTable$isolationWindowLowerOffset)
    isolationWindowUpperOffset <- as.numeric(scanTable$isolationWindowUpperOffset)
    if (length(!is.na(isolationWindowTargetMZ)) == 0 | length(!is.na(isolationWindowLowerOffset)) == 0 | length(!is.na(isolationWindowUpperOffset)) == 0) {
      isolationWindowOffsetCheck <- FALSE
    } else {
      isolationWindowOffsetCheck <- TRUE
      ##
      isolationWindowLowerLimit <- isolationWindowTargetMZ - isolationWindowLowerOffset
      isolationWindowUpperLimit <- isolationWindowTargetMZ + isolationWindowUpperOffset
    }
    ##
    scanTable <- NULL
    ############################################################################
    scanNumberStartPL <- peaklist[, 1]
    scanNumberEndPL <- peaklist[, 2]
    mz12CIPA <- peaklist[, 8]
    RTIPA <- peaklist[, 3]
    IntIPA <- peaklist[, 4]
    ############################################################################
    call_DIA_peaklist <- function(i) {
      DIA_list <- NULL
      ##
      if (isolationWindowOffsetCheck) {
        x_pl_MS2 <- x_MS2[x_MS2 %in% ((x_MS1[scanNumberStartPL[i]] + 1):(x_MS1[scanNumberEndPL[i]] - 1))]
        ##
        x_precursor_MS2 <- which(mz12CIPA[i] >= isolationWindowLowerLimit[x_pl_MS2] & mz12CIPA[i] <= isolationWindowUpperLimit[x_pl_MS2])
        if (length(x_precursor_MS2) > 2) {
          DIA2check <- TRUE
        } else {
          DIA2check <- FALSE
        }
      } else {
        DIA2check <- TRUE
      }
      ##########################################################################
      if (DIA2check) {
        ##
        scanNumberApex1 <- which.min(abs(RetentionTime_MS1 - RTIPA[i]))
        scanNumberStart1 <- scanNumberStartPL[i] - (scanTolerance + 1)
        if (scanNumberStart1 < 1) {
          scanNumberStart1 <- 1
        }
        scanNumberEnd1 <- scanNumberEndPL[i] + (scanTolerance + 1)
        if (scanNumberEnd1 > n_RT_MS1) {
          scanNumberEnd1 <- n_RT_MS1
        }
        chromatogramMatrixMS1 <- XIC(aggregatedSpectraListMS1, scanNumberStart1, scanNumberEnd1, mz12CIPA[i], massError)
        if (!is.null(chromatogramMatrixMS1)) {
          chromatogramMatrixMS1 <- cbind(chromatogramMatrixMS1[, 1], chromatogramMatrixMS1[, 3], chromatogramMatrixMS1[, 3])
          SZC <- nrow(chromatogramMatrixMS1)
          ##
          chromatogramMatrixMS1 <- data.frame(chromatogramMatrixMS1)
          colnames(chromatogramMatrixMS1) <- c("scanNumber", "smoothChromatogram", "rawChromatogram")
          loess_SZC <- loess(smoothChromatogram ~ scanNumber, data = chromatogramMatrixMS1, span = smoothingWindowMS1/SZC, control = loess.control(surface = "direct"))
          chromatogramMatrixMS1[, 2] <- predict(loess_SZC)
          x_neg <- which(chromatogramMatrixMS1[, 2] < 0)
          chromatogramMatrixMS1[x_neg, 2] <- 0
          ##
          if (chromatogramMatrixMS1[1, 1] == 1) {
            chromatogramMatrixMS1[1, 2] <- 0
          }
          if (chromatogramMatrixMS1[SZC, 1] == n_RT_MS1) {
            chromatogramMatrixMS1[SZC, 2] <- 0
          }
          ##
          Segment <- chromatographicPeakDetector(chromatogramMatrixMS1[, 2])
          if (!is.null(Segment)) {
            Segment1 <- Segment + chromatogramMatrixMS1[1, 1] - 1
            x_seg_apex <- which(Segment1[, 1] <= scanNumberApex1 & Segment1[, 2] >= scanNumberApex1)
            L_x_seg_apex <- length(x_seg_apex)
            if (L_x_seg_apex > 0) {
              if (L_x_seg_apex > 1) {
                s_x <- do.call(c, lapply(x_seg_apex, function(s) {
                  x_sh <- which.max(chromatogramMatrixMS1[Segment[s, 1]:Segment[s, 2], 3])
                  x_sh[1] + Segment[s, 1] - 1
                }))
                x_max <- which.max(chromatogramMatrixMS1[s_x, 3])
                x_seg_apex <- x_seg_apex[x_max[1]]
              }
              ##
              chromatogramMatrixMS1 <- chromatogramMatrixMS1[Segment[x_seg_apex, 1]:Segment[x_seg_apex, 2], ]
              ##
              RT_chrom_MS1 <- RetentionTime_MS1[chromatogramMatrixMS1[, 1]]
              Int_chrom_MS1 <- chromatogramMatrixMS1[, 2]
              ##
              W_MS1 <- spline(RT_chrom_MS1, Int_chrom_MS1 , n = nSpline, method = "fmm",
                              xmin = RT_chrom_MS1[1], xmax = RT_chrom_MS1[length(RT_chrom_MS1)], ties = mean) # To smooth the curve for derivative calculations
              RT_spline_MS1 <- W_MS1[[1]]
              Int_spline_MS1 <- W_MS1[[2]]
              #
              x_topRatioPeakHeight <- which(Int_spline_MS1/max(Int_spline_MS1) >= (1 - topRatioPeakHeight))
              RT_spline_MS1 <- RT_spline_MS1[x_topRatioPeakHeight]
              Int_spline_MS1 <- Int_spline_MS1[x_topRatioPeakHeight]
              ######################### MS level = 2 ###########################
              if (isolationWindowOffsetCheck) {
                ## update the `x_precursor_MS2` parameter
                x_precursor_MS2 <- x_pl_MS2[x_precursor_MS2]
                x_precursor_MS2 <- which(x_MS2 %in% x_precursor_MS2)
                ##
                scanNumberApex2 <- which.min(abs(RetentionTime_MS2[x_precursor_MS2] - RTIPA[i]))
                scanNumberApex2 <- x_precursor_MS2[scanNumberApex2]
              } else {
                scanNumberApex2 <- which.min(abs(RetentionTime_MS2 - RTIPA[i]))
              }
              x_scanNumberStart2 <- which.min(abs(x_MS1[scanNumberStartPL[i]] - x_MS2))[1]
              scanNumberStart2 <- x_scanNumberStart2 - (scanTolerance + 1)
              if (scanNumberStart2 < 1) {
                scanNumberStart2 <- 1
              }
              ##
              x_scanNumberEnd2 <- which.min(abs(x_MS1[scanNumberEndPL[i]] - x_MS2))[1]
              scanNumberEnd2 <- x_scanNumberEnd2 + (scanTolerance + 1)
              if (scanNumberEnd2 > n_RT_MS2) {
                scanNumberEnd2 <- n_RT_MS2
              }
              ##################################################################
              peaks_MS2 <- spectraList_MS2[[scanNumberApex2]]
              ##
              x_mz_MS2 <- which((peaks_MS2[, 2] >= intensityThresholdFragment) & (peaks_MS2[, 1] <= (mz12CIPA[i] + 5))) # 5 was added to include the isotope envelope of the precursor mass
              L_mz_MS2 <- length(x_mz_MS2)
              if (L_mz_MS2 > 1) {
                mz_MS2 <- peaks_MS2[x_mz_MS2, 1]
                ##
                DIA_EICs <- lapply(1:L_mz_MS2, function(k) {
                  ##
                  chromatogramMatrixMS2 <- XIC(aggregatedSpectraListMS2, scanNumberStart2, scanNumberEnd2, mz_MS2[k], massError)
                  ##
                  if (!is.null(chromatogramMatrixMS2)) {
                    if (isolationWindowOffsetCheck) {
                      chromatogramMatrixMS2[which(!(chromatogramMatrixMS2[, 1] %in% x_precursor_MS2)), c(2, 3)] <- 0
                      if (length(which(chromatogramMatrixMS2[, 3] == 0)) < 3) {
                        chromatogramMatrixMS2 <- NULL
                      }
                    }
                    ##
                    if (!is.null(chromatogramMatrixMS2)) {
                      chromatogramMatrixMS2 <- cbind(chromatogramMatrixMS2[, 1], chromatogramMatrixMS2[, 3], chromatogramMatrixMS2[, 3])
                      ##
                      Top_ScN <- (scanNumberStart2 - smoothingWindowMS2 - 1):(scanNumberStart2 - 1)
                      x_Top <- which(Top_ScN > 0)
                      L_Top <- length(x_Top)
                      if (L_Top > 0) {
                        Top_chrom_builder <- cbind(Top_ScN[x_Top], rep(0, L_Top), rep(0, L_Top))
                      } else {
                        Top_chrom_builder <- NULL
                      }
                      Bottom_ScN <- (scanNumberEnd2 + 1):(scanNumberEnd2 + smoothingWindowMS2 + 1)
                      x_Bottom <- which(Bottom_ScN <= n_RT_MS2)
                      L_Bottom <- length(x_Bottom)
                      if (L_Bottom > 0) {
                        Bottom_chrom_builder <- cbind(Bottom_ScN[x_Bottom], rep(0, L_Bottom), rep(0, L_Bottom))
                      } else {
                        Bottom_chrom_builder <- NULL
                      }
                      chromatogramMatrixMS2 <- rbind(Top_chrom_builder, chromatogramMatrixMS2, Bottom_chrom_builder)
                      SZC <- nrow(chromatogramMatrixMS2)
                      ## Smoothing the chromatogram trace over a smoothing window
                      chromatogramMatrixMS2 <- data.frame(chromatogramMatrixMS2)
                      colnames(chromatogramMatrixMS2) <- c("scan_number", "smooth_chrom", "raw_chrom")
                      loess_SZC <- loess(smooth_chrom ~ scan_number, data = chromatogramMatrixMS2, span = smoothingWindowMS2/SZC, control = loess.control(surface = "direct"))
                      chromatogramMatrixMS2[, 2] <- predict(loess_SZC)
                      x_neg <- which(chromatogramMatrixMS2[, 2] < 0)
                      chromatogramMatrixMS2[x_neg, 2] <- 0
                      ##
                      if (chromatogramMatrixMS2[1, 1] == 1) {
                        chromatogramMatrixMS2[1, 2] <- 0
                      }
                      if (chromatogramMatrixMS2[SZC, 1] == n_RT_MS2) {
                        chromatogramMatrixMS2[SZC, 2] <- 0
                      }
                      ## Peak detection module
                      Segment <- chromatographicPeakDetector(chromatogramMatrixMS2[, 2])
                      if (!is.null(Segment)) {
                        Segment2 <- Segment + chromatogramMatrixMS2[1, 1] - 1
                        x_seg_apex <- which(Segment2[, 1] <= scanNumberApex2 & Segment2[, 2] >= scanNumberApex2)
                        L_x_seg_apex <- length(x_seg_apex)
                        if (L_x_seg_apex > 0) {
                          if (L_x_seg_apex > 1) {
                            s_x <- do.call(c, lapply(x_seg_apex, function(s) {
                              x_sh <- which.max(chromatogramMatrixMS2[Segment[s, 1]:Segment[s, 2], 3])
                              x_sh[1] + Segment[s, 1] - 1
                            }))
                            x_max <- which.max(chromatogramMatrixMS2[s_x, 3])
                            x_seg_apex <- x_seg_apex[x_max[1]]
                          }
                          ##
                          chromatogramMatrixMS2 <- chromatogramMatrixMS2[Segment[x_seg_apex, 1]:Segment[x_seg_apex, 2], ]
                          ##
                          MS2_height <- max(chromatogramMatrixMS2[, 3]) # raw intensity
                          ##
                          if (MS2_height > 0) {
                            ##
                            RT_chrom_MS2 <- RetentionTime_MS2[chromatogramMatrixMS2[, 1]]
                            Int_chrom_MS2 <- chromatogramMatrixMS2[, 2] # smooth chromatogram
                            #
                            Int_spline_MS2 <- approx(RT_chrom_MS2, Int_chrom_MS2, RT_spline_MS1, method = "linear", 0, 0, rule = 2, f = 0, ties = mean)[[2]]
                            #
                            pearsonRHO <- suppressWarnings(cor(Int_spline_MS1, Int_spline_MS2, method = "pearson"))
                            ##
                            if (!is.na(pearsonRHO)) {
                              if (pearsonRHO >= pearsonRHOthreshold) {
                                if (plotEICcheck) {
                                  L_chrom_MS2 <- length(Int_chrom_MS2)
                                  DIAEICdata <- cbind(rep(k, L_chrom_MS2), RT_chrom_MS2, Int_chrom_MS2)
                                } else {
                                  DIAEICdata <- NULL
                                }
                                ##
                                xTopRatioPeakHeight <- which(chromatogramMatrixMS2[, 3]/MS2_height >= (1 - topRatioPeakHeight))
                                MS2_height <- sum(chromatogramMatrixMS2[xTopRatioPeakHeight, 3]) # to use an integrated intensity of the raw chromatogram
                                ##
                                DIA_fragments <- c(mz_MS2[k], MS2_height, pearsonRHO)
                                list(DIAEICdata, DIA_fragments)
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                })
                ################################################################
                x_fragment_dia <- do.call(c, lapply(1:L_mz_MS2, function(j) {
                  if (!is.null(DIA_EICs[[j]])) {
                    j
                  }
                }))
                ##
                if (length(x_fragment_dia) > 1) {
                  ##
                  DIAfragmentationList <- do.call(rbind, lapply(x_fragment_dia, function(j) {
                    DIA_EICs[[j]][[2]]
                  }))
                  orderMSP <- order(DIAfragmentationList[, 2], decreasing = TRUE)
                  DIAfragmentationList <- matrix(DIAfragmentationList[orderMSP, ], ncol = 3)
                  ##
                  spectralEntropy <- round(spectral_entropy_calculator(DIAfragmentationList[, 1:2], allowedWeightedSpectralEntropy = TRUE, noiseRemovalRatio = 1e-16)[[1]], 5)
                  ##
                  L_fragments <- nrow(DIAfragmentationList)
                  DIA_list <- cbind(rep(i, L_fragments), rep(mz12CIPA[i], L_fragments), rep(RTIPA[i], L_fragments), rep(IntIPA[i], L_fragments), DIAfragmentationList[, 1:2],
                                    rep(spectralEntropy, L_fragments), rep(MS2polarity[scanNumberApex2], L_fragments), rep(precursorCE_MS2[scanNumberApex2], L_fragments), DIAfragmentationList[, 3])
                  ##############################################################
                  if (plotEICcheck) {
                    DIAEICdata <- do.call(rbind, lapply(x_fragment_dia, function(j) {
                      DIA_EICs[[j]][[1]]
                    }))
                    ##
                    yMaxLimPlot <- max(c(Int_chrom_MS1, DIAEICdata[, 3]))
                    ##
                    xLinesDiff <- c(0, which(abs(diff(DIAEICdata[, 1])) > 0), nrow(DIAEICdata))
                    nLines <- length(xLinesDiff)
                    ##
                    legText <- rep("", nLines)
                    legText[1] <- paste0("(MS1) m/z = ", round(mz12CIPA[i], 5))
                    colors <- c("black", rainbow(nLines, alpha = 1))
                    legCol <- rep("", nLines)
                    legCol[1] <- colors[1]
                    ##
                    orderLines <- do.call(rbind, lapply(2:nLines, function(p) {
                      c((xLinesDiff[p - 1] + 1), xLinesDiff[p])
                    }))
                    orderLines <- matrix(orderLines[orderMSP, ], ncol = 2)
                    ##
                    alignedEICfilename <- paste0(outputDIAeic, "/", i, "_MS1_", mz12CIPA[i], "_RT_", RTIPA[i], ".png")
                    png(alignedEICfilename, width = 16, height = 8, units = "in", res = 100)
                    ##
                    par(mar = c(5.1, 4.1, 4.1, 13.8))
                    plot(RT_chrom_MS1, Int_chrom_MS1, type = "l", ylim = c(0, yMaxLimPlot*1.01), lwd = 4, col = colors[1], cex = 4, yaxt = "n", xlab = "", ylab = "")
                    ##
                    pCounter <- 1
                    for (p in 1:(nLines - 1)) {
                      pCounter <- pCounter + 1
                      ##
                      xLines <- seq(orderLines[p, 1], orderLines[p, 2], 1)
                      ##
                      lines(DIAEICdata[xLines, 2], DIAEICdata[xLines, 3], lwd = 2, col = colors[pCounter], cex = 4)
                      ##
                      legText[pCounter] <- paste0("(MS2) m/z = ", round(mz_MS2[DIAEICdata[xLines[1], 1]], 5))
                      legCol[pCounter] <- colors[pCounter]
                    }
                    ##
                    mtext(text = paste0("S = ", spectralEntropy), side = 3, adj = 1, line = 0.25, cex = 1.0)
                    mtext("Retention time (min)", side = 1, adj = 0.5, line = 2, cex = 1.35)
                    mtext("Intensity", side = 2, adj = 0.50, line = 1, cex = 1.35)
                    mtext(DIA_hrms_file, side = 3, adj = 0, line = 0.25, cex = 1.0)
                    legend(x = "topright", inset = c(-0.22, 0), legend = legText, lwd = c(4, rep(2, nLines)), cex = 1.125, bty = "n",
                           col = legCol, seg.len = 1, x.intersp = 0.5, y.intersp = 0.9, xpd = TRUE)
                    ##
                    dev.off()
                    ############################################################
                  }
                }
              }
            }
          }
        }
      }
      return(DIA_list)
    }
    ############################################################################
    if (number_processing_threads == 1) {
      DIA_peaklist <- do.call(rbind, lapply(selectedIPApeaks, function(i) {
        call_DIA_peaklist(i)
      }))
    } else {
      osType <- Sys.info()[['sysname']]
      if (osType == "Windows") {
        ##
        clust <- makeCluster(number_processing_threads)
        clusterExport(clust, setdiff(ls(), c("clust", "selectedIPApeaks")), envir = environment())
        ##
        DIA_peaklist <- do.call(rbind, parLapply(clust, selectedIPApeaks, function(i) {
          call_DIA_peaklist(i)
        }))
        ##
        stopCluster(clust)
        ##
      } else {
        ##
        DIA_peaklist <- do.call(rbind, mclapply(selectedIPApeaks, function(i) {
          call_DIA_peaklist(i)
        }, mc.cores = number_processing_threads))
        ##
        closeAllConnections()
        ##
      }
    }
    ##
    if (!is.null(DIA_peaklist)) {
      DIA_peaklist <- data.frame(DIA_peaklist)
      DIA_peaklist[, 1] <- as.numeric(DIA_peaklist[, 1])
      DIA_peaklist[, 5] <- round(as.numeric(DIA_peaklist[, 5]), 5)
      DIA_peaklist[, 6] <- round(as.numeric(DIA_peaklist[, 6]), 0)
      DIA_peaklist[, 10] <- round(as.numeric(DIA_peaklist[, 10]), 2)
    } else {
      DIA_peaklist <- data.frame(matrix(rep(0, 10), ncol = 10))
    }
  }
  ##
  rownames(DIA_peaklist) <- NULL
  colnames(DIA_peaklist) <- c("IDSL.IPA_PeakID", "PrecursorMZ", "Precursor_RT", "Precursor_Intensity", "CSA_mz_fragment", "CSA_int_fragment", "Weighted_spectral_entropy_0noiseRemoval", "Ion_mode", "Collision_energy", "Pearson_rho")
  ##
  return(DIA_peaklist)
}