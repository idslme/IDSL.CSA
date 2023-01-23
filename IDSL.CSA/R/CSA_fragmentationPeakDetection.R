CSA_fragmentationPeakDetection <- function(CSA_hrms_address, CSA_hrms_file, tempAlignedTableSubsetsFolder = NULL,
                                           peaklist, selectedIPApeaks = NULL, RTtolerance, massError, minSNRbaseline,
                                           smoothingWindowMS1, scanTolerance, nSpline, topRatioPeakHeight,
                                           minIonRangeDifference, minNumCSApeaks, pearsonRHOthreshold, outputCSAeic = NULL) {
  ##
  minNumCSApeaks <- minNumCSApeaks - 1 # To count for the seed ion
  ##
  ##############################################################################
  ##
  if (is.null(outputCSAeic)) {
    plotEICcheck <- FALSE
  } else {
    ##
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    ##
    plotEICcheck <- TRUE
    ##
    FSA_dir.create(outputCSAeic, allowedUnlink = FALSE)
  }
  ##
  ################################ MS level = 1 ################################
  ##
  p2l <- IDSL.MXP::peak2list(CSA_hrms_address, CSA_hrms_file)
  scanTable <- p2l[["scanTable"]] # this gets table of details for each spectra
  spectraList <- p2l[["spectraList"]] # this gets the spectra values
  p2l <- NULL
  ##
  xMS1 <- which(scanTable$peaksCount > 0 & scanTable$msLevel == 1) # peaks from soft ionization channel ## some files may not have data in the re-calibration column.
  spectraList <- spectraList[xMS1]
  scanTable <- scanTable[xMS1, ]
  aggregatedSpectraList <- IPA_spectraListAggregator(spectraList)
  spectraList <- NULL
  ##
  retentionTime <- scanTable$retentionTime
  ##
  precursorCE <- as.matrix(scanTable$collisionEnergy) # Collision energy
  MS1polarity <- as.matrix(scanTable$polarity)
  ##
  LretentionTime <- length(retentionTime)
  ##
  ##############################################################################
  ##
  if (!is.null(tempAlignedTableSubsetsFolder)) {
    alignedTableCheck <- TRUE
    ##
    xNon0Xcol <- IDSL.IPA::loadRdata(paste0(tempAlignedTableSubsetsFolder, "/", CSA_hrms_file, "/xNon0Xcol.Rdata"))
    peakXcolID <- IDSL.IPA::loadRdata(paste0(tempAlignedTableSubsetsFolder, "/", CSA_hrms_file, "/peakXcolSubset.Rdata"))
    correlationList <- IDSL.IPA::loadRdata(paste0(tempAlignedTableSubsetsFolder, "/", CSA_hrms_file, "/correlationListSubset.Rdata"))
    ##
    IPApeakID <- peakXcolID[xNon0Xcol]
    orderIPApeakID <- order(peaklist[IPApeakID, 4], decreasing = TRUE)
    correlationList <- correlationList[orderIPApeakID]
    IPApeakID <- IPApeakID[orderIPApeakID]
    nPeaks <- length(IPApeakID)
  } else {
    alignedTableCheck <- FALSE
    ##
    IPApeakID <- order(peaklist[, 4], decreasing = TRUE)
    nPeaks <- dim(peaklist)[1]
  }
  scanNumberStartPL <- peaklist[IPApeakID, 1]
  scanNumberEndPL <- peaklist[IPApeakID, 2]
  mz12CIPA <- peaklist[IPApeakID, 8]
  RTIPA <- peaklist[IPApeakID, 3]
  IntIPA <- peaklist[IPApeakID, 4]
  mz13CIPA <- peaklist[IPApeakID, 10]
  SNRbaseline <- peaklist[IPApeakID, 21]
  ##############################################################################
  i <- 1
  counterCSAblock <- 0
  listCSAalignedTable <- vector(mode = "list", nPeaks)
  while (i < nPeaks) {
    if ((IPApeakID[i] != 0) & (SNRbaseline[i] >= minSNRbaseline)) {
      ##
      if (alignedTableCheck) {
        listCL <- correlationList[[i]]
        listCL <- peakXcolID[listCL]
        listCL <- listCL[listCL != 0]
        xCL <- do.call(c, lapply(listCL, function(j) {
          which(IPApeakID == j)
        }))
        ##
        if (length(xCL) > 0) {
          xPLretentionTime <- which((abs(RTIPA[xCL] - RTIPA[i]) <= RTtolerance) & (IPApeakID[xCL] != 0))
          LxPLretentionTime <- length(xPLretentionTime)
          if (LxPLretentionTime > 0) {
            xPLretentionTime <- xCL[xPLretentionTime]
          } else {
            LxPLretentionTime <- 0
          }
        } else {
          LxPLretentionTime <- 0
        }
      } else {
        xPLretentionTime <- which((abs(RTIPA - RTIPA[i]) <= RTtolerance) & (IPApeakID != 0))
        xPLretentionTime <- setdiff(xPLretentionTime, i)
        LxPLretentionTime <- length(xPLretentionTime)
      }
      ##
      if (LxPLretentionTime >= minNumCSApeaks) {
        mz_fragment <- c(mz12CIPA[xPLretentionTime], mz13CIPA[i], mz13CIPA[xPLretentionTime])
        peakID <- c(xPLretentionTime, rep(0, (LxPLretentionTime + 1)))
        ########################################################################
        nCSAfragmentIons <- length(mz_fragment)
        if (nCSAfragmentIons >= minNumCSApeaks) {
          ######################################################################
          iMzCSAions <- c(mz12CIPA[i], mz_fragment)
          ionRangeDifference <- max(iMzCSAions) - min(iMzCSAions)
          if (ionRangeDifference >= minIonRangeDifference) {
            iMzCSA12Cions <- mz12CIPA[c(i, xPLretentionTime)]
            ##
            nMzCSA12Cions <- length(which(diff(iMzCSA12Cions[order(iMzCSA12Cions, decreasing = FALSE)]) >= 2.006709670672))
            if (nMzCSA12Cions >= minNumCSApeaks) {
              ##################################################################
              scanNumberApex <- which.min(abs(retentionTime - RTIPA[i]))
              scanNumberStart <- scanNumberStartPL[i] - (scanTolerance + 1)
              if (scanNumberStart < 1) {
                scanNumberStart <- 1
              }
              scanNumberEnd <- scanNumberEndPL[i] + (scanTolerance + 1)
              if (scanNumberEnd > LretentionTime) {
                scanNumberEnd <- LretentionTime
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
                xNeg <- which(chromatogramMatrixPrecursor[, 2] < 0)
                chromatogramMatrixPrecursor[xNeg, 2] <- 0
                ##
                if (chromatogramMatrixPrecursor[1, 1] == 1) {
                  chromatogramMatrixPrecursor[1, 2] <- 0
                }
                if (chromatogramMatrixPrecursor[SZC, 1] == LretentionTime) {
                  chromatogramMatrixPrecursor[SZC, 2] <- 0
                }
                ##
                Segment <- chromatographicPeakDetector(chromatogramMatrixPrecursor[, 2])
                if (!is.null(Segment)) {
                  Segment1 <- Segment + chromatogramMatrixPrecursor[1, 1] - 1
                  xSegApex <- which(Segment1[, 1] <= scanNumberApex & Segment1[, 2] >= scanNumberApex)
                  LxSegApex <- length(xSegApex)
                  if (LxSegApex > 0) {
                    if (LxSegApex > 1) {
                      sX <- do.call(c, lapply(xSegApex, function(s) {
                        xSH <- which.max(chromatogramMatrixPrecursor[Segment[s, 1]:Segment[s, 2], 3])
                        xSH[1] + Segment[s, 1] - 1
                      }))
                      xMax <- which.max(chromatogramMatrixPrecursor[sX, 3])
                      xSegApex <- xSegApex[xMax[1]]
                    }
                    ##
                    chromatogramMatrixPrecursor <- chromatogramMatrixPrecursor[Segment[xSegApex, 1]:Segment[xSegApex, 2], ]
                    ##
                    RT_chrom_precursor <- retentionTime[chromatogramMatrixPrecursor[, 1]]
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
                    ############################################################
                    CSA_EICs <- vector(mode = "list", nCSAfragmentIons)
                    indexCSAfragmentIons <- rep(0, nCSAfragmentIons)
                    numCSApeaks <- 0
                    ##
                    L_12 <- (nCSAfragmentIons + 1)/2
                    k <- 1
                    mz13CIPAcheck <- TRUE
                    minNumCSApeaksPass <- FALSE
                    ##
                    while (k <= nCSAfragmentIons) {
                      ##
                      if (mz13CIPAcheck) {
                        if (k == L_12) {
                          mz13CIPAcheck <- FALSE
                          ##
                          if (numCSApeaks >= minNumCSApeaks) {
                            ##
                            index12C <- c(1, (indexCSAfragmentIons[1:numCSApeaks] + 1))
                            orderMZ12C <- order(iMzCSA12Cions[index12C], decreasing = FALSE)
                            index12C <- index12C[orderMZ12C]
                            nMzCSA12Cions <- length(which(diff(iMzCSA12Cions[index12C]) >= 2.006709670672))
                            if (nMzCSA12Cions >= minNumCSApeaks) {
                              minNumCSApeaksPass <- TRUE
                            }
                          }
                          ##
                          if (!minNumCSApeaksPass) {
                            k <- nCSAfragmentIons
                            peakID[k] <- NA
                          }
                        }
                      }
                      ##
                      if (!is.na(peakID[k])) {
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
                          x_Bottom <- which(Bottom_ScN <= LretentionTime)
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
                          colnames(chromatogramMatrixFragment) <- c("scanNumber", "smoothChromatogram", "rawChromatogram")
                          loess_SZC <- loess(smoothChromatogram ~ scanNumber, data = chromatogramMatrixFragment, span = smoothingWindowMS1/SZC, control = loess.control(surface = "direct"))
                          chromatogramMatrixFragment[, 2] <- predict(loess_SZC)
                          xNeg <- which(chromatogramMatrixFragment[, 2] < 0)
                          chromatogramMatrixFragment[xNeg, 2] <- 0
                          ##
                          if (chromatogramMatrixFragment[1, 1] == 1) {
                            chromatogramMatrixFragment[1, 2] <- 0
                          }
                          if (chromatogramMatrixFragment[SZC, 1] == LretentionTime) {
                            chromatogramMatrixFragment[SZC, 2] <- 0
                          }
                          ## Peak detection module
                          Segment <- chromatographicPeakDetector(chromatogramMatrixFragment[, 2])
                          if (!is.null(Segment)) {
                            Segment2 <- Segment + chromatogramMatrixFragment[1, 1] - 1
                            xSegApex <- which(Segment2[, 1] <= scanNumberApex & Segment2[, 2] >= scanNumberApex)
                            LxSegApex <- length(xSegApex)
                            if (LxSegApex > 0) {
                              if (LxSegApex > 1) {
                                sX <- do.call(c, lapply(xSegApex, function(s) {
                                  xSH <- which.max(chromatogramMatrixFragment[Segment[s, 1]:Segment[s, 2], 3])
                                  xSH[1] + Segment[s, 1] - 1
                                }))
                                xMax <- which.max(chromatogramMatrixFragment[sX, 3])
                                xSegApex <- xSegApex[xMax[1]]
                              }
                              ##
                              chromatogramMatrixFragment <- chromatogramMatrixFragment[Segment[xSegApex, 1]:Segment[xSegApex, 2], ]
                              ##
                              height_fragment <- max(chromatogramMatrixFragment[, 3]) # raw intensity
                              if (height_fragment > 0) {
                                ##
                                RT_chrom_fragment <- retentionTime[chromatogramMatrixFragment[, 1]]
                                Int_chrom_fragment <- chromatogramMatrixFragment[, 2] # smooth chromatogram
                                ##
                                Int_spline_fragment <- approx(RT_chrom_fragment, Int_chrom_fragment, RT_spline_precursor, method = "linear", 0, 0, rule = 2, f = 0, ties = mean)[[2]]
                                ##
                                pearsonRHO <- suppressWarnings(cor(Int_spline_precursor, Int_spline_fragment, method = "pearson"))
                                ##
                                if (!is.na(pearsonRHO)) {
                                  if (pearsonRHO >= pearsonRHOthreshold) {
                                    if (plotEICcheck) {
                                      L_chrom_fragment <- length(Int_chrom_fragment)
                                      CSAEICdata <- cbind(rep(k, L_chrom_fragment), RT_chrom_fragment, Int_chrom_fragment)
                                    } else {
                                      CSAEICdata <- NULL
                                    }
                                    ##
                                    CSA_fragments <- c(mz_fragment[k], height_fragment, pearsonRHO)
                                    CSA_EICs[[k]] <- list(CSAEICdata, CSA_fragments, peakID[k])
                                  }
                                  ##
                                  if (is.null(CSA_EICs[[k]])) {
                                    if (k < L_12) {
                                      k1 <- k + L_12
                                      peakID[k1] <- NA
                                    }
                                  } else {
                                    numCSApeaks <- numCSApeaks + 1
                                    indexCSAfragmentIons[numCSApeaks] <- k
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                      k <- k + 1
                    }
                    ############################################################
                    if (minNumCSApeaksPass) {
                      indexCSAfragmentIons <- indexCSAfragmentIons[indexCSAfragmentIons != 0]
                      ##
                      CSA_fragments <- do.call(rbind, lapply(indexCSAfragmentIons, function(j) {
                        CSA_EICs[[j]][[2]]
                      }))
                      ##
                      IPA12Cmz <- c(mz12CIPA[i], IntIPA[i], 1)
                      CSAfragmentationList <- rbind(IPA12Cmz, CSA_fragments)
                      ionRangeDifference <- max(CSAfragmentationList[, 1]) - min(CSAfragmentationList[, 1])
                      ##
                      if (ionRangeDifference >= minIonRangeDifference) {
                        ########################################################
                        jIPApeakID <- do.call(c, lapply(indexCSAfragmentIons, function(j) {
                          if (CSA_EICs[[j]][[3]] != 0) {
                            IPApeakID[CSA_EICs[[j]][[3]]]
                          } else {
                            0
                          }
                        }))
                        jIPApeakID <- c(IPApeakID[i], jIPApeakID)
                        CSAfragmentationList <- cbind(CSAfragmentationList, jIPApeakID)
                        ##
                        CSAfragmentationList <- matrix(CSAfragmentationList[order(CSAfragmentationList[, 2], decreasing = TRUE), ], ncol = 4)
                        ##
                        spectralEntropy <- round(spectral_entropy_calculator(CSAfragmentationList[, 1:2], allowedWeightedSpectralEntropy = TRUE, noiseRemovalRatio = 1e-16)[[1]], 5)
                        ##
                        numCSApeaks <- numCSApeaks + 1
                        counterCSAblock <- counterCSAblock + 1
                        ##
                        CSA_list <- cbind(rep(counterCSAblock, numCSApeaks), rep(mz12CIPA[i], numCSApeaks), rep(RTIPA[i], numCSApeaks),
                                          rep(IntIPA[i], numCSApeaks), CSAfragmentationList[, 1:2], rep(spectralEntropy, numCSApeaks),
                                          rep(MS1polarity[scanNumberApex], numCSApeaks), rep(precursorCE[scanNumberApex], numCSApeaks), CSAfragmentationList[, 3:4])
                        ##
                        listCSAalignedTable[[counterCSAblock]] <- CSA_list
                        ########################################################
                        kIPApeakID <- do.call(c, lapply(indexCSAfragmentIons, function(k) {CSA_EICs[[k]][[3]]}))
                        kIPApeakID <- kIPApeakID[kIPApeakID != 0]
                        kIPApeakID <- c(i, kIPApeakID)
                        detectedIPApeaks <- IPApeakID[kIPApeakID]
                        IPApeakID[kIPApeakID] <- 0
                        ########################################################
                        allowedSelectedIPApeaks <- TRUE
                        if (!is.null(selectedIPApeaks)) {
                          xSelectedIPApeaks <- which(selectedIPApeaks %in% detectedIPApeaks)
                          if (length(xSelectedIPApeaks) > 0) {
                            selectedIPApeaks[xSelectedIPApeaks] <- 0
                          } else {
                            allowedSelectedIPApeaks <- FALSE
                          }
                        }
                        ##
                        if (allowedSelectedIPApeaks) {
                          ######################################################
                          if (plotEICcheck) {
                            CSAEICdata <- do.call(rbind, lapply(indexCSAfragmentIons, function(j) {
                              CSA_EICs[[j]][[1]]
                            }))
                            ##
                            yMaxLimPlot <- max(c(Int_chrom_precursor, CSAEICdata[, 3]))
                            ##
                            xLinesDiff <- c(0, which(abs(diff(CSAEICdata[, 1])) > 0), nrow(CSAEICdata))
                            nLines <- length(xLinesDiff)
                            nLines1 <- nLines - 1
                            ##
                            legText <- rep("", nLines)
                            legText[1] <- paste0("* m/z = ", round(mz12CIPA[i], 5))
                            colors <- c("black", rainbow(nLines1, alpha = 1))
                            legCol <- rep("", nLines)
                            legCol[1] <- colors[1]
                            ##
                            orderMSP <- order(CSA_fragments[, 2], decreasing = TRUE)
                            orderLines <- do.call(rbind, lapply(2:nLines, function(p) {
                              c((xLinesDiff[p - 1] + 1), xLinesDiff[p])
                            }))
                            orderLines <- matrix(orderLines[orderMSP, ], ncol = 2)
                            ##
                            alignedEICfilename <- paste0(outputCSAeic, "/CSApeakGrouping_ID_", counterCSAblock, "_RT_", RTIPA[i], ".png")
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
                              lines(CSAEICdata[xLines, 2], CSAEICdata[xLines, 3], lwd = 2, col = colors[pCounter], cex = 4)
                              ##
                              legText[pCounter] <- paste0("m/z = ", round(mz_fragment[CSAEICdata[xLines[1], 1]], 5))
                              legCol[pCounter] <- colors[pCounter]
                            }
                            ##
                            mtext(text = paste0("S = ", spectralEntropy), side = 3, adj = 1, line = 0.25, cex = 1.0)
                            mtext("Retention time (min)", side = 1, adj = 0.5, line = 2, cex = 1.35)
                            mtext("Intensity", side = 2, adj = 0.5, line = 2, cex = 1.35)
                            mtext(CSA_hrms_file, side = 3, adj = 0, line = 0.25, cex = 1.4)
                            legend(x = "topright", inset = c(-0.20, 0), legend = legText, lwd = c(4, rep(2, nLines1)), cex = 1.125, bty = "n",
                                   col = legCol, seg.len = 1, x.intersp = 0.5, y.intersp = 0.9, xpd = TRUE)
                            ##
                            dev.off()
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ##
    if (!is.null(selectedIPApeaks)) {
      if (length(which(selectedIPApeaks != 0)) == 0) {
        break
      }
    }
    ##
    x_i <- which(IPApeakID[(i + 1):nPeaks] != 0)
    if (length(x_i) > 0) {
      i <- i + x_i[1]
    } else {
      i <- nPeaks
    }
  }
  ##
  if (counterCSAblock > 0) {
    CSA_peaklist <- do.call(rbind, lapply(listCSAalignedTable, function(i) {i}))
    ##
    CSA_peaklist <- data.frame(CSA_peaklist)
    CSA_peaklist[, 1] <- as.numeric(CSA_peaklist[, 1])
    CSA_peaklist[, 5] <- round(as.numeric(CSA_peaklist[, 5]), 5)
    CSA_peaklist[, 6] <- round(as.numeric(CSA_peaklist[, 6]), 0)
    CSA_peaklist[, 10] <- round(as.numeric(CSA_peaklist[, 10]), 2)
    CSA_peaklist[, 11] <- as.numeric(CSA_peaklist[, 11])
    ##
    if (alignedTableCheck) {
      ##
      xXcolID <- do.call(c, lapply(CSA_peaklist[, 11], function(i) {
        if (i != 0) {
          which(peakXcolID == i)
        } else {
          0
        }
      }))
      ##
      CSA_peaklist <- cbind(CSA_peaklist, xXcolID)
    }
    ##
  } else {
    if (alignedTableCheck) {
      CSA_peaklist <- data.frame(matrix(rep(0, 12), nrow = 1))
    } else {
      CSA_peaklist <- data.frame(matrix(rep(0, 11), nrow = 1))
    }
  }
  ##
  rownames(CSA_peaklist) <- NULL
  ##
  if (alignedTableCheck) {
    colnames(CSA_peaklist) <- c("CSApeakGrouping_ID", "PrecursorMZ", "Precursor_RT", "Precursor_Intensity",
                                "CSA_mz_fragment", "CSA_int_fragment", "Weighted_spectral_entropy_0noiseRemoval",
                                "Ion_mode", "Collision_energy", "Pearson_rho", "IDSL.IPA_PeakID", "IDSL.IPA_AlignedTable_PeakIDs")
  } else {
    colnames(CSA_peaklist) <- c("CSApeakGrouping_ID", "PrecursorMZ", "Precursor_RT", "Precursor_Intensity",
                                "CSA_mz_fragment", "CSA_int_fragment", "Weighted_spectral_entropy_0noiseRemoval",
                                "Ion_mode", "Collision_energy", "Pearson_rho", "IDSL.IPA_PeakID")
  }
  ##
  return(CSA_peaklist)
}