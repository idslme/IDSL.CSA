CSA_alignedPeaksTanimotoCoefficientCalculator <- function(address_input_msp, peakXcol, minPercenetageDetection = 5, minNumberFragments = 2,
                                                          minTanimotoCoefficient = 0.1, RTtolerance = 0.05, number_processing_threads = 1) {
  ##
  FSA_logRecorder("Initiated clustering collected IDSL.IPA peaks from the CSA workflow on the aligned peak table!")
  ##
  file_name_sample_msp <- dir(path = address_input_msp)
  file_name_sample_msp <- file_name_sample_msp[grepl(".msp$", file_name_sample_msp, ignore.case = TRUE)]
  L_MSP <- length(file_name_sample_msp)
  ##
  L_peakXcol <- dim(peakXcol)[1]
  rep0LpeakXcol <- rep(0, L_peakXcol)
  Lsamples <- dim(peakXcol)[2]
  ##############################################################################
  ##############################################################################
  ColPL <- colnames(peakXcol)[4:Lsamples]
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
    Lsamples <- Lsamples - nMissedPL
  }
  ##
  seqXcolMSP <- matrix(seqXcolMSP[order(seqXcolMSP[, 2], decreasing = FALSE), ], ncol = 2)
  orderSeqMSP <- seqXcolMSP[, 1]
  seqMSP <- seqXcolMSP[, 2]
  peakXcol <- peakXcol[, c(seq(1, 3, 1), (orderSeqMSP + 3))]
  ##############################################################################
  ##############################################################################
  CSAmsp2FSdb <- function(path, mspFileName) {
    ##
    msp <- readLines(paste0(path, "/", mspFileName), warn = FALSE)
    ##
    loc_collectivePeakIDs <- grep("IDSL.IPA_Collective_PeakIDs: ", msp, ignore.case = TRUE)
    ##
    if (length(loc_collectivePeakIDs) == 0) {
      FSA_logRecorder(paste0("WARNING!!! 'IDSL.IPA_Collective_PeakIDs' rows are not available in ", mspFileName))
    }
    ##
    strCollectivePeakIDs <- gsub("IDSL.IPA_Collective_PeakIDs: " , "", msp[loc_collectivePeakIDs])
    ##
    lapply(strCollectivePeakIDs, function(i) {
      IDj <- eval(parse(text = paste0("c(", i, ")")))
      IDj[IDj != 0]
    })
  }
  ##############################################################################
  ##############################################################################
  call_peakXmsp <- function(i) {
    ##
    peak_table_id <- peakXcol[, (i + 3)]
    ##
    peakIDj <- tryCatch(CSAmsp2FSdb(address_input_msp, file_name_sample_msp[i]),
                        warning = function(w) {stop(message(paste0("problem with ", file_name_sample_msp[i]),
                                                            "/nProbabely 'IDSL.IPA_Collective_PeakIDs' rows are not available in the .msp files!"))})
    ##
    codetectedIDs <- rep(0, L_peakXcol)
    counter <- 0
    for (IDj in peakIDj) {
      ##
      x_IDj <- which(peak_table_id %in% IDj)
      if (length(x_IDj) > 0) {
        counter <- counter + 1
        codetectedIDs[x_IDj] <- counter
      }
    }
    return(codetectedIDs)
  }
  ##############################################################################
  ##############################################################################
  call_IDTC <- function(i) {
    x_rt <- which(abs(RTX - RTX[i]) <= RTtolerance)
    ##
    if (length(x_rt) >= minNumberFragments) {
      x_rt <- setdiff(x_rt, i)
      ##
      do.call(rbind, lapply(x_rt, function(j) {
        A <- length(which(peakXlist[[i]] %in% peakXlist[[j]]))
        if (A > 0) {
          B <- SX[i] - A
          C <- SX[j] - A
          TanimotoCoefficient <- A/(A + B + C)
          if (TanimotoCoefficient >= minTanimotoCoefficient) {
            c(IDX[i], IDX[j], TanimotoCoefficient)
          }
        }
      }))
    }
  }
  ##############################################################################
  ##############################################################################
  if (number_processing_threads == 1) {
    ##
    ############################################################################
    ##
    peakXmsp <- do.call(cbind, lapply(1:L_MSP, function(i) {
      iCheck <- i %in% seqMSP
      if (iCheck) {
        call_peakXmsp(i)
      } else {
        rep0LpeakXcol
      }
    }))
    ##
    ############################################################################
    ##
    l_non0 <- max(peakXmsp[, 1])
    for (i in 2:L_MSP) {
      x_non0 <- which(peakXmsp[, i] > 0)
      peakXmsp[x_non0, i] <- peakXmsp[x_non0, i] + l_non0 + 1
      l_non0 <- max(c(peakXmsp[, i], l_non0))
    }
    ##
    ############################################################################
    if (l_non0 > 0) {
      ##
      ##########################################################################
      ##
      numberCSAdetFreq <- do.call(c, lapply(1:L_peakXcol, function(i) {
        length(which(peakXmsp[i, ] > 0))
      }))
      ##
      xDetFreq <- which(numberCSAdetFreq/L_MSP >= minPercenetageDetection/100)
      ##
      peakXmsp <- cbind(xDetFreq, peakXcol[xDetFreq, 2:3], peakXmsp[xDetFreq, ])
      peakXmsp <- peakXmsp[order(peakXmsp[, 1], decreasing = TRUE), ]
      peakXmsp <- peakXmsp[order(peakXmsp[, 2], decreasing = TRUE), ]
      ##
      peakXcol <- NULL
      LXDF <- dim(peakXmsp)[1]
      ##
      ##########################################################################
      ##
      IDX <- peakXmsp[, 1]
      RTX <- peakXmsp[, 2]
      SX <- peakXmsp[, 3]
      ##
      ##########################################################################
      ##
      peakXlist <- lapply(1:LXDF, function(i) {
        px <- peakXmsp[i, 4:Lsamples]
        px[px != 0]
      })
      ##
      names(peakXlist) <- as.character(IDX)
      ##
      peakXmsp <- NULL
      ##
      ##########################################################################
      ##
      FSA_logRecorder("Initiated calculating Tanimoto coefficients!")
      ##
      progressBARboundaries <- txtProgressBar(min = 0, max = LXDF, initial = 0, style = 3)
      ##
      IDTC <- do.call(rbind, lapply(1:LXDF, function(i) {
        ##
        setTxtProgressBar(progressBARboundaries, i)
        ##
        call_IDTC(i)
      }))
      ##
      close(progressBARboundaries)
      ##
    } else {
      ##
      IDTC <- NULL
    }
    ##
    ############################################################################
    ##
  } else {
    ##
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Windows") {
      ##
      ##########################################################################
      ##
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, setdiff(ls(), c("clust", "L_MSP")), envir = environment())
      ##
      peakXmsp <- do.call(cbind, parLapply(clust, 1:L_MSP, function(i) {
        iCheck <- i %in% seqMSP
        if (iCheck) {
          call_peakXmsp(i)
        } else {
          rep0LpeakXcol
        }
      }))
      ##
      stopCluster(clust)
      ##
      ##########################################################################
      ##
      l_non0 <- max(peakXmsp[, 1])
      for (i in 2:L_MSP) {
        x_non0 <- which(peakXmsp[, i] > 0)
        peakXmsp[x_non0, i] <- peakXmsp[x_non0, i] + l_non0 + 1
        l_non0 <- max(c(peakXmsp[, i], l_non0))
      }
      ##
      ##########################################################################
      if (l_non0 > 0) {
        ##
        ########################################################################
        ##
        clust <- makeCluster(number_processing_threads)
        clusterExport(clust, c("peakXmsp"), envir = environment())
        ##
        numberCSAdetFreq <- do.call(c, parLapply(clust, 1:L_peakXcol, function(i) {
          length(which(peakXmsp[i, ] > 0))
        }))
        ##
        stopCluster(clust)
        ##
        xDetFreq <- which(numberCSAdetFreq/L_MSP >= minPercenetageDetection/100)
        ##
        peakXmsp <- cbind(xDetFreq, peakXcol[xDetFreq, 2:3], peakXmsp[xDetFreq, ])
        peakXmsp <- peakXmsp[order(peakXmsp[, 1], decreasing = TRUE), ]
        peakXmsp <- peakXmsp[order(peakXmsp[, 2], decreasing = TRUE), ]
        ##      
        peakXcol <- NULL
        LXDF <- dim(peakXmsp)[1]
        ##
        ########################################################################
        ##
        IDX <- peakXmsp[, 1]
        RTX <- peakXmsp[, 2]
        SX <- peakXmsp[, 3]
        ##
        ########################################################################
        ##
        clust <- makeCluster(number_processing_threads)
        clusterExport(clust, c("peakXmsp", "Lsamples"), envir = environment())
        ##
        peakXlist <- parLapply(clust, 1:LXDF, function(i) {
          px <- peakXmsp[i, 4:Lsamples]
          px[px != 0]
        })
        ##
        stopCluster(clust)
        ##
        names(peakXlist) <- as.character(IDX)
        ##
        peakXmsp <- NULL
        ##
        ########################################################################
        ##
        FSA_logRecorder("Initiated calculating Tanimoto coefficients!")
        ##
        clust <- makeCluster(number_processing_threads)
        clusterExport(clust, setdiff(ls(), c("clust", "LXDF")), envir = environment())
        ##
        IDTC <- do.call(rbind, parLapply(clust, 1:LXDF, function(i) {
          call_IDTC(i)
        }))
        ##
        stopCluster(clust)
        ##
      } else {
        ##
        IDTC <- NULL
      }
      ##
      ##########################################################################
      ##
    } else {
      ##
      ##########################################################################
      ##
      peakXmsp <- do.call(cbind, mclapply(seqMSP, function(i) {
        iCheck <- i %in% seqMSP
        if (iCheck) {
          call_peakXmsp(i)
        } else {
          rep0LpeakXcol
        }
      }, mc.cores = number_processing_threads))
      ##
      ##########################################################################
      ##
      l_non0 <- max(peakXmsp[, 1])
      for (i in 2:L_MSP) {
        x_non0 <- which(peakXmsp[, i] > 0)
        peakXmsp[x_non0, i] <- peakXmsp[x_non0, i] + l_non0 + 1
        l_non0 <- max(c(peakXmsp[, i], l_non0))
      }
      ##
      ##########################################################################
      if (l_non0 > 0) {
        ##
        ########################################################################
        ##
        numberCSAdetFreq <- do.call(c, mclapply(1:L_peakXcol, function(i) {
          length(which(peakXmsp[i, ] > 0))
        }, mc.cores = number_processing_threads))
        ##
        xDetFreq <- which(numberCSAdetFreq/L_MSP >= minPercenetageDetection/100)
        ##
        peakXmsp <- cbind(xDetFreq, peakXcol[xDetFreq, 2:3], peakXmsp[xDetFreq, ])
        peakXmsp <- peakXmsp[order(peakXmsp[, 1], decreasing = TRUE), ]
        peakXmsp <- peakXmsp[order(peakXmsp[, 2], decreasing = TRUE), ]
        ##
        peakXcol <- NULL
        LXDF <- dim(peakXmsp)[1]
        ##
        ########################################################################
        ##
        IDX <- peakXmsp[, 1]
        RTX <- peakXmsp[, 2]
        SX <- peakXmsp[, 3]
        ##
        ########################################################################
        ##
        peakXlist <- mclapply(1:LXDF, function(i) {
          px <- peakXmsp[i, 4:Lsamples]
          px[px != 0]
        }, mc.cores = number_processing_threads)
        ##
        names(peakXlist) <- as.character(IDX)
        ##
        peakXmsp <- NULL
        ##
        ########################################################################
        ##
        FSA_logRecorder("Initiated calculating Tanimoto coefficients!")
        ##
        IDTC <- do.call(rbind, mclapply(1:LXDF, function(i) {
          call_IDTC(i)
        }, mc.cores = number_processing_threads))
        ##
      } else {
        ##
        IDTC <- NULL
      }
      ##########################################################################
      ##
      closeAllConnections()
      ##
    }
  }
  ##############################################################################
  ##############################################################################
  ##############################################################################
  if (!is.null(IDTC)) {
    ##
    xRemoveMinFreq <- 0
    ##
    while (length(xRemoveMinFreq) > 0) {
      ##
      if (length(IDTC) > 0) {
        tIDTC <- table(c(IDTC[, 1], IDTC[, 2]))
        ##
        tMinFreqIDTC <- as.numeric(names(tIDTC[tIDTC < minNumberFragments]))
        xRemoveMinFreq <- which((IDTC[, 1] %in% tMinFreqIDTC) |
                                  (IDTC[, 2] %in% tMinFreqIDTC))
        if (length(xRemoveMinFreq) > 0) {
          IDTC <- matrix(IDTC[-xRemoveMinFreq, ], ncol = 3)
        }
      } else {
        IDTC <- NULL
        xRemoveMinFreq <- NULL
      }
    }
    ##
    FSA_logRecorder("Completed calculating Tanimoto coefficients!")
  }
  ##############################################################################
  groupedID_peakXcol <- matrix(rep(0, 3*L_peakXcol), ncol = 3)
  ##  
  if (!is.null(IDTC)) {
    ##
    ############################################################################
    ##
    IDTC <- IDTC[order(IDTC[, 3], decreasing = TRUE), ]
    tIDTC <- sort(table(c(IDTC[, 1], IDTC[, 2])), decreasing = TRUE)
    idtIDTC <- as.numeric(names(tIDTC))
    LidtIDTC <- length(idtIDTC)
    ##
    ############################################################################
    ##
    call_listXIDTC <- function(i) {
      x_i1 <- which(idtIDTC[i] == IDTC[, 1])
      x_i2 <- which(idtIDTC[i] == IDTC[, 2])
      list(x_i1, x_i2)
    }
    ##
    ############################################################################
    ##
    FSA_logRecorder("Initiated clustering aligned peaks!")
    ##
    if (number_processing_threads == 1) {
      ##
      listXIDTC <- lapply(1:LidtIDTC, function(i) {
        call_listXIDTC(i)
      })
      ##
    } else {
      ##
      if (osType == "Windows") {
        ##
        clust <- makeCluster(number_processing_threads)
        clusterExport(clust, c("call_listXIDTC", "idtIDTC", "IDTC"), envir = environment())
        ##
        listXIDTC <- parLapply(clust, 1:LidtIDTC, function(i) {
          call_listXIDTC(i)
        })
        ##
        stopCluster(clust)
        ##
      } else {
        listXIDTC <- mclapply(1:LidtIDTC, function(i) {
          call_listXIDTC(i)
        }, mc.cores = number_processing_threads)
        ##
        closeAllConnections()
        ##
      }
    }
    ##
    ############################################################################
    ##
    progressBARboundaries <- txtProgressBar(min = 0, max = LidtIDTC, initial = 0, style = 3)
    ##
    rootPowerTC <- minTanimotoCoefficient^(1/minNumberFragments)
    ##
    for (i in 1:LidtIDTC) {
      ##
      x_i1 <- listXIDTC[[i]][[1]]
      x_i2 <- listXIDTC[[i]][[2]]
      x_i <- c(x_i1, x_i2)
      ##
      jCheck <- TRUE
      ##
      for (j in x_i) {
        for (k in 1:2) {
          if (IDTC[j, 3] > (groupedID_peakXcol[IDTC[j, k], 2]/rootPowerTC)) {
            groupedID_peakXcol[IDTC[j, k], 1] <- i
            groupedID_peakXcol[IDTC[j, k], 2] <- IDTC[j, 3]
            if (jCheck) {
              iCluster <- peakXlist[[as.character(IDTC[j, k])]]
              groupedID_peakXcol[IDTC[j, k], 3] <- ceiling(length(iCluster)*IDTC[j, 3])
              ##
              jCheck <- FALSE
            } else {
              jMatch <- which((iCluster %in% peakXlist[[as.character(IDTC[j, k])]]))
              groupedID_peakXcol[IDTC[j, k], 3] <- ceiling(length(jMatch)*IDTC[j, 3])
            }
          }
        }
      }
      ##
      setTxtProgressBar(progressBARboundaries, i)
    }
    ##
    close(progressBARboundaries)
    ##
    FSA_logRecorder("Completed clustering aligned peaks!")
    ############################################################################
    tIDTC <- table(groupedID_peakXcol[, 1])
    tMinFreqIDTC <- tryCatch(as.numeric(names(tIDTC[tIDTC < minNumberFragments])), warning = function(w){0})
    x0 <- groupedID_peakXcol[, 1] %in% tMinFreqIDTC
    groupedID_peakXcol[x0, ] <- 0
  }
  ##############################################################################
  ##############################################################################
  groupedID_peakXcol[, 2] <- round(groupedID_peakXcol[, 2], digits = 3)
  groupedID_peakXcol[which(groupedID_peakXcol[, 1] == 0), 1] <- NA
  groupedID_peakXcol <- data.frame(groupedID_peakXcol, stringsAsFactors = FALSE)
  colnames(groupedID_peakXcol) <- c("coDetectedGroupingID", "TanimotoCoefficient", "CSAdetectionFrequency")
  ##
  FSA_logRecorder("Completed clustering collected IDSL.IPA peaks from the CSA workflow on the aligned peak table!")
  ##
  return(groupedID_peakXcol)
}