IDSL.CSA_MSPgenerator <- function(CSA_peaklist, msLevel, spectral_search_mode = "dda", spectral_search_mode_option = NA, number_processing_threads = 1) {
  ##
  DDAcheck <- FALSE
  CSAcheck <- FALSE
  if ((spectral_search_mode == "dda")) {
    MSP_MODE <- "DDA"
    DDAcheck <- TRUE
  } else if (spectral_search_mode == "dia") {
    MSP_MODE <- "DIA"
  } else if (spectral_search_mode == "csa") {
    MSP_MODE <- "CSA"
    CSAcheck <- TRUE
  }
  ##
  rawddaspectraCheck <- FALSE
  if (!is.na(spectral_search_mode_option)) {
    if (spectral_search_mode_option == "rawddaspectra") {
      rawddaspectraCheck <- TRUE
    }
  }
  ##
  if (CSAcheck) {
    alignedTableCheck <- FALSE
    if (!is.na(spectral_search_mode_option)) {
      if (spectral_search_mode_option == "alignedtable") {
        alignedTableCheck <- TRUE
      }
    }
  }
  ###########################################
  x_1 <- which(CSA_peaklist$Ion_mode == "1")
  if (length(x_1) > 0) {
    CSA_peaklist$Ion_mode[x_1] <- "Positive"
  }
  ##
  x_0 <- which(CSA_peaklist$Ion_mode == "0")
  if (length(x_0) > 0) {
    CSA_peaklist$Ion_mode[x_0] <- "Negative"
  }
  ###########################################
  ##
  ID <- as.numeric(CSA_peaklist[, 1])
  x_diffID <- c(0, which(abs(diff(ID)) > 0), length(ID))
  ##
  call_MSPvector <- function(i) {
    x_ID <- seq((x_diffID[i] + 1), x_diffID[i + 1], 1)
    ID_i <-  CSA_peaklist[x_ID[1], 1]
    ##
    if (rawddaspectraCheck) {
      MSPid <- paste0("Name: Precursor_scan_number_", ID_i, "_mz_", CSA_peaklist[x_ID[1], 2], "_RT_", CSA_peaklist[x_ID[1], 3], "\n")
    } else if (CSAcheck) {
      MSPid <- paste0("Name: CSApeakGrouping_ID_", ID_i, "_RT_", CSA_peaklist[x_ID[1], 3], "\n")
      MSPid <- paste0(MSPid, "CSApeakGrouping_ID: ", ID_i, "\n")
      MSPid <- paste0(MSPid, "Retention_time: ", CSA_peaklist[x_ID[1], 3], "\n")
    } else {
      MSPid <- paste0("Name: IDSL.IPA_PeakID_", ID_i, "_mz_", CSA_peaklist[x_ID[1], 2], "_RT_", CSA_peaklist[x_ID[1], 3], "\n")
    }
    ##
    MSPid <- paste0(MSPid, "MSP_mode: ", MSP_MODE, "\n")
    MSPid <- paste0(MSPid, "MS_level: ", msLevel, "\n")
    ##
    if (CSAcheck) {
      MSPid <- paste0(MSPid, "IDSL.IPA_Collective_PeakIDs: ", paste0(CSA_peaklist[x_ID, 11], collapse = ","), "\n")
      ##
      if (alignedTableCheck) {
        MSPid <- paste0(MSPid, "IDSL.IPA_AlignedTable_PeakIDs: ", paste0(CSA_peaklist[x_ID, 12], collapse = ","), "\n")
      }
      ##
    } else {
      MSPid <- paste0(MSPid, "IDSL.IPA_PeakID: ", ID_i, "\n")
      MSPid <- paste0(MSPid, "PrecursorMZ: ", CSA_peaklist[x_ID[1], 2], "\n")
      MSPid <- paste0(MSPid, "Precursor_RT: ", CSA_peaklist[x_ID[1], 3], "\n")
      MSPid <- paste0(MSPid, "Precursor_Intensity: ", CSA_peaklist[x_ID[1], 4], "\n")
    }
    ##
    MSPid <- paste0(MSPid, "basePeakMZ: ", CSA_peaklist[x_ID[1], 5], "\n")
    MSPid <- paste0(MSPid, "basePeakIntensity: ", CSA_peaklist[x_ID[1], 6], "\n")
    ##
    if (DDAcheck) {
      MSPid <- paste0(MSPid, "Count_DDA_scan: ", CSA_peaklist[x_ID[1], 10], "\n")
    } else {
      MSPid <- paste0(MSPid, "Pearson_rho: ", paste0(CSA_peaklist[x_ID, 10], collapse = ","), "\n")
    }
    ##
    MSPid <- paste0(MSPid, "Collision_energy: ", CSA_peaklist[x_ID[1], 9], "\n")
    ##
    MSPid <- paste0(MSPid, "Ion_mode: ", CSA_peaklist[x_ID[1], 8], "\n")    
    ##
    MSPid <- paste0(MSPid, "Weighted_spectral_entropy_0noiseRemoval: ", CSA_peaklist[x_ID[1], 7], "\n")
    ##
    MSPid <- paste0(MSPid, "Num Peaks: ", length(x_ID), "\n")
    ##
    MSPid_mz_int <- paste0(CSA_peaklist[x_ID, 5], " ", CSA_peaklist[x_ID, 6], "\n", collapse = "")
    ##
    paste0(MSPid, MSPid_mz_int, "\n")
  }
  ##############################################################################
  if (number_processing_threads == 1) {
    MSPvector <- do.call(c, lapply(1:(length(x_diffID) - 1), function(i) {
      call_MSPvector(i)
    }))
  } else {
    ##
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      MSPvector <- do.call(c, mclapply(1:(length(x_diffID) - 1), function(i) {
        call_MSPvector(i)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      MSPvector <- foreach(i = 1:(length(x_diffID) - 1), .combine = 'c', .verbose = FALSE) %dopar% {
        call_MSPvector(i)
      }
      ##
      stopCluster(clust)
    }
  }
  ##
  return(MSPvector)
}