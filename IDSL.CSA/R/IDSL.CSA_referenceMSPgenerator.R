IDSL.CSA_referenceMSPgenerator <- function(REF_peaklist, refTable, selectedIPApeaks_IDref, msLevel, spectral_search_mode = "dda", spectral_search_mode_option = NA) {
  ##
  DIAcheck <- FALSE
  CSAcheck <- FALSE
  if ((spectral_search_mode == "dda")) {
    MSP_MODE <- "DDA"
  } else if (spectral_search_mode == "dia") {
    MSP_MODE <- "DIA"
    DIAcheck <- TRUE
  } else if (spectral_search_mode == "csa") {
    MSP_MODE <- "CSA"
    CSAcheck <- TRUE
    ##
    alignedTableCheck <- FALSE
    if (!is.na(spectral_search_mode_option)) {
      if (spectral_search_mode_option == "alignedtable") {
        alignedTableCheck <- TRUE
      }
    }
  }
  ##
  rawDDAspectraCheck <- FALSE
  if (!is.na(spectral_search_mode_option)) {
    if (spectral_search_mode_option == "rawddaspectra") {
      rawDDAspectraCheck <- TRUE
    }
  }
  ##
  x_1 <- which(REF_peaklist$Ion_mode == "1")
  if (length(x_1) > 0) {
    REF_peaklist$Ion_mode[x_1] <- "Positive"
  }
  ##
  x_0 <- which(REF_peaklist$Ion_mode == "0")
  if (length(x_0) > 0) {
    REF_peaklist$Ion_mode[x_0] <- "Negative"
  }
  ##
  ##############################################################################
  ##
  if (rawDDAspectraCheck) {
    REF_peaklist <- do.call(rbind, lapply(1:nrow(selectedIPApeaks_IDref), function(j) {
      x_j <- which(REF_peaklist$PrecursorScanNumber == selectedIPApeaks_IDref[j, 1])
      l_x_j <- length(x_j)
      if (l_x_j > 0) {
        do.call(rbind, lapply(x_j, function(rowx_j) {
          cbind(j, REF_peaklist[rowx_j, ], refTable[selectedIPApeaks_IDref[j, 2], ])
        }))
      }
    }))
  } else {
    REF_peaklist <- do.call(rbind, lapply(1:nrow(selectedIPApeaks_IDref), function(j) {
      x_j <- which(REF_peaklist$IDSL.IPA_PeakID == selectedIPApeaks_IDref[j, 1])
      l_x_j <- length(x_j)
      if (l_x_j > 0) {
        do.call(rbind, lapply(x_j, function(rowx_j) {
          cbind(j, REF_peaklist[rowx_j, ], refTable[selectedIPApeaks_IDref[j, 2], ])
        }))
      }
    }))
  }
  ##
  colnames(REF_peaklist)[1] <- 'CSA_reference_ID'
  ##
  ##############################################################################
  ##
  CSAcolNames <- colnames(REF_peaklist)
  x_IDcol <- which(CSAcolNames == 'CSA_reference_ID')
  x_mzCol <- which(CSAcolNames == 'CSA_mz_fragment')
  x_intCol <- which(CSAcolNames == 'CSA_int_fragment')
  x_nameCol <- which(CSAcolNames == 'Name')
  ##
  if ((length(x_IDcol) != 1) | (length(x_mzCol) != 1) | (length(x_intCol) != 1) | (length(x_nameCol) != 1)) {
    stop("The REF_peaklist file should have only one column for the following headers (case-senstive):
    'CSA_reference_ID'
    'CSA_mz_fragment'
    'CSA_int_fragment'
    'Name'")
  }
  ##
  CSAcolumns <- setdiff(CSAcolNames, c('CSA_reference_ID', 'CSA_mz_fragment', 'CSA_int_fragment', 'Name',
                                       'IDSL.IPA_Collective_PeakIDs', 'IDSL.IPA_AlignedTable_PeakIDs', 'Pearson_rho'))
  if (length(CSAcolumns) > 0) {
    CSAcolumnsCheck <- TRUE
    x_CSAcolumns <- do.call(c, lapply(CSAcolumns, function(i) {
      which(CSAcolNames == i)
    }))
  } else {
    CSAcolumnsCheck <- FALSE
  }
  ##
  CSA_reference_ID <- as.numeric(REF_peaklist[, x_IDcol])
  x_diffID <- c(0, which(abs(diff(CSA_reference_ID)) > 0), length(CSA_reference_ID))
  ##
  MSPvector <- do.call(c, lapply(1:(length(x_diffID) - 1), function(i) {
    x_ID <- seq((x_diffID[i] + 1), x_diffID[i + 1], 1)
    ID_i <-  x_diffID[i] + 1
    ##
    MSP1 <- paste0("Name: ", REF_peaklist[ID_i, x_nameCol], "\n")
    MSP1 <- paste0(MSP1, "MSP_mode: ", MSP_MODE, "\n")
    MSP1 <- paste0(MSP1, "MS_level: ", msLevel, "\n")
    MSP1 <- paste0(MSP1, "basePeakMZ: ", REF_peaklist[x_ID[1], x_mzCol], "\n")
    MSP1 <- paste0(MSP1, "basePeakIntensity: ", REF_peaklist[x_ID[1], x_intCol], "\n")
    ##
    if (CSAcolumnsCheck) {
      MSPid <- do.call(paste0, lapply(x_CSAcolumns, function(j) {
        paste0(CSAcolNames[j], ": ", REF_peaklist[ID_i, j], "\n")
      }))
    } else {
      MSPid <- NULL
    }
    ##
    if (CSAcheck | DIAcheck) {
      MSPid <- paste0(MSPid, "Pearson_rho: ", paste0(REF_peaklist$`Pearson_rho`[x_ID], collapse = ","), "\n")
    }
    ##
    if (CSAcheck) {
      MSPid <- paste0(MSPid, "IDSL.IPA_Collective_PeakIDs: ", paste0(REF_peaklist$IDSL.IPA_Collective_PeakIDs[x_ID], collapse = ","), "\n")
      ##
      if (alignedTableCheck) {
        MSPid <- paste0(MSPid, "IDSL.IPA_AlignedTable_PeakIDs: ", paste0(REF_peaklist$IDSL.IPA_AlignedTable_PeakIDs[x_ID], collapse = ","), "\n")
      }
      ##
    }
    ##
    MSPid <- paste0(MSPid, "Num Peaks: ", length(x_ID), "\n")
    ##
    MSPid_mz_int <- paste0(REF_peaklist[x_ID, x_mzCol], " ", REF_peaklist[x_ID, x_intCol], "\n", collapse = "")
    ##
    paste0(MSP1, MSPid, MSPid_mz_int, "\n")
  }))
  ##
  return(MSPvector)
}