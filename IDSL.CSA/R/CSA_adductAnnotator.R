CSA_adductAnnotator <- function(IPApeakList, CSA_peaklist, massError) {
  ##
  system_path <- system.file("data", package = "IDSL.CSA")
  ploarity <- CSA_peaklist[1, 8]
  if ((ploarity == "1") | (ploarity == 1)) {
    REFadduct <- IDSL.IPA::loadRdata(paste0(system_path, "/positiveAdducts.rda"))
  } else {
    REFadduct <- IDSL.IPA::loadRdata(paste0(system_path, "/negativeAdducts.rda"))
  }
  ##
  nIPApeakList <- dim(IPApeakList)[1]
  nCSAid <- max(CSA_peaklist[, 1])
  xDiff <- which(abs(diff(CSA_peaklist[, 1])) > 0)
  #
  xCSAid <- matrix(rep(0, 2*nCSAid), ncol = 2)
  #
  uCSAid <- seq(1, nCSAid, 1)
  xCSAid[uCSAid, 1] <- c(1, (xDiff + 1))
  xCSAid[uCSAid, 2] <- c(xDiff, dim(CSA_peaklist)[1])
  ##
  CSA_spectra_peaklist <- matrix(rep(0, 2*nIPApeakList), ncol = 2)
  strCSA_spectra_peaklist <- rep("", nIPApeakList)
  for (j in 1:nCSAid) {
    IPApeakIDs <- CSA_peaklist[xCSAid[j, 1]:xCSAid[j, 2], 11]
    IPApeakIDs <- IPApeakIDs[IPApeakIDs != 0]
    nIPApeakIDs <- length(IPApeakIDs)
    jIPApeakIDs <- rep(j, nIPApeakIDs)
    CSA_spectra_peaklist[IPApeakIDs, ] <- cbind(IPApeakIDs, jIPApeakIDs)
    strCSA_spectra_peaklist[IPApeakIDs] <- paste0("CSA_ID_", jIPApeakIDs, "_RT_", rep(CSA_peaklist[xCSAid[j, 1], 3], nIPApeakIDs))
  }
  ##
  CSA_spectra_IDs <- CSA_spectra_peaklist[, 2]
  CSA_spectra_peaklist <- CSA_spectra_peaklist[CSA_spectra_peaklist[, 1] != 0, ]
  CSA_spectra_peaklist <- CSA_spectra_peaklist[order(CSA_spectra_peaklist[, 2], decreasing = FALSE), ]
  ##############################################################################
  xDiff <- which(abs(diff(CSA_spectra_peaklist[, 2])) > 0)
  #
  xCSAid <- matrix(rep(0, 2*nCSAid), ncol = 2)
  #
  uCSAid <- unique(CSA_spectra_peaklist[, 2])
  xCSAid[uCSAid, 1] <- c(1, (xDiff + 1))
  xCSAid[uCSAid, 2] <- c(xDiff, dim(CSA_spectra_peaklist)[1])
  ##
  matched_adducts <- do.call(rbind, lapply(1:(length(xDiff) + 1), function(j) {
    IPApeakIDs <- CSA_spectra_peaklist[xCSAid[j, 1]:xCSAid[j, 2], 1]
    IPApeakIDs <- sort(IPApeakIDs, decreasing = FALSE)
    ##
    adduct_ids <- do.call(rbind, lapply(IPApeakIDs, function(parent) {
      do.call(rbind, lapply(setdiff(IPApeakIDs, parent), function(ion) {
        mass12Cdiff <- IPApeakList[parent, 8] - IPApeakList[ion, 8]
        ##
        x <- which(abs(mass12Cdiff - REFadduct[, 2]) <= massError)
        ##
        if (length(x) > 0) {
          c(parent, ion, paste0(REFadduct[x, 1], collapse = ", "))
        }
      }))
    }))
    ##
    if (!is.null(adduct_ids)) {
      adduct_ids <- data.frame(adduct_ids)
      adduct_ids[, 1] <- as.numeric(adduct_ids[, 1])
      ##
      xDiff <- which(abs(diff(adduct_ids[, 1])) > 0)
      nIPAid <- length(unique(adduct_ids[, 1]))
      uIPAid <- seq(1, nIPAid, 1)
      #
      xIPAid <- matrix(rep(0, 2*nIPAid), ncol = 2)
      #
      xIPAid[uIPAid, 1] <- c(1, (xDiff + 1))
      xIPAid[uIPAid, 2] <- c(xDiff, dim(adduct_ids)[1])
      ##
      do.call(rbind, lapply(uIPAid, function(id) {
        c(adduct_ids[xIPAid[id, 1], 1], paste0(adduct_ids[xIPAid[id, 1]:xIPAid[id, 2], 2], collapse = ", "), paste0(adduct_ids[xIPAid[id, 1]:xIPAid[id, 2], 3], collapse = ", "))
      }))
    }
  }))
  ##############################################################################
  if (!is.null(matched_adducts)) {
    matched_adducts <- data.frame(matched_adducts)
    matched_adducts[, 1] <- as.numeric(matched_adducts[, 1])
    ##
    annex_CSA <- matrix(rep("", 2*nIPApeakList), ncol = 2)
    annex_CSA[matched_adducts[, 1], 1] <- matched_adducts[, 2]
    annex_CSA[matched_adducts[, 1], 2] <- matched_adducts[, 3]
    ##
    annex_CSA <- cbind(strCSA_spectra_peaklist, annex_CSA)
    annex_CSA <- data.frame(annex_CSA, stringsAsFactors = FALSE)
    colnames(annex_CSA) <- c("CSA_spectra_IDs", "corrolated_IPA_peak_ID", "adducts")
    ##
    IPApeakList <- cbind(IPApeakList, annex_CSA)
    ##
    for (j in 1:24) {
      IPApeakList[, j] <- as.numeric(IPApeakList[, j])
    }
  }
  ##
  return(IPApeakList)
}