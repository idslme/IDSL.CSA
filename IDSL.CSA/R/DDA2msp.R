DDA2msp <- function(input_path_hrms, file_name_hrms = NULL, number_processing_threads = 1) {
  ##
  if (is.null(file_name_hrms)) {
    file_name_hrms <- dir(path = input_path_hrms, pattern = ".mzML$|.mzXML$|.CDF$", ignore.case = TRUE)
  }
  ##
  call_DDA2msp <- function(input_path_hrms, iHRMSfilename) {
    DDA_peaklist <- DDA_rawSpectraDeconvolution(input_path_hrms, iHRMSfilename, rawDDAspectraVar = NULL, number_processing_threads = NPT)
    ##
    DDA_MSP <- IDSL.CSA_MSPgenerator(DDA_peaklist, msLevel = 2, spectral_search_mode = "dda", spectral_search_mode_option = "rawddaspectra", number_processing_threads = NPT)
    write.table(DDA_MSP, file = paste0(input_path_hrms, "/DDA_MSP_", iHRMSfilename, ".msp"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
    ##
    return()
  }
  ##
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    NPT <- 1
    ##
    iCounter <- 0
    progressBARboundaries <- txtProgressBar(min = 0, max = length(file_name_hrms), initial = 0, style = 3)
    for (i in file_name_hrms) {
      ##
      null_variable <- tryCatch(call_DDA2msp(input_path_hrms, iHRMSfilename = i),
                                error = function(e) {FSA_message(paste0("Problem with `", i,"`!"))})
      ##
      iCounter <- iCounter + 1
      setTxtProgressBar(progressBARboundaries, iCounter)
    }
    close(progressBARboundaries)
    ##
  } else {
    NPT0 <- number_processing_threads
    NPT <- 1
    ##
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      null_variable <- mclapply(file_name_hrms, function(i) {
        ##
        tryCatch(call_DDA2msp(input_path_hrms, iHRMSfilename = i),
                 error = function(e) {FSA_message(paste0("Problem with `", i,"`!"))})
      }, mc.cores = NPT0)
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      ##
      clust <- makeCluster(NPT0)
      registerDoParallel(clust)
      ##
      null_variable <- foreach(i = file_name_hrms, .verbose = FALSE) %dopar% {
        ##
        tryCatch(call_DDA2msp(input_path_hrms, iHRMSfilename = i),
                 error = function(e) {FSA_message(paste0("Problem with `", i,"`!"))})
      }
      ##
      stopCluster(clust)
      ##
    }
  }
  ##
  ##############################################################################
  ##
}