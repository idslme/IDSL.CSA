# IDSL.CSA <img src='CSA_educational_files/Figures/IDSL.CSA-logo.PNG' width="250px" align="right" />

<!-- badges: start -->
[![Maintainer](https://img.shields.io/badge/maintainer-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
<!-- badges: end -->

[**Composite Spectra Analysis (CSA)**](https://www.csa.idsl.me/) by the [**Integrated Data Science Laboratory for Metabolomics and Exposomics (IDSL.ME)**](https://www.idsl.me/) is an R package designed to deconvolute fragmentation spectra from Composite Spectra Analysis (CSA), Data Dependent Acquisition (DDA), and various Data-Independent Acquisition (DIA) methods such as MS<sup>E</sup>, and All-Ion Fragmentation (AIF).

	install.packages("IDSL.CSA") # IDSL.CSA package is set to release on CRAN by the end of October

## Workflow
Prior to processing your mass spectrometry data (**mzXML**, **mzML**, **netCDF**) using the IDSL.CSA workflow, mass spectrometry data should be processed using the [IDSL.IPA](https://github.com/idslme/IDSL.IPA) workflow to acquire chromatographic information of the peaks (***m/z-RT***). When the chromatographic information of individual and aggregated aligned peaklists were generated using the [IDSL.IPA](https://github.com/idslme/IDSL.IPA) workflow, download the [CSA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.CSA/main/CSA_parameters.xlsx) and select the parameters accordingly and then use this spreadsheet as the input for the IDSL.CSA workflow:

	library(IDSL.CSA)
	IDSL.CSA_workflow("Address of the CSA parameter spreadsheet")


Visit https://csa.idsl.me/ for the detailed documentation and tutorial.