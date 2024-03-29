# IDSL.CSA <img src='CSA_educational_files/Figures/IDSL.CSA-logo.PNG' width="250px" align="right" />

<!-- badges: start -->
[![Developed-by](https://img.shields.io/badge/Developed_by-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
[![CRAN status](https://www.r-pkg.org/badges/version/IDSL.CSA)](https://cran.r-project.org/package=IDSL.CSA)
![](http://cranlogs.r-pkg.org/badges/IDSL.CSA?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/IDSL.CSA?color=brightgreen)
[![Dependencies](https://tinyverse.netlify.com/badge/IDSL.CSA)](https://cran.r-project.org/package=IDSL.CSA)
<!-- badges: end -->

The **Composite Spectra Analysis (IDSL.CSA)** R package for the analysis of mass spectrometry data has been developed by the [**Integrated Data Science Laboratory for Metabolomics and Exposomics (IDSL.ME)**](https://www.idsl.me/). This package can be used for the deconvolution of fragmentation spectra obtained through various analytical methods such as MS1-only Composite Spectra deconvolution Analysis (**CSA**), Data Dependent Acquisition (**DDA**), and a various Data-Independent Acquisition (**DIA**) methods including MS<sup>E</sup>, All-Ion Fragmentation (AIF), and SWATH-MS analyses. The aim of the **IDSL.CSA** package is to assist in streamlining the data analysis process and improving the overall chemical structure annotation in the fields of metabolomics and exposomics.

## <img src='CSA_educational_files/Figures/IDSL.CSA-TOC_Art.png' align="right" />

## Table of Contents

- [Features of IDSL.CSA](https://github.com/idslme/IDSL.CSA#features-of-idslcsa)
- [Installation](https://github.com/idslme/IDSL.CSA#installation)
- [Workflow](https://github.com/idslme/IDSL.CSA#workflow)
- [Quick Batch Example](https://github.com/idslme/IDSL.CSA#quick-batch-example)
- [Wiki](https://github.com/idslme/IDSL.CSA#wiki)
- [Citation](https://github.com/idslme/IDSL.CSA#citation)

## Features of IDSL.CSA

1) Parameter selection through a user-friendly and well-described [parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.CSA/main/CSA_parameters.xlsx)
2) Peak detection and chromatogram deconvolution for various fragmentation data analyses including [Composite Spectra Analysis (CSA)](https://github.com/idslme/IDSL.CSA/wiki/CSA-analysis-by-IDSL.CSA), [Data Dependent Acquisition (DDA)](https://github.com/idslme/IDSL.CSA/wiki/DDA-analysis-by-IDSL.CSA), and [Data-Independent Acquisition (DIA)](https://github.com/idslme/IDSL.CSA/wiki/DIA-analysis-by-IDSL.CSA)
3) Analyzing population size untargeted studies (n > 500)
4) Aggregating annotated chemical structures on the aligned peak table using meta-variables such as InChIKey, SMILES, precursor type, molecular formula,... depending on the information in the reference library. This is a very unique feature that is only presented by IDSL.CSA. To familiarize with this statistical mass spectrometry feature, try **PARAM0006** in the `Start` tab in the [IDSL.CSA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.CSA/main/CSA_parameters.xlsx).
5) Generating batch untargeted aligned extracted ion chromatograms (EIC) figures for the DIA and CSA analyses in addition to generating batch DDA spectra figures.
6) Parallel processing in Windows and Linux environments
7) Integration with [IDSL.FSA](https://github.com/idslme/IDSL.FSA) workflow to annotate various types of MSP files and generating fragmentation libraries.

## Installation

	install.packages("IDSL.CSA")

## Workflow

Prior to processing your mass spectrometry data (**mzXML**, **mzML**, **netCDF**) using the IDSL.CSA workflow, mass spectrometry data should be processed using the [IDSL.IPA](https://github.com/idslme/IDSL.IPA) workflow to acquire chromatographic information of the peaks (***m/z-RT***). When the chromatographic information of individual and aggregated aligned peaklists were generated using the [IDSL.IPA](https://github.com/idslme/IDSL.IPA) workflow, download the [IDSL.CSA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.CSA/main/CSA_parameters.xlsx) and select the parameters accordingly and then use this spreadsheet as the input for the IDSL.CSA workflow:

	library(IDSL.CSA)
	IDSL.CSA_workflow("Address of the CSA parameter spreadsheet")

## Quick Batch Example

Follow these steps for a quick case study (n = 33) [ST002263](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST002263&DataMode=AllData&ResultType=1) which has Thermo Q Exactive HF hybrid Orbitrap data collected in the HILIC-ESI-POS/NEG modes. 

1. Process raw mass spectrometry data and chromatographic information using the method described for [IDSL.IPA](https://github.com/idslme/IDSL.IPA#quick-batch-example)

2. The **Composite Spectra Analysis** requires 39 parameters distributed into 5 separate sections for a full scale analysis. For this study, use default parameter values presented in the [IDSL.CSA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.CSA/main/CSA_parameters.xlsx). Next, provide information for 
	
	2.1. Select **YES** for **PARAM0001** in the `Start` tab to only process **CSA** workflow.
	
	2.2. **CSA0005** for *HRMS data location address (MS1 level HRMS data)*
	
	2.3. **CSA0008** for *Address of the `peaklists` directory generated by the IDSL.IPA workflow*
	
	2.4. **CSA0009** for *Address of the `peak_alignment` directory generated by the IDSL.IPA workflow*
	
	2.5. **CSA0011** for *Output location (.msp files and EICs)*
	
	2.6. You may also increase the number of processing threads using **CSA0004** according to your computational power

3. Run this command in R/Rstudio console or terminal:

```
library(IDSL.CSA)
IDSL.CSA_workflow("Address of the CSA parameter spreadsheet")
```

4. You may parse the results at the address you provided for **CSA0011**.
	
	4.1. *CSA_MSP* includes ***.msp*** file
	
	4.2. *CSA_adduct_annotation* includes peaklists with potential adduct information
	
	4.3. *peak_alignment_subset* includes subsets of aligned peak tables for the major ions in each CSA cluster
	
	4.4. *aligned_spectra_table* includes information for the CSA aggregation on the aligned table

## [**Wiki**](https://github.com/idslme/IDSL.CSA/wiki)

1. [**CSA analysis by IDSL.CSA**](https://github.com/idslme/IDSL.CSA/wiki/CSA-analysis-by-IDSL.CSA)
2. [**DDA analysis by IDSL.CSA**](https://github.com/idslme/IDSL.CSA/wiki/DDA-analysis-by-IDSL.CSA)
3. [**DIA analysis by IDSL.CSA**](https://github.com/idslme/IDSL.CSA/wiki/DIA-analysis-by-IDSL.CSA)
4. [**Unique spectra aggregation**](https://github.com/idslme/IDSL.CSA/wiki/Unique-spectra-aggregation)

## Citation

[1] Fakouri Baygi, S., Kumar, Y. Barupal, D.K. [IDSL.CSA: Composite Spectra Analysis for Chemical Annotation of Untargeted Metabolomics Datasets](https://doi.org/10.1021/acs.analchem.3c00376). *Analytical Chemistry*, **2023**, *95(25)*, 9480�9487.

[2] Fakouri Baygi, S., Kumar, Y. Barupal, D.K. [IDSL. IPA characterizes the organic chemical space in untargeted LC/HRMS datasets](https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00120). *Journal of proteome research*, **2022**, *21(6)*, 1485-1494.