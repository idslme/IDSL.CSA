\name{CSA_adductAnnotator}
\alias{CSA_adductAnnotator}
\title{
CSA Adduct Annotator
}
\description{
This function updates IDSL.IPA peaklists with IDSL.CSA grouping
}
\usage{
CSA_adductAnnotator(IPApeakList, CSA_peaklist, massError)
}
\arguments{
  \item{IPApeakList}{
IDSL.IPA peaklist
}
  \item{CSA_peaklist}{
A dataframe peaklist of co-detected CSA analysis.
}
  \item{massError}{
Mass accuracy in Da
}
}
\value{
IDSL.IPA peaklists with IDSL.CSA grouping
}