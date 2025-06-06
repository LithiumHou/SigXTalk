% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Crosstalk_analysis.R
\name{Calculate_Pathway_Specificity}
\alias{Calculate_Pathway_Specificity}
\title{Calculate the L-R-T pathway specificity out of all the pathways with the same Signal and SSC}
\usage{
Calculate_Pathway_Specificity(CC_results, KeyRec, KeySSC)
}
\arguments{
\item{CC_results}{The dataframe of PRS values.}

\item{KeyRec}{The signal (receptor) to be focused on.}

\item{KeySSC}{The SSC to be focused on}
}
\value{
The array containing the fidelity of the signal, out of all the pathways with the same Signal and SSC.
}
\description{
Calculate the L-R-T pathway specificity out of all the pathways with the same Signal and SSC
}
