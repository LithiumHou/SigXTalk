% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Crosstalk_analysis.R
\name{Calculate_Pathway_Fidelity}
\alias{Calculate_Pathway_Fidelity}
\title{Calculate the L-R-T pathway fidelity out of all the pathways with the same SSC and Target}
\usage{
Calculate_Pathway_Fidelity(CC_results, KeyTG, KeySSC)
}
\arguments{
\item{CC_results}{The dataframe of PRS values.}

\item{KeyTG}{The target to be focused on.}

\item{KeySSC}{The SSC to be focused on}
}
\value{
The array containing the fidelity of the signal, out of all the pathways with the same SSC and Target.
}
\description{
Calculate the L-R-T pathway fidelity out of all the pathways with the same SSC and Target
}
