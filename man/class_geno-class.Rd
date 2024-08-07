% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_geno.R
\docType{class}
\name{class_geno-class}
\alias{class_geno-class}
\alias{class_geno}
\title{S4 class for marker genotype data}
\description{
S4 class for marker genotype data
}
\section{Slots}{

\describe{
\item{\code{ploidy}}{If mixed ploidy, then a vector equal to pop size; otherwise a single integer}

\item{\code{map}}{Marker map positions}

\item{\code{coeff}}{Coefficients of the marker effects (dim: indiv x marker)}

\item{\code{scale}}{Scaling factor between markers and indiv, vector of length equal to pop size}

\item{\code{G}}{Additive relationship matrix (from markers and potentially also pedigree)}

\item{\code{eigen.G}}{list of eigenvalues and eigenvectors}
}}

