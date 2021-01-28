\name{fabisearch}
\alias{fabisearch}
\title{FaBiSearch: a Change Point Detection Method in the Network Structure of Multivariate High-Dimensional Time Series}

\description{The \code{fabisearch} package implements the Factorized Binary Search (FaBiSearch) methodology for the estimation of the number and location of multiple change points in the network (or clustering) structure of multivariate high-dimensional time series. The method is motivated by the detection of change points in functional connectivity networks for functional magnetic resonance imaging (fMRI) data. FaBiSearch uses non-negative matrix factorization (NMF), an unsupervised dimension reduction technique and a new binary search algorithm to identify multiple change points.  It also requires minimal assumptions. The  main routines of the package are \link{detect.cps}, for multiple change point detection, \link{est.net}, for estimating a network from stationary multivariate time series, \link{net.3dplot}, for plotting the estimated functional connectivity networks, and \link{opt.rank}, for finding the optimal rank to use for a given data set. The functions have been extensively tested on high-dimensional simulated data and fMRI data.}

\author{Martin Ondrus, \email{mondrus@ualberta.ca}, Ivor Cribben, \email{cribben@ualberta.ca}}

\section{References}{"Factorized Binary Search: a novel technique for change point detection in multivariate high-dimensional time series networks", Ondrus et al (2021), preprint.}

\section{See Also}{\code{\link{detect.cps}, \link{est.net}, \link{net.3dplot}, \link{opt.rank}}}
