\name{fabisearch}
\alias{fabisearch}
\title{FaBiSearch: a change-point detection method in the network structure of high-dimensional multivariate time series}

\description{The FaBiSearch package implements the Factorized Binary Search (FaBiSearch) methodology for the estimation of the number and location of multiple change-points in the network (or clustering) structure of multivariate high-dimensional time series. The method is motivated by the detection of change points in functional connectivity networks for functional magnetic resonance imaging (fMRI) data. FaBiSearch uses non-negative matrix factorization (NMF), an unsupervised dimension reduction technique and a new binary search algorithm to identify multiple change points.  It also requires minimal assumptions. The  main routines of the package are detect.cps, for change point detection, est.net, for estimating a network from stationary multivariate time series, plot.net, for plotting the estimated functional connectivity networks, and opt.rank, for finding the optimal rank to use for a given dataset. The functions have been extensively tested on high-dimensional simulated data and fMRI data.}

\author{Martin Ondrus, \email{mondrus@ualberta.ca}, Ivor Cribben, \email{cribben@ualberta.ca}}

\section{References}{"Factorized Binary Search: a novel technique for change point detection in multivariate high-dimensional time series networks", Ondrus et al (2021), preprint.}

\section{See Also}{\code{\link{detect.cps}, \link{est.net}, \link{plot.net}, \link{opt.rank}}}
