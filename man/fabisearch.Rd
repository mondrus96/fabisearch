\name{fabisearch}
\alias{fabisearch}
\title{Change Point Detection in High-Dimensional Time Series Networks}

\description{Implementation of the Factorized Binary Search (FaBiSearch) methodology for the estimation of the number and location of multiple change points in the network (or clustering) structure of multivariate high-dimensional time series. The method is motivated by the detection of change points in functional connectivity networks for functional magnetic resonance imaging (fMRI) data. FaBiSearch uses non-negative matrix factorization (NMF), an unsupervised dimension reduction technique, and a new binary search algorithm to identify multiple change points.  It also requires minimal assumptions. The  main routines of the package are \code{\link{detect.cps}()}, for multiple change point detection, \code{\link{est.net}()}, for estimating a network between stationary multivariate time series, \code{\link{net.3dplot}()}, for plotting the estimated functional connectivity networks, and \code{\link{opt.rank}()}, for finding the optimal rank in NMF for a given data set. The functions have been extensively tested on simulated multivariate high-dimensional time series data and fMRI data. For details on the FaBiSearch methodology, please see Ondrus et al. (2021).}

\value{No return value, called for side effects}

\references{"Factorized Binary Search: a novel technique for change point detection in multivariate high-dimensional time series networks", Ondrus et al (2021), preprint.}

\author{Martin Ondrus, \email{mondrus@ualberta.ca}, Ivor Cribben, \email{cribben@ualberta.ca}}

\section{See Also}{\code{\link{detect.cps}, \link{est.net}, \link{net.3dplot}, \link{opt.rank}}}
