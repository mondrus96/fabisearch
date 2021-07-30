\name{fabisearch}
\alias{fabisearch}
\title{Change Point Detection in High-Dimensional Time Series Networks}

\description{Implementation of the Factorized Binary Search (FaBiSearch) methodology for the estimation of the number and the location of multiple change points in the network (or clustering) structure of multivariate high-dimensional time series. The method is motivated by the detection of change points in functional connectivity networks for functional magnetic resonance imaging (fMRI) data. FaBiSearch uses non-negative matrix factorization (NMF), an unsupervised dimension reduction technique, and a new binary search algorithm to identify multiple change points.  It requires minimal assumptions. Lastly, we provide interactive, 3-dimensional, brain-specific network visualization capability in a flexible, stand-alone function. This function can be conveniently used with any node coordinate atlas, and nodes can be color coded according to community membership, if applicable. The output is an elegantly displayed network laid over a cortical surface, which can berotated in the 3-dimensional space. The main routines of the package are \code{\link{detect.cps}()}, for multiple change point detection, \code{\link{est.net}()}, for estimating a network between stationary multivariate time series, \code{\link{net.3dplot}()}, for plotting the estimated functional connectivity networks, and \code{\link{opt.rank}()}, for finding the optimal rank in NMF for a given data set. The functions have been extensively tested on simulated multivariate high-dimensional time series data and fMRI data. For details on the FaBiSearch methodology, please see Ondrus et al. (2021) <arXiv:2103.06347>. For a more detailed explanation and applied examples of the \code{fabisearch} package, please see Ondrus and Cribben (2021), preprint.}

\value{No return value, called for side effects}

\references{"Factorized Binary Search: a novel technique for change point detection in multivariate high-dimensional time series networks", Ondrus et al (2021), <arXiv:2103.06347>.}

\author{Martin Ondrus, \email{mondrus@ualberta.ca}, Ivor Cribben, \email{cribben@ualberta.ca}}

\section{See Also}{\code{\link{detect.cps}, \link{est.net}, \link{net.3dplot}, \link{opt.rank}}}
