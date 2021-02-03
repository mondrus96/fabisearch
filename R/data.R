#===========================================================================
# Documentation for data sets

#' A simulated data set (see simulation 2 from Ondrus et al., 2021)
#'
#' A simulated data set (see simulation 2 from Ondrus et al., 2021). The data is generated from a multivariate Gaussian distribution with 2
#' clusters, where the correlation between nodes in the same cluster is 0.75 and between different clusters is 0.2 for the first 100 time
#' points. The vertex labels are randomly reshuflled for the second 100 time points. Hence, there is one change point at t=100.
#'
#' @format A matrix with 200 rows and 80 columns/variables.
"sim2"

#' Data from the NYU test-retest resting state fMRI data set
#'
#' A data matrix from the first stationary block of the second scan from the first subject from the NYU test-restest resting-state fMRI data set.
#'
#' @format A data matrix with 35 rows and 333 columns/variables, where each column corresponds to an ROI from the Gordon atlas.
#'
#' @source \url{https://www.nitrc.org/projects/nyu_trt}
"fmridata"

#' Adjacency matrix for the NYU test-restest resting-state fMRI data set
#'
#' The adjacency matrix calculated from the "fmridata" data set, using the Gordon atlas.
#'
#' @format A 333 * 333 matrix, where each entry takes a value 1 (0) if two nodes are (not) connected by an edge, using the Gordon atlas.
#'
#' @source \url{https://www.nitrc.org/projects/nyu_trt}
"adjmatrix"
