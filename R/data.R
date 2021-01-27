#===========================================================================
# Documentation for datasets

#' Simulation 2 data set
#'
#' A data set of Simulation 2, which is generated from a multivariate Gaussian distribution with 2 clusters where the correlation between nodes in the same
#' cluster is 0.75 and between different clusters is 0.2. This data set contains a change point at t=100 where vertex labels are randomly reshuffled.
#'
#' @format A matrix with 200 rows and 80 columns/variables.
"sim2"

#' fMRI data set
#'
#' A data set of the first subject and their corresponding first stationary block from the second fMRI scan.
#'
#' @format A matrix with 35 rows and 333 columns/variables, where each column corresponds to an ROI from the Gordon atlas.
#'
#' @source \url{https://www.nitrc.org/projects/nyu_trt}
"fmridata"

#' Adjacency matrix of fMRI data set
#'
#' The adjacency matrix calculated from the "fmridata" data set, using the Gordon atlas.
#'
#' @format A matrix with 333 rows and 333 columns, where each entry denotes whether two nodes in the Gordon atlas are clustered together.
#'
#' @source \url{https://www.nitrc.org/projects/nyu_trt}
"adjmatrix"
