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

#' 333 ROI data from the NYU test-retest resting state fMRI data set
#'
#' A data matrix of the second scan from the first subject from the NYU test-restest resting-state fMRI data set.
#'
#' @format A data matrix with 197 rows and 333 columns/variables, where each column corresponds to an ROI from the Gordon atlas.
#'
#' @source \url{https://www.nitrc.org/projects/nyu_trt}
"gordfmri"

#' 90 ROI data from the NYU test-retest resting state fMRI data set
#'
#' A data matrix of the second scan from the first subject from the NYU test-restest resting-state fMRI data set.
#'
#' @format A data matrix with 197 rows and 90 columns/variables, where each column corresponds to an ROI from the AAL atlas.
#'
#' @source \url{https://www.nitrc.org/projects/nyu_trt}
"AALfmri"

#' Adjacency matrix for the NYU test-restest resting-state fMRI data set
#'
#' The adjacency matrix calculated from the "fmridata" data set, using the Gordon atlas.
#'
#' @format A 333 * 333 matrix, where each entry takes a value 1 (0) if two nodes are (not) connected by an edge, using the Gordon atlas.
#'
#' @source \url{https://www.nitrc.org/projects/nyu_trt}
"adjmatrix"

#' Gordon atlas coordinates
#'
#' A dataframe of the Gordon et al. (2016) atlas to use with the \code{\link{net.3dplot}()} function. Each row corresponds to a region of interest
#' (ROI) to be plotted using the Montreal Neurological Institute (MNI) space. The first column corresponds to the community labels as a string,
#' the second, third, and fourth columns correspond to the X, Y, and Z coordinates of the ROIs in MNI space. See Gordon et al. (2016)
#' <doi:10.1093/cercor/bhu239> for more details.
#'
#' @format A dataframe with 333 rows and 4 columns/variables.
#'
#' @source \url{https://doi.org/10.1093/cercor/bhu239}
"gordatlas"

#' Automated Anatomical Labeling (AAL) atlas coordinates
#'
#' A dataframe of the Automated Anatomical Labeling from the work of Tzourio-Mazoyer et al. (2002) atlas to use with the \code{\link{net.3dplot}()}
#' function. Each row corresponds to a region of interest (ROI) to be plotted using the Montreal Neurological Institute (MNI) space. The first
#' column corresponds to the community labels (in this atlas, there are none, therefore this column is filled with NA), the second, third, and
#' fourth columns correspond to the X, Y, and Z coordinates of the ROIs in MNI space. See Tzourio-Mazoyer et al. (2002) <doi:10.1006/nimg.2001.0978>
#' for more details.
#'
#' @format A dataframe with 90 rows and 4 columns/variables.
#'
#' @source \url{https://doi.org/10.1006/nimg.2001.0978}
"AALatlas"
