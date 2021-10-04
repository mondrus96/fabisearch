#===========================================================================
# The function that plots the adjacency matrix in 3D w/ brain, utilizes the Gordon atlas

#' 3D network plot of an adjacency matrix between pairs of change points
#' @description This function takes an adjacency matrix of a brain network and returns a 3D plot of it.
#'
#' @importFrom rgl par3d mfrow3d plot3d lines3d legend3d
#' @importFrom reshape2 melt
#'
#' @param A An adjacency matrix to be plotted (in numerical matrix format).
#' @param ROIs Either a vector of character strings specifying the communities to plot, or a vector of integers specifying which ROIs to plot by their ID. By
#' default it is set to NULL, and all communities and ROIs are plotted. Communities available for the Gordon atlas are: "Default", "SMhand", "SMmouth",
#' "Visual", "FrontoParietal", "Auditory", "None", "CinguloParietal", "RetrosplenialTemporal", "CinguloOperc", "VentralAttn", "Salience", and "DorsalAttn".
#' @param colors A vector of character strings specifying the hex codes for node colors to distinguish each community. By default, each community is given
#' a predefined, unique color.
#' @param coordROIs A dataframe of community tags and Montreal Neurological Institute (MNI) coordinates for regions of interest (ROIs) to plot, which is by
#' default set to \code{NULL} and uses the Gordon atlas. See ?gordon.atlas for an example using the Gordon atlas. Format of the dataframe is as follows: first column
#' is a string of community labels, then the subsequent three columns are the x, y, and z coordinates, respectively. See \code{AALatlas} and \code{gordatlas}
#' for examples.
#'
#' @return A 3D network plot of an adjacency matrix between pairs of change points, or for data without change points.
#' @export
#'
#' @examples
#' ## Plotting a 333 * 333 adjacency matrix "adjmatrix" with default settings
#' \donttest{net.3dplot(adjmatrix)}
#'
#' ## Plotting a 333 * 333 adjacency matrix "adjmatrix" with default colours but only
#' ## the "Visual", "FrontoParietal", and "Auditory" communities
#' comms = c("Visual", "FrontoParietal", "Auditory")
#' \donttest{net.3dplot(adjmatrix, ROIs = comms)}
#'
#' ## Plotting a 333 * 333 adjacency matrix "adjmatrix" with red, blue, and green
#' ## nodes to denote the "Default", "SMhand", and "Visual" communities
#' comms = c("Default", "SMhand", "Visual")
#' colrs = c("#FF0000", "#00FF00", "#0000FF")
#' \donttest{net.3dplot(adjmatrix, ROIs = comms, colors = colrs)}
#'
#' ## The default color palette is defined as follows
#' ## c("#D32F2F", "#303F9F", "#388E3C", "#FFEB3B", "#03A9F4", "#FF9800", "#673AB7",
#' ## "#CDDC39", "#9C27B0", "#795548", "#212121", "#009688", "#FFC0CB")
#'
#' @author Martin Ondrus, \email{mondrus@ualberta.ca}, Ivor Cribben, \email{cribben@ualberta.ca}
#' @references "Factorized Binary Search: a novel technique for change point detection in multivariate high-dimensional time series networks", Ondrus et al.
#' (2021), <arXiv:2103.06347>.

net.3dplot = function(A, ROIs = NULL, colors = NULL, coordROIs = NULL){

  # If colors are null, define a color palette
  if(is.null(colors)){
    colors = c("#D32F2F",
      "#303F9F",
      "#388E3C",
      "#FFEB3B",
      "#03A9F4",
      "#FF9800",
      "#673AB7",
      "#CDDC39",
      "#9C27B0",
      "#795548",
      "#212121",
      "#009688",
      "#FFC0CB")
  }

  # Get coordinates for the main brain frame
  coord = rbind(lcoord, rcoord)

  # Plot the main brain frame
  par3d(windowRect = c(0, 0, 800, 800),zoom=0.7)
  mfrow3d(1,1,sharedMouse = T)
  plot3d(coord,col='grey',size=0.1,alpha=0.7,
         box=F,axes=F,xlab='',ylab='',zlab='')

  # If input coordinates is not available, assume Gordon atlas
  if(is.null(coordROIs)){
    coordROIs = sys.gordatlas
  }

  # If ROIs is null, plot all ROIs
  if(is.null(ROIs)){
    ROIs = coordROIs
  } else if (class(ROIs) == "character"){
    ROIs = coordROIs[coordROIs[,1] %in% ROIs, ]
  } else if (is.numeric(ROIs)){
    ROIs = coordROIs[1:nrow(coordROIs) %in% ROIs, ]
  }

  # Prepare the adjacency matrix for plotting
  colnames(A) = rownames(A) = NULL
  A[!lower.tri(A)] = NA
  ma3d = melt(A, na.rm = TRUE)

  # Remove any edges which connect nodes to themselves, keep only entries where there is a connection
  ma3d = ma3d[!ma3d[,1] == ma3d[,2],]
  ma3d = ma3d[ma3d[,3] == 1,]

  # Loop through and plot specified communities
  for(i in 1:length(unique(ROIs[,1]))){
    # Define the current community
    curr.comm = unique(ROIs[,1])[i]

    # Find the coordinates of this community and the relevant nodes
    if(is.na(curr.comm)){
      coord.comm = ROIs[,2:4]
    } else if (!is.na(curr.comm)){
      coord.comm = ROIs[ROIs[,1] == curr.comm, 2:4]
    }

    # Plot these coordinates as nodes
    plot3d(coord.comm, col = colors[i], size=12, add=T)
  }

  # Narrow down ma3d to only include the edges for nodes that were specified
  ROI.vals = as.numeric(rownames(ROIs))
  ma3d = ma3d[ma3d[,1] %in% ROI.vals & ma3d[,2] %in% ROI.vals,]

  # Add a legend to the plot to denote the node communities, if communities exist
  if(length(unique(ROIs[,1])) > 1){
    communities = cbind(as.vector(unique(ROIs[,1])), colors[1:length(unique(ROIs[,1]))])
    legend3d("topright", pch = 16, legend = communities[,1], col = communities[,2], cex=1, inset=c(0.02))
  }

  # Plot the edges in ma3d
  for (i in 1:dim(ma3d)[1]) {
    lines3d(coordROIs[unlist(ma3d[i,1:2]), 2:4],
            size=2,
            add=T,
            col="black",
            alpha=0.4)
  }
}
