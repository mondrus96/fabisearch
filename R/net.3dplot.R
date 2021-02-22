#===========================================================================
# The function that plots the adjacency matrix in 3D w/ brain, utilizes the Gordon atlas

#' 3D network plot of an adjacency matrix between pairs of change points
#' @description This function takes a Gordon atlas defined adjacency matrix and returns a 3D plot of it.
#'
#' @importFrom rgl par3d mfrow3d plot3d lines3d legend3d
#' @importFrom reshape2 melt
#'
#' @param adjmatrix A numerical matrix of dimension 333*333 (Gordon atlas). This is the adjacency matrix to be plotted.
#' @param communities A vector of character strings specifying the communities to plot. By default, all communities are plotted. Communities available are:
#' "Default", "SMhand", "SMmouth", "Visual", "FrontoParietal", "Auditory", "None", "CinguloParietal", "RetrosplenialTemporal", "CinguloOperc",
#' "VentralAttn", "Salience", and "DorsalAttn".
#' @param colors A vector of character strings specifying the hex codes for node colors to distinguish each community. By default, each community is given
#' a predefined, unique color.
#'
#' @return A 3D network plot of an adjacency matrix between pairs of change points.
#' @export
#'
#' @examples
#' ## Plotting a 333 * 333 adjacency matrix "adjmatrix" with default settings
#' \donttest{net.3dplot(adjmatrix)}
#'
#' ## Plotting a 333 * 333 adjacency matrix "adjmatrix" with default colours but only
#' ## the "Visual", "FrontoParietal", and "Auditory" communities
#' comms = c("Visual", "FrontoParietal", "Auditory")
#' \donttest{net.3dplot(adjmatrix, communities = comms)}
#'
#' ## Plotting a 333 * 333 adjacency matrix "adjmatrix" with red, blue, and green
#' ## nodes to denote the "Default", "SMhand", and "Visual" communities
#' comms = c("Default", "SMhand", "Visual")
#' colrs = c("#FF0000", "#00FF00", "#0000FF")
#' \donttest{net.3dplot(adjmatrix, communities = comms, colors = colrs)}
#'
#' ## The default color palette is defined as follows
#' ## c("#D32F2F", "#303F9F", "#388E3C", "#FFEB3B", "#03A9F4", "#FF9800", "#673AB7",
#' ## "#CDDC39", "#9C27B0", "#795548", "#212121", "#009688", "#FFC0CB")
#'
#' @author Martin Ondrus, \email{mondrus@ualberta.ca}, Ivor Cribben, \email{cribben@ualberta.ca}
#' @references "Factorized Binary Search: a novel technique for change point detection in multivariate high-dimensional time series networks", Ondrus et al.
#' (2021), preprint.

net.3dplot = function(adjmatrix, communities = NULL, colors = NULL){

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
         box=F,axes=F,xlab='',ylab='',zlab='',
         mar = c(0, 0, 0, 0))

  # Get the coordinates for the Gordon atlas regions
  coord333 = as.matrix(gordon.atlas[,c('x.mni','y.mni','z.mni')])
  rownames(coord333) = paste0(1:333)
  name.netwk = as.matrix(gordon.atlas[,c('Community')])

  # If communities is null, plot all communities
  if(is.null(communities)){
    communities = unique(name.netwk)
  }

  # Prepare the adjacency matrix for plotting
  colnames(adjmatrix) = rownames(adjmatrix) = NULL
  adjmatrix[!lower.tri(adjmatrix)] = NA
  ma3d = melt(adjmatrix, na.rm = TRUE)

  # Remove any edges which connect nodes to themselves, keep only entries where there is a connection
  ma3d = ma3d[!ma3d[,1] == ma3d[,2],]
  ma3d = ma3d[ma3d[,3] == 1,]

  # Loop through and plot specified communities
  for(i in 1:length(communities)){
    # Define the current community
    curr.netwk = communities[i]

    # Find the coordinates of this community
    coord.comm = coord333[name.netwk == curr.netwk,]

    # Plot these coordinates as nodes
    plot3d(coord.comm, col = colors[i], size=12, add=T)
  }

  # Narrow down ma3d to only include the edges for communities that were specified
  ROI.vals = (1:333)[name.netwk %in% communities]
  ma3d = ma3d[ma3d[,1] %in% ROI.vals & ma3d[,2] %in% ROI.vals,]

  # Add a legend to the plot to denote the node communities
  communities = cbind(communities, colors[1:length(communities)])
  legend3d("topright", pch = 16, legend = communities[,1], col = communities[,2], cex=1, inset=c(0.02))

  # Plot the edges in ma3d
  for (i in 1:dim(ma3d)[1]) {
    lines3d(coord333[unlist(ma3d[i,1:2]),],
            size=2,
            add=T,
            col="black",
            alpha=0.4)
  }
}
