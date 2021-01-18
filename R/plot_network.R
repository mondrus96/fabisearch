#===========================================================================
# The function that plots the adjacency matrix in 3D w/ brain, utilizes the Gordon atlas

#' plot_network
#' @description This function uses a Gordon atlas defined adjacency matrix and returns a 3D plot of the estimated stationary network of this adjacency matrix
#'
#' @importFrom rgl par3d mfrow3d plot3d lines3d
#' @importFrom reshape2 melt
#' @importFrom utils read.csv
#'
#' @param adj.matrix A 333 by 333 adjacency matrix
#' @param communities A vector of strings related to communities to plot - by default, all communities are plotted - (any combination of "Default", "SMhand", "SMmouth", "Visual", "FrontoParietal", "Auditory", "None", "CinguloParietal", "RetrosplenialTemporal", "CinguloOperc", "VentralAttn", "Salience", "DorsalAttn)
#' @param colors A vector of hex codes for node colors to distinguish each community - by default, each community is given a predefined color
#'
#' @return a 3D plot of the estimated stationary network from the adjacency matrix
#' @export
#'
#' @examples

plot_network = function(adj.matrix, communities=NULL, colors=NULL){

  # adj.matrix  = input adjacency matrix to plot from Gordon atlas
  # communities = a vector of communities from the Gordon atlas to include in the plot
  # colors      = a vector of colors of the selected communities - must be of equal length to communities

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
  lcoord = FaBiSearch::lcoord
  rcoord = FaBiSearch::rcoord
  coord = rbind(lcoord, rcoord)

  # Plot the main brain frame
  par3d(windowRect = c(0, 0, 800, 800),zoom=0.7)
  mfrow3d(1,1,sharedMouse = T)
  plot3d(coord,col='grey',size=0.1,alpha=0.7,
         box=F,axes=F,xlab='',ylab='',zlab='',
         mar = c(0, 0, 0, 0))

  # Get the coordinates for the Gordon atlas regions
  gordon.atlas = FaBiSearch::gordon.atlas
  coord333 = as.matrix(gordon.atlas[,c('x.mni','y.mni','z.mni')])
  rownames(coord333) = paste0(1:333)
  name.netwk = as.matrix(gordon.atlas[,c('Community')])

  # If communities is null, plot all communities
  if(is.null(communities)){
    communities = unique(name.netwk)
  }

  # Prepare the adjacency matrix for plotting
  adj.matrix[!lower.tri(adj.matrix)] = NA
  ma3d = melt(adj.matrix, na.rm = TRUE)

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
