#===========================================================================
# This function splits the time course into partitions exhaustively based on NMF

#' @importFrom NMF nmf

split_all = function(data, split.index, lower, upper, x, min.dist, n.runs, n.rank, alg.type){

  # curr.subj     = data for current subject
  # split.index   = for saving results
  # lower, upper  = lower limit of split and upper limit of split respectively
  # x             = time series of data
  # min.dist      = smallest gap between splits
  # n.runs        = n.runs to try the NMF algorithm for, the higher the more exhaustive the search for an optimal matrix is
  # n.rank        = value of rank to use in the NMF function
  # alg.type      = algorithm type -> check ?nmf for details, under "method"

  # Define the lower and upper boundaries for split
  low.s = lower + min.dist
  upp.s = upper - min.dist

  # Check if the lower boundary is below the upper
  if (low.s <= upp.s) {
    # Define the indices
    indices = low.s:upp.s

    # While loop to go through the binary search
    while (length(indices) > 2){
      print(paste(min(indices) - 1, ":", max(indices)))
      # Define how to split the time series
      block.size = length(indices)/2

      # Define the left and right block
      if (block.size%%1==0){
        # Define the left block and right block
        l.block = (min(indices) - min.dist):indices[block.size]
        r.block = indices[block.size]:(max(indices) + min.dist - 1)
      } else if (block.size%%1!=0){
        # Round block.size up
        block.size = ceiling(block.size)

        # Define the left block and right block
        l.block = (min(indices) - min.dist):indices[block.size]
        r.block = indices[block.size]:(max(indices) + min.dist)
      }

      # Calculate the left and right NMF
      l.NMF = nmf(data[l.block,], rank = n.rank, nrun = n.runs, method = alg.type)
      r.NMF = nmf(data[r.block,], rank = n.rank, nrun = n.runs, method = alg.type)

      # Calculate the metrics for the left and right
      l.metric = l.NMF@residuals
      r.metric = r.NMF@residuals

      # Redefine indices based on which side has the higher metric
      if (l.metric > r.metric){
        indices = min(indices):indices[block.size]
      } else if (r.metric > l.metric){
        indices = indices[block.size]:max(indices)
      }
    }

    # For when indices is less than or equal to two
    if (length(indices) == 2){
      print(paste(min(indices) - 1, ":", max(indices)))
      block.size = length(indices)/2

      # Define the left block and right block
      l.block = (min(indices) - min.dist):indices[block.size]
      r.block = indices[block.size]:(max(indices) + min.dist - 1)

      # Calculate the left and right NMF
      l.NMF = nmf(data[l.block,], rank = n.rank, nrun = n.runs, method = alg.type)
      r.NMF = nmf(data[r.block,], rank = n.rank, nrun = n.runs, method = alg.type)

      # Calculate the metrics for the left and right
      l.metric = l.NMF@residuals
      r.metric = r.NMF@residuals

      # Redefine indices based on which side has the greater metric
      if (l.metric > r.metric){
        indices = indices[1]
      } else if (r.metric > l.metric){
        indices = indices[2]
      }
    }
    if (length(indices) == 1){
      print(paste(indices - 1, ":", indices))
      # Define the left block and right block
      l.block = (indices - min.dist):indices
      r.block = indices:(indices + min.dist)

      # Calculate the left and right NMF
      l.NMF = nmf(data[l.block,], rank = n.rank, nrun = n.runs, method = alg.type)
      r.NMF = nmf(data[r.block,], rank = n.rank, nrun = n.runs, method = alg.type)

      # Calculate the metrics for the left and right
      l.metric = l.NMF@residuals
      r.metric = r.NMF@residuals

      # Find the change point
      if (l.metric > r.metric){
        T.split = indices - 1
      } else if (r.metric > l.metric){
        T.split = indices
      }

      # Define the start and end points of the block to be evaluated
      T.start = T.split - min.dist + 1
      T.end   = T.split + min.dist

      # Fit the unsplit data as orig.loss and the split data as split.loss
      orig.loss  = nmf(data[which(x>=T.start & x<=T.end),], rank = n.rank, nrun = n.runs, method = alg.type)@residuals
      split.loss = nmf(data[which(x>=T.start & x<=T.split),], rank = n.rank, nrun = n.runs, method = alg.type)@residuals +
        nmf(data[which(x>=(T.split+1) & x<=T.end),], rank = n.rank, nrun = n.runs, method = alg.type)@residuals

      # Calculate the reduction in loss by splitting
      chg.loss = split.loss - orig.loss
    }

    # Return the change point found
    print(paste("Change Point At:",T.split,", Delta Loss:",chg.loss))
    split.index = rbind(split.index, data.frame(T.split, chg.loss), make.row.names = FALSE)

    # Apply the function exhaustively
    split.index = Recall(data, split.index, lower, T.split, x, min.dist, n.runs, n.rank, alg.type)
    split.index = Recall(data, split.index, T.split + 1, upper, x, min.dist, n.runs, n.rank, alg.type)
  }
  return(split.index)
}
