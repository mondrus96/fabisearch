#===========================================================================
# Recalculates the NMF basis matrix difference between each pair of splits

#' refit_splits
#' @description This function is embedded in the main FaBiSearch function. It serves to recalculate the loss at each candidate change point considering the fact that
#' multiple change points might have been detected.
#'
#' @importFrom NMF nmf
#'
#' @param orig.splits Matrix of candidate change points (\code{T.split}) and the change in loss (\code{chg.loss})
#' @param curr.subj Multivariate time series data from the current subject
#' @param x Index set of time series data
#' @param n.rep Number of repetitions for bootstrapping procedure
#' @param n.rank Rank value for NMF function
#' @param alg.type Type of algorithm for NMF function -> check ?nmf for details, under "method"
#'
#' @return List where each element is a unique change point and contains a vector of loss values from the refitted NMF calculations
#' @export
#'
#' @examples

refit_splits = function(orig.splits, curr.subj, T, x, n.rep, n.rank, alg.type){

  # orig.splits = splitting times for subjects, contains the ranks for left and right used as well
  # curr.subj   = data for current subject
  # x           = time series of time
  # T           = number of time points in the data set
  # n.rep       = number of times to run the NMF algorithm for statistical inference
  # n.rank      = value of rank to use in the NMF function
  # alg.type    = algorithm type -> check ?nmf for details, under "method"

  # Output will be saved as refit.results
  refit.results = list()

  # Find the splits where the loss metric is reduced
  reduced.splits = orig.splits[orig.splits$chg.loss < 0, ]

  # Order the vector from smallest to largest, add 0 to beggining and T to end
  split.times = c(0,sort(reduced.splits$T.split),T)

  for (ij in 1:(length(split.times)-2)){

    # Print the current split being evaluated
    print(paste("Refitting split at", split.times[ij+1]))

    # Define the lower and upper parts of the split, as well as the ind.interval of time points
    lower1  = split.times[ij]
    T.split = split.times[ij+1]
    upper1  = split.times[ij+2]

    # Loop through to find the sum of the left and right sides for each run
    curr.results = c()
    for(i in 1:n.rep){
      # Fit NMF to the left and right sides
      l.NMF = nmf(curr.subj[which(x<=T.split & x>lower1),], rank=n.rank, method=alg.type)
      r.NMF = nmf(curr.subj[which(x<=upper1 & x>T.split),], rank=n.rank, method=alg.type)

      curr.results = c(curr.results, sum(l.NMF@residuals + r.NMF@residuals), use.names = FALSE)
    }

    # Compile results into refit.results matrix
    refit.results[[ij]] = curr.results
  }
  return(refit.results)
}
