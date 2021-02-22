#===========================================================================
# The main function that calls all other functions in fabisearch

#' Multiple change point detection in the network (or clustering) structure of multivariate high-dimensional time series
#' @description This function detects multiple change points in the network (or clustering) structure of multivariate high-dimensional time series using
#' non-negative matrix factorization and a binary search.
#'
#' @param Y A numerical matrix representing the multivariate time series, with the columns representing its components.
#' @param mindist A positive integer with default value equal to 35. It is used to define the minimum distance acceptable between detected change points.
#' @param nruns A positive integer with default value equal to 50. It is used to define the number of runs in the NMF function.
#' @param nreps A positive integer with default value equal to 100. It is used to define the number of permutations for the statistical inference procedure.
#' @param alpha A character string or a positive real number with default value equal to 0.05. If alpha = a positive integer value, say 0.05, then it is
#' used to define the significance level for inference on the change points. If alpha = "p-value", then the p-value calculated for inference on the change
#' points is returned.
#' @param rank A character string or a positive integer, which defines the rank used in the optimization procedure to detect the change points.
#' If rank = "optimal", which is also the default value, then the optimal rank is used. If rank = a positive integer value, say 4, then a predetermined
#' rank is used.
#' @param algtype A character string, which defines the algorithm to be used in the NMF function. By default it is set to "brunet". See the "Algorithms" section of
#' \code{\link[NMF]{nmf}} for more information on the available algorithms.
#'
#' @return A list with the following components :\cr
#' \code{rank}: The rank used in the optimization procedure for change point detection.\cr
#' \code{change_points}: A table of the detected change points where column "T" is the time of the change point and "stat_test" is the result (either a boolean value if alpha = a positive real number, or the p-value if alpha = "p-value") of the t-test.\cr
#' \code{compute_time}: The computational time, saved as a "difftime" object.\cr
#' @export
#'
#' @examples
#' ## Change point detection for a multivariate data set, sim2, using the default settings
#' \donttest{detect.cps(sim2)}
#'
#' ## Change point detection for a multivariate data set, sim2, with an alpha value of 0.05
#' \donttest{detect.cps(sim2, alpha = 0.05)}
#'
#' ## Change point detection for a multivariate data set, sim2, with a prespecified rank of 6
#' \donttest{detect.cps(sim2, rank = 6)}
#'
#' ## Change point detection for a multivariate data set, sim2, with non-default values
#' \donttest{detect.cps(sim2, mindist = 50, nruns = 100, nreps = 1000,
#' alpha = 0.001, rank = 7, algtype = "snmf/l")}
#'
#' ## Example output from the detect.cps() function
#' \donttest{detect.cps(sim2, mindist = 50, nruns = 20)}
#'
#' # $rank
#' # [1] 5
#' #
#' # $change_points
#' #    T stat_test
#' # 1  99      TRUE
#' # 2 148     FALSE
#' #
#' # $compute_time
#' # Time difference of 15.8113 mins
#'
#' @author Martin Ondrus, \email{mondrus@ualberta.ca}, Ivor Cribben, \email{cribben@ualberta.ca}
#' @references "Factorized Binary Search: a novel technique for change point detection in multivariate high-dimensional time series networks", Ondrus et al.
#' (2021), preprint.

detect.cps = function(Y, mindist = 35, nruns = 50, nreps = 100, alpha = 0.05, rank = "optimal", algtype = "brunet"){

  # Find T as the number of rows in the input matrix
  T = nrow(Y)

  # Initialize lower, upper, and define time series -> required for Recall function to work correctly
  lower = 1
  upper = T
  x = 1:T

  # Start timer for finding how long it took to compute
  compute.T.start = Sys.time()

  # Define the Y as a matrix, put rownames into the matrix
  Y = as.matrix(Y)
  rownames(Y) = x

  # If rank has not been specified, then it must be found
  if (rank == "optimal"){
    n.rank = opt.rank(Y, nruns, algtype)
  } else {
    n.rank = rank
    print(paste("User defined rank:", n.rank))
  }

  # Define split.index and optimal.ranks, need to define outside of function so "Recall" works inside the function
  split.index   = c()
  # Define the original splits
  orig.splits = split_all(Y, split.index, lower, upper, x, mindist, nruns, n.rank, algtype)

  # Define the refitted splits
  refit.splits = refit_splits(orig.splits, Y, T, x, nreps, n.rank, algtype)

  # Define the permutation distribution to compare with refitted splits
  perm.distr = perm_distr(orig.splits, Y, T, x, nreps, n.rank, algtype)

  # Determine which splits are significant
  sign.splits = sign_splits(orig.splits, refit.splits, perm.distr, alpha)

  # End timer
  compute.T.end = Sys.time()

  # Define the variables for the final output
  cpt.time = difftime(compute.T.end, compute.T.start, units="mins")

  # Save the final output as a list and return from the function
  final.output = list(rank = n.rank, change_points = sign.splits, compute_time = cpt.time)
  return(final.output)
}
