#===========================================================================
# The main function that calls all other functions in fabisearch

#' detect.cps
#' @description This function takes a multivariate time series, \eqn{Y}, and returns the change points detected. Utilizes non-negative
#' factorization (NMF) to detect changes in clustering structure.
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#'
#' @param data A matrix with time points in rows and variables in columns. This input is the multivariate time series, \eqn{Y}, to be analyzed.
#' @param mindist A positive integer, by default is set to 35. Defines the minimum distance between change points.
#' @param nruns A positive integer, by default is set to 5. Defines the number of runs to use for NMF function.
#' @param nreps A positive integer, by default is set to 100. Defines the number of repetitions for bootstrapping procedure.
#' @param alpha A positive real number denoting the significance level cutoff to determine whether to keep change points or not. By default is set to
#' "p-value" and returns the associated p-value of the t-test.
#' @param rank A positive integer denoting the value of rank to use in the algorithm. By default is set to "optimal" so the function finds the optimal rank
#' to use.
#' @param algtype A character string denoting the type of algorithm for NMF function, by default is set to "brunet".
#'
#' @return Outputs a list where:\cr
#' \code{$rank} is the rank used for change point detection \cr
#' \code{$change.points} is a table of the change points detected where column "T" is the time of the change point and "stat.test" is the result of the t-test\cr
#' \code{$compute.time} is the compute time for the algorithm\cr
#' @export
#'
#' @examples
#' ## Estimating the change points for a multivariate dataset, "data", using the default settings
#' detect.cps(data)
#'
#' ## Estimating the change points for a multivariate dataset, "data", with an alpha value of 0.05
#' detect.cps(data, alpha = 0.05)
#'
#' ## Estimating the change points for a multivariate dataset, "data", with a prespecified rank of 6
#' detect.cps(data, rank = 6)
#'
#' ## Estimating the change points for a multivariate dataset, "data", with non-default values
#' detect.cps(data, mindist = 50, nruns = 100, nreps = 1000, alpha = 0.001, rank = 7, algtype = "ls-nmf")
#'
#' ## Example output from the detect.cps() function
#' $rank
#' [1] 2
#'
#' $change.points
#' T stat.test
#' 1 130 0.2292454
#'
#' $compute.time
#' Time difference of 1.655591 mins
#' @author Martin Ondrus, \email{mondrus@ualberta.ca}, Ivor Cribben, \email{cribben@ualberta.ca}
#' @references "Factorized Binary Search: a novel technique for change point detection in multivariate high-dimensional time series networks", Ondrus et al
#' (2021), preprint.

detect.cps = function(data, mindist = 35, nruns = 50, nreps = 100, alpha = "p-value", rank = "optimal", algtype = "brunet"){

  # Parallelization setup
  n.cores = detectCores()
  registerDoParallel(n.cores)
  print(paste("Number of cores:", n.cores))

  # Find T as the number of rows in the input matrix
  T = nrow(data)

  # Initialize lower, upper, and define time series -> required for Recall function to work correctly
  lower = 1
  upper = T
  x = 1:T

  # Set seed
  set.seed(123)

  # Start timer for finding how long it took to compute
  compute.T.start = Sys.time()

  # Define the data as a matrix
  data = as.matrix(data)

  # If rank has not been specified, then it must be found
  if (rank == "optimal"){
    n.rank = fabisearch:::opt.rank(data, nruns, algtype)
  } else {
    n.rank = rank
    print(paste("User defined rank:", n.rank))
  }

  # Define split.index and optimal.ranks, need to define outside of function so "Recall" works inside the function
  split.index   = c()
  # Define the original splits
  orig.splits = fabisearch:::split_all(data, split.index, lower, upper, x, mindist, nruns, n.rank, algtype)

  # Define the refitted splits
  refit.splits = fabisearch:::refit_splits(orig.splits, data, T, x, nreps, n.rank, algtype)

  # Define the permutation distribution to compare with refitted splits
  perm.distr = fabisearch:::perm_distr(orig.splits, data, T, x, nreps, n.rank, algtype)

  # Determine which splits are significant
  sign.splits = fabisearch:::sign_splits(orig.splits, refit.splits, perm.distr, alpha)

  # End timer
  compute.T.end = Sys.time()

  # Define the variables for the final output
  cpt.time = difftime(compute.T.end, compute.T.start, units="mins")

  # Save the final output as a list and return from the function
  final.output = list(rank = n.rank, change.points = sign.splits, compute.time = cpt.time)
  return(final.output)
}
