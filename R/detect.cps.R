#===========================================================================
# The main function that calls all other functions in fabisearch

#' detect.cps
#' @description This function takes a multivariate time series, \eqn{Y}, and returns the change points detected. Utilizes non-negative
#' factorization (NMF) to detect changes in clustering structure.
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#'
#' @param data Multivariate time series, \eqn{Y}, to be analyzed, should be in a matrix format with time points in rows and variables in columns
#' @param mindist Minimum distance between change points, by default is set to =35
#' @param nruns Number of runs to use for NMF function, by default is set to =50
#' @param nreps Number of repetitions for bootstrapping procedure, by default is set to =100
#' @param alpha Significance level cutoff for deciding whether to retain change points or not, by default is set to ="p-value" and returns the associated
#' p-value of the p-value test
#' @param rank Specifies whether to find rank using optimization procedure (specify "optimal") or use a predetermined rank (integer value, e.g., 4),
#' by default is set to ="optimal"
#' @param algtype Type of algorithm for NMF function, by default is set to ="brunet"
#'
#' @return Outputs a list where:\cr
#' \code{$rank} is the rank used for change point detection \cr
#' \code{$change.points} is a table of the change points detected where column "T" is the time of the change point and "stat.test" is the result of the t-test\cr
#' \code{$compute.time} is the compute time for the algorithm\cr
#' @export

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
