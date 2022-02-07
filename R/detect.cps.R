#===========================================================================
# The main function that calls all other functions in fabisearch

#' Multiple change point detection in the network (or clustering) structure of multivariate high-dimensional time series
#' @description This function detects multiple change points in the network (or clustering) structure of multivariate high-dimensional time series using
#' non-negative matrix factorization and a binary search.
#'
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom parallel detectCores
#' @importFrom pkgmaker isCHECK
#'
#' @param Y An input multivariate time series in matrix format, with variables organized in columns and time points in rows. All entries in Y must be positive.
#' @param mindist A positive integer with default value equal to 35. It is used to define the minimum distance acceptable between detected change points.
#' @param nruns A positive integer with default value equal to 50. It is used to define the number of runs in the NMF function.
#' @param nreps A positive integer with default value equal to 100. It is used to define the number of permutations for the statistical inference procedure.
#' @param alpha A positive real number with default value set to NULL. When alpha = NULL, then the p-value calculated for inference on the change
#' points is returned. If alpha = a positive integer value, say 0.05, then it is used to define the significance level for inference on the change points.
#' @param rank A positive integer, which defines the rank used in the optimization procedure to detect the change points. If rank = NULL, which is also the
#' default value, then the optimal rank is computed. If rank = a positive integer value, say 4, then a predetermined rank is used.
#' @param algtype A character string, which defines the algorithm to be used in the NMF function. By default it is set to "brunet". See the "Algorithms" section of
#' \code{\link[NMF]{nmf}} for more information on the available algorithms.
#' @param testtype A character string, which defines the type of statistical test to use during the inference procedure. By default it is set to "t-test". The
#' other options are "ks" and "wilcox" which correspond to the Kolmogorov-Smirnov and Wilcoxon tests, respectively.
#'
#' @return A list with the following components :\cr
#' \code{rank}: The rank used in the optimization procedure for change point detection.\cr
#' \code{change_points}: A table of the detected change points where column "T" is the time of the change point and "stat_test" is the result (either a boolean value if alpha = a positive real number, or the p-value if alpha = NULL) of the t-test.\cr
#' \code{compute_time}: The computational time, saved as a "difftime" object.\cr
#' @export
#'
#' @examples
#' \donttest{
#' ## Change point detection for a multivariate data set, sim2, using settings:
#' ## rank = 3, mindist = 99, nruns = 2, and nreps = 2
#' set.seed(123)
#' detect.cps(sim2, rank = 3, mindist = 99, nruns = 2, nreps = 2)
#' }
#'
#' # $rank
#' # [1] 3
#' #
#' # $change_points
#' #     T stat_test
#' # 1 101 0.3867274
#' #
#' # $compute_time
#' # Time difference of 0.741534 mins
#'
#' @author Martin Ondrus, \email{mondrus@ualberta.ca}, Ivor Cribben, \email{cribben@ualberta.ca}
#' @references "Factorized Binary Search: a novel technique for change point detection in multivariate high-dimensional time series networks", Ondrus et al.
#' (2021), <arXiv:2103.06347>.

detect.cps = function(Y, mindist = 35, nruns = 50, nreps = 100, alpha = NULL, rank = NULL, algtype = "brunet", testtype = "t-test"){

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
  if (is.null(rank)){
    n.rank = opt.rank(Y, nruns, algtype)
  } else {
    n.rank = rank
    print(paste("User defined rank:", n.rank))
  }

  # Define split.index and optimal.ranks, need to define outside of function so "Recall" works inside the function
  split.index   = c()
  # Define the original splits
  orig.splits = split_all(Y, split.index, lower, upper, x, mindist, nruns, n.rank, algtype)

  # Check whether any of the change points found have a negative change in loss, otherwise do not run rest of procedures
  if(any(orig.splits$chg.loss < 0)){
    # Register parallel backend, use maximum number of cores unless running a CRAN check
    if(isCHECK()){
      registerDoParallel(min(detectCores(), 2))
    } else if (!isCHECK()){
      registerDoParallel(detectCores())
    }

    # Define the refitted splits
    refit.splits = refit_splits(orig.splits, Y, T, x, nreps, n.rank, algtype)

    # Define the permutation distribution to compare with refitted splits
    perm.distr = perm_distr(orig.splits, Y, T, x, nreps, n.rank, algtype)

    # Determine which splits are significant
    sign.splits = sign_splits(orig.splits, refit.splits, perm.distr, alpha, testtype)
  } else {
    # Initialize an empty dataframe
    sign.splits = data.frame(T = double(), stat_test = logical())
  }

  # End timer
  compute.T.end = Sys.time()

  # Define the variables for the final output
  cpt.time = difftime(compute.T.end, compute.T.start, units="mins")

  # Close the parallel backend
  stopImplicitCluster()

  # Save the final output as a list and return from the function
  final.output = list(rank = n.rank, change_points = sign.splits, compute_time = cpt.time)
  return(final.output)
}
