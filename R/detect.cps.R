#===========================================================================
# The main function that calls all other functions in fabisearch

#' detect.cps
#' @description This function takes a multivariate time series, \eqn{Y}, and returns the change points detected. Utilizes non-negative
#' factorization (NMF) to detect changes in clustering structure.
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#' @importFrom Rdpack reprompt
#'
#' @param Y A numerical matrix representing the multivariate time series, with the columns representing its components.
#' @param mindist A positive integer, by default is set to 35. It is used to define the minimum distance acceptable between detected change-points.
#' @param nruns A positive integer, by default is set to 50. It is used to define the number of runs in the NMF function.
#' @param nreps A positive integer, by default is set to 100. It is used to define the number of resamplings in the bootstrap procedure.
#' @param alpha A character string or a positive real number with default value equal to 0.05. If alpha = a positive integer value, say 0.05, then it is
#' used to define the significance level for inference on the change points. If alpha = "p-value", then the p-value calculated for inference on the change
#' points is returned.
#' @param rank A character string or a positive integer, which defines the rank be used in the optimization procedure to detect the change points.
#' If rank = “optimal”, which is also the default value, then the optimal rank is used. If rank = a positive integer value, say 4, then a predetermined
#' rank is used.
#' @param algtype A character string, which defines the algorithm to be used in the NMF function. By default is set to “brunet” - please see the "Algorithms"
#' section for more information on the available algorithms.
#'
#' @return A list where:\cr
#' \code{$rank} is the rank used for change point detection \cr
#' \code{$change.points} is a table of the change points detected where column "T" is the time of the change point and "stat.test" is the result of the t-test\cr
#' \code{$compute.time} is the compute time for the algorithm\cr
#' @export
#'
#' @section Algorithms:
#' All algorithms available are presented below, please note the "fabisearch" package builds upon the algorithms available in the "NMF" package
#' \insertCite{Gaujoux2010}{fabisearch}:\cr
#'
#' \code{"brunet"} - This algorithm is based on Kullback-Leibler divergence, from \insertCite{Brunet2004a}{fabisearch}. It uses multiplicative updates from
#' \insertCite{NIPS2000_1861}{fabisearch} with some small enhancements.\cr
#'
#' \code{"lee"} - This algorithm is based on Euclidian distances from \insertCite{NIPS2000_1861}{fabisearch}, and uses simple multiplicative updates.\cr
#'
#' \code{"ls-nmf"} - This is the least-squares NMF method from \insertCite{Wang2006}{fabisearch}. This algorithm uses an altered version of the Euclidian
#' distance based, multiplicative updates from \insertCite{NIPS2000_1861}{fabisearch}. It incorporates weights on each entry of the target matrix.\cr
#'
#' \code{"nsNMF"} - This is the nonsmooth NMF method from \insertCite{Pascual-Montano2006}{fabisearch}. This algorithm uses an altered version of the
#' Kullback-Leibler based, multiplicative updates from \insertCite{NIPS2000_1861}{fabisearch}. It includes an intermediate "smoothing" matrix, which is
#' intended to produce sparser factors.\cr
#'
#' \code{"offset"} - This is the offset NMF method from \insertCite{Badea2008}{fabisearch}. This algorithm uses an altered version of the Euclidian
#' distance based, multiplicative updates from \insertCite{NIPS2000_1861}{fabisearch}. It incorporates an intercept which is intended to reflect a common
#' pattern or baseline amongst components.\cr
#'
#' \code{"pe-nmf"} - This is the pattern-expression NMF from \insertCite{Zhang2008}{fabisearch}. This algorithm utilizes multiplicative updates to minimize
#' a Euclidian distance based objective function. It is further regularized such that the basis vectors effectively express patterns.\cr
#'
#' \code{"snmf/r","snmf/l"} - This is the alternating least-squares (ALS) approach from \insertCite{10.1093/bioinformatics/btm134}{fabisearch}. It uses the
#' non-negative, least-squares algorithm from \insertCite{VanBenthem2004}{fabisearch} to alternatingly estimate the basis and coefficent matrices. It
#' utilizes an Euclidian distance based objective function, and is regularized to promote either sparse basis ("snmf/l") or coefficent ("snmf/r") matrices \cr
#'
#' @examples
#' ## Estimating the change points for a multivariate dataset, Y, using the default settings
#' detect.cps(Y)
#'
#' ## Estimating the change points for a multivariate dataset, Y, with an alpha value of 0.05
#' detect.cps(Y, alpha = 0.05)
#'
#' ## Estimating the change points for a multivariate dataset, Y, with a prespecified rank of 6
#' detect.cps(Y, rank = 6)
#'
#' ## Estimating the change points for a multivariate dataset, Y, with non-default values
#' detect.cps(Y, mindist = 50, nruns = 100, nreps = 1000,
#'    alpha = 0.001, rank = 7, algtype = "ls-nmf")
#'
#' ## Example output from the detect.cps() function
#' $rank
#' [1] 2
#'
#' $change.points
#'     T stat.test
#' 1 130 0.2292454
#'
#' $compute.time
#' Time difference of 1.655591 mins
#'
#' @author Martin Ondrus, \email{mondrus@ualberta.ca}, Ivor Cribben, \email{cribben@ualberta.ca}
#' @references "Factorized Binary Search: a novel technique for change point detection in multivariate high-dimensional time series networks", Ondrus et al
#' (2021), preprint.
#'
#' \insertAllCited{}

detect.cps = function(Y, mindist = 35, nruns = 50, nreps = 100, alpha = 0.05, rank = "optimal", algtype = "brunet"){

  # Parallelization setup
  n.cores = detectCores()
  registerDoParallel(n.cores)
  print(paste("Number of cores:", n.cores))

  # Find T as the number of rows in the input matrix
  T = nrow(Y)

  # Initialize lower, upper, and define time series -> required for Recall function to work correctly
  lower = 1
  upper = T
  x = 1:T

  # Set seed
  set.seed(123)

  # Start timer for finding how long it took to compute
  compute.T.start = Sys.time()

  # Define the Y as a matrix
  Y = as.matrix(Y)

  # If rank has not been specified, then it must be found
  if (rank == "optimal"){
    n.rank = fabisearch:::opt.rank(Y, nruns, algtype)
  } else {
    n.rank = rank
    print(paste("User defined rank:", n.rank))
  }

  # Define split.index and optimal.ranks, need to define outside of function so "Recall" works inside the function
  split.index   = c()
  # Define the original splits
  orig.splits = fabisearch:::split_all(Y, split.index, lower, upper, x, mindist, nruns, n.rank, algtype)

  # Define the refitted splits
  refit.splits = fabisearch:::refit_splits(orig.splits, Y, T, x, nreps, n.rank, algtype)

  # Define the permutation distribution to compare with refitted splits
  perm.distr = fabisearch:::perm_distr(orig.splits, Y, T, x, nreps, n.rank, algtype)

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
