#===========================================================================
# The main function that calls all other functions FaBi Search

#' FaBiSearch
#' @description This is the main FaBiSearch function which takes a multivariate time series, \eqn{Y}, and returns change points detected. Utilizes non-negative
#' factorization (NMF) to detect changes in clustering structure.
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom utils write.csv
#' @importFrom parallel detectCores
#'
#' @param output.name What to call the output - note, NO NEED TO PUT FILE EXTENSION NAME (e.g. ".csv")
#' @param data Multivariate time series, \eqn{Y}, to be analyzed, should be in a matrix format with time points in rows and variables in columns
#' @param which.subj Which subjects to analyze, vector (e.g., 1:4). If left unspecified, then all subjects are considered
#' @param n.subj Total number of subjects in the dataset
#' @param min.dist Minimum distance between change points
#' @param n.rep Number of repetitions for bootstrapping procedure
#' @param alpha Significance level cutoff for deciding whether to retain change points or not (e.g., 0.05 or 0.01)
#' @param method.rank Specifies whether to find rank using optimization procedure (specify "optimal") or use a predetermined rank (integer value, e.g., 4)
#' @param n.runs Number of runs to use for NMF function
#' @param alg.type Type of algorithm for NMF function -> check ?nmf for details, under "method"
#' @param test.type Type of statistical test to run, "t-test" for simple t-test (using \code{alpha}, returns \code{TRUE/FALSE}), or "pval.t-test" to
#' return the probability value for the t-test (no cutoff)
#'
#' @return A \code{.csv} file of change points detected. Each row is a unique change point detected in the input multivariate time series, and contains the following columns:\cr
#' \code{subject}: subject number\cr
#' \code{rank.selection}: method used to determine rank\cr
#' \code{n.rank}: rank used for NMF and the change point estimation process\cr
#' \code{number.runs}: number of runs used for NMF and the change point estimation process\cr
#' \code{T.split}: time point where change point was detected\cr
#' \code{significant}: either the \eqn{p} value of the change point, or a boolean \code{TRUE/FALSE} whether this change point is kept using \code{alpha} as the
#' cutoff value\cr
#' \code{repetitions}: number of bootstrap repetitions for inference step\cr
#' \code{algorithm.type}: type of algorithm used by the NMF function\cr
#' \code{subject.compute.time}: time, in minutes, to compute change points for the corresponding subject\cr
#' @export
#'
#' @examples

FaBiSearch = function(output.name, data, which.subj, n.subj=NULL, T=NULL, min.dist, n.rep, alpha, method.rank, n.runs, alg.type, test.type){

  # output.name = what to call the output - note, NO NEED TO PUT FILE EXTENSION NAME (e.g. ".csv")
  # data        = data set to be analyzed
  # which.subj  = which subjects to analyze, vector (e.g. 1:4)
  # n.subj      = number of subjects in the data set
  # T           = number of time points in the data set
  # min.dist    = the minimum distance between change points
  # n.rep       = number of replications in the bootstrap procedures
  # alpha       = level of significance for statistical inference
  # method.rank = method type for determining rank value, can be an integer, vector, or string depending on
  # n.runs      = n.runs to try the NMF algorithm for, the higher the more exhaustive the search for an optimal matrix is
  # alg.type    = algorithm type -> check ?nmf for details, under "method"
  # test.type   = type of statistical test to run, "t-test" for simple t-test (using alpha, returns TRUE/FALSE), or "pval.t-test" to return the probability value for the t-test (no cutoff)

  # Parallelization setup
  n.cores = detectCores()
  registerDoParallel(n.cores)
  print(n.cores)

  # If one of n.subj or T is NULL, we can use the other to calculate it
  if (is.null(n.subj)){
    n.subj = nrow(data)/T
  }
  if (is.null(T)){
    T = nrow(data)/n.subj
  }

  # If which.subj is NULL, assume all subjects to be calculated
  if (is.null(which.subj)){
    which.subj = 1:n.subj
  }

  # Define data.list usine data_setup function
  data.list = data_setup(data, n.subj, T)

  # Initialize lower, upper, and define time series -> required for Recall function to work correctly
  lower = 1
  upper = T
  x = 1:T

  # Create the ALL.SPLITS variable which will hold all results
  ALL.SPLITS = list()

  # Main loop that goes through selected subjects
  for (j in which.subj) {
    # Setting seed using the subject iterator
    set.seed(123*j)

    # Start timer for finding how long it took to compute everything for this subject
    compute.T.start = Sys.time()

    # Define the current subject to be evaluated in this loop
    curr.subj = as.matrix(data.list[[j]])

    # If rank has not been specified, then it must be found
    if (method.rank == "optimal"){
      n.rank = optimal_rank(curr.subj, n.runs, alg.type)
    } else {
      n.rank = method.rank
      print(paste("User defined rank:", n.rank))
    }

    # Define split.index and optimal.ranks, need to define outside of function so "Recall" works inside the function
    split.index   = c()
    # Define the original splits
    orig.splits = split_all(curr.subj, split.index, lower, upper, x, min.dist, n.runs, n.rank, alg.type)

    # Define the refitted splits
    refit.splits = refit_splits(orig.splits, curr.subj, T, x, n.rep, n.rank, alg.type)

    # Define the permutation distribution to compare with refitted splits
    perm.distr = perm_distr(orig.splits, curr.subj, T, x, n.rep, n.rank, alg.type)

    # Determine which splits are significant
    sign.splits = sign_splits(orig.splits, refit.splits, perm.distr, alpha, test.type)

    # End timer
    compute.T.end = Sys.time()

    # Define rest of variables for final output, define and use num.rows to determine length
    num.rows = nrow(sign.splits)
    if(is.numeric(method.rank[1])){
      rank.selection = "user input"
    } else {
      rank.selection = method.rank[1]
    }

    rank.selection       = rep(rank.selection, num.rows)
    rank                 = rep(n.rank, num.rows)
    algorithm.type       = rep(alg.type, num.rows)
    number.runs          = rep(n.runs, num.rows)
    repetitions          = rep(n.rep, num.rows)
    subject              = rep(j, num.rows)
    subject.compute.time = rep(difftime(compute.T.end, compute.T.start, units="mins"), num.rows)

    # Save all results for current subject as final.output
    final.output = cbind(subject, rank.selection, n.rank, number.runs, sign.splits, repetitions, algorithm.type, subject.compute.time)

    # Save and print all current results
    ALL.SPLITS = rbind(ALL.SPLITS, final.output)
    print(ALL.SPLITS)

    # Save the output as a csv
    write.csv(ALL.SPLITS, paste(output.name, ".csv", sep=""))
  }
}
