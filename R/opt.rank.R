#===========================================================================
# Helper function to find the optimal rank of a data set

#' Finds the optimal rank for non-negative matrix factorization (NMF)
#' @description This function finds the optimal rank for non-negative matrix factorization (NMF).
#'
#' @importFrom NMF nmf
#'
#' @param Y An input multivariate time series in matrix format, with variables organized in columns and time points in rows. All entries in Y must be positive.
#' @param nruns A positive integer with default value equal to 50. It is used to define the number of runs in the NMF function.
#' @param algtype A character string, which defines the algorithm to be used in the NMF function. By default it is set to "brunet". See the "Algorithms" section of
#' \code{\link[NMF]{nmf}} for more information on the available algorithms.
#'
#' @return A positive integer representing the optimal rank.
#' @export
#'
#' @examples
#' \donttest{
#' ## Finding the optimal rank for an input data set "sim2" with the default settings
#' opt.rank(sim2, nruns = 4)
#' # [1] 2
#' }
#'
#' @author Martin Ondrus, \email{mondrus@ualberta.ca}, Ivor Cribben, \email{cribben@ualberta.ca}
#' @references "Factorized Binary Search: a novel technique for change point detection in multivariate high-dimensional time series networks", Ondrus et al.
#' (2021), <arXiv:2103.06347>.

opt.rank = function(Y, nruns = 50, algtype = "brunet"){

  print("Finding optimal rank")
  Y = as.matrix(Y)

  # Create a permuted data set which will be compared with the original Y
  perm.subj = sample(as.vector(Y))
  perm.subj = matrix(perm.subj, ncol = ncol(Y))

  # Calculate the losses for original and permuted Y for first two rank values
  results.df = c()
  for (k in 1:2){
    # Fit NMF to the original and permuted Y
    orig.loss = nmf(Y, rank = k, nrun = nruns, method = algtype)@residuals
    perm.loss = nmf(perm.subj, rank = k, nrun = nruns, method = algtype)@residuals

    # Add these results to the results dataframe
    results.df = rbind(results.df, data.frame(k, orig.loss, perm.loss))
  }

  # Find the change in the original loss and permuted loss
  results.df[2,4] = results.df[2,2] - results.df[1,2]
  results.df[2,5] = results.df[2,3] - results.df[1,3]

  # Adjust the column names in the results dataframe
  colnames(results.df)[c(1,4,5)] = c("rank", "orig.change", "perm.change")

  # Loop which continues increasing rank while the decrease in loss for the original data is greater than the permuted
  k = 2
  while (results.df[k,4] < results.df[k,5]){
    # Add to the iterator
    k = k + 1

    # Fit NMF to the original and permuted data
    orig.loss = nmf(Y, rank = k, nrun = nruns, method = algtype)@residuals
    perm.loss = nmf(perm.subj, rank = k, nrun = nruns, method = algtype)@residuals

    # Find the change in loss for original and permuted Y
    orig.change = orig.loss - results.df[k-1,2]
    perm.change = perm.loss - results.df[k-1,3]

    # Add these results to the results dataframe
    results.df = rbind(results.df, c(k, orig.loss, perm.loss, orig.change, perm.change))
  }

  # Print the results and return the optimal rank
  print(paste("Optimal rank:", k))
  return(k)
}
