#===========================================================================
# Helper function to find the optimal rank of a dataset

#' opt.rank
#' @description This function finds the optimal rank to use on a dataset using non-negative matrix factorization
#'
#' @importFrom NMF nmf
#'
#' @param data Multivariate time series, \eqn{Y}, to be analyzed, should be in a matrix format with time points in rows and variables in columns
#' @param nruns Number of runs to use for NMF function, by default is set to =50
#' @param algtype Type of algorithm for NMF function, by default is set to ="brunet"
#'
#' @return Integer denoting the optimal rank found
#' @export

opt.rank = function(data, nruns = 50, algtype = "brunet"){

  print("Finding optimal rank")
  data = as.matrix(data)

  # Create a permuted dataset which will be compared with the original data
  perm.subj = sample(as.vector(data))
  perm.subj = matrix(perm.subj, ncol = ncol(data))

  # Calculate the losses for original and permuted data for first two rank values
  results.df = c()
  for (k in 1:2){
    # Fit NMF to the original and permuted data
    orig.loss = nmf(data, rank = k, nrun = nruns, method = algtype)@residuals
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
    orig.loss = nmf(data, rank = k, nrun = nruns, method = algtype)@residuals
    perm.loss = nmf(perm.subj, rank = k, nrun = nruns, method = algtype)@residuals

    # Find the change in loss for original and permuted data
    orig.change = orig.loss - results.df[k-1,2]
    perm.change = perm.loss - results.df[k-1,3]

    # Add these results to the results dataframe
    results.df = rbind(results.df, c(k, orig.loss, perm.loss, orig.change, perm.change))
  }

  # Print the results and return the optimal rank
  print(paste("Optimal rank:", k))
  return(k)
}
