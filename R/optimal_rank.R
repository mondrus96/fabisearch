#===========================================================================
# This function finds the optimal rank for the given subject using the change in loss method

#' optimal_rank
#' @description This function is embedded in the main FaBiSearch function. It serves to find the optimal rank for the given subject using the change in loss method.
#'
#' @importFrom NMF nmf
#'
#' @param curr.subj Multivariate time series data from the current subject
#' @param n.runs Number of runs to use for NMF function
#' @param alg.type Type of algorithm for NMF function -> check ?nmf for details, under "method"
#'
#' @return Optimal rank for the subject, an integer
#' @export
#'
#' @examples

optimal_rank = function(curr.subj, n.runs, alg.type){

  # curr.subj   = multivariate time series data from the current subject
  # n.runs      = n.runs to try the NMF algorithm for, the higher the more exhaustive the search for an optimal matrix is
  # alg.type    = algorithm type -> check ?nmf for details, under "method

  print("Finding optimal rank")

  # Create a permuted dataset which will be compared with the original data
  perm.subj = sample(as.vector(curr.subj))
  perm.subj = matrix(perm.subj, ncol = ncol(curr.subj))

  # Calculate the losses for original and permuted data for first two rank values
  results.df = c()
  for (k in 1:2){
    # Fit NMF to the original and permuted data
    orig.loss = nmf(curr.subj, rank = k, nrun = n.runs, method = alg.type)@residuals
    perm.loss = nmf(perm.subj, rank = k, nrun = n.runs, method = alg.type)@residuals

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
    orig.loss = nmf(curr.subj, rank = k, nrun = n.runs, method = alg.type)@residuals
    perm.loss = nmf(perm.subj, rank = k, nrun = n.runs, method = alg.type)@residuals

    # Find the change in loss for original and permuted data
    orig.change = orig.loss - results.df[k-1,2]
    perm.change = perm.loss - results.df[k-1,3]

    # Add these results to the results dataframe
    results.df = rbind(results.df, c(k, orig.loss, perm.loss, orig.change, perm.change))
  }

  # Print the results and return the optimal rank
  print(results.df)
  print(paste("Optimal rank:", k))
  return(k)
}
