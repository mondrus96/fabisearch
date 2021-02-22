#===========================================================================
# Builds the permutation distrubtion

#' @importFrom NMF nmf

perm_distr = function(orig.splits, curr.subj, T, x, n.rep, n.rank, alg.type) {

  # orig.splits = splitting times for subjects, contains the ranks for left and right used as well
  # curr.subj   = data for current subject
  # x           = time series of time
  # T           = number of time points in the data set
  # n.rep       = number of times to run the NMF algorithm for statistical inference
  # n.rank      = value of rank to use in the NMF function
  # alg.type    = algorithm type -> check ?nmf for details, under "method"

  # Output will be saved as distr.results
  perm.results = list()

  # Find the splits where the loss metric is reduced
  reduced.splits = orig.splits[orig.splits$chg.loss < 0, ]

  # Order the vector from smallest to largest, add 0 to beggining and T to end
  split.times = c(0,sort(reduced.splits$T.split),T)

  for (ik in 1:(length(split.times)-2)){

    # Print the current split being evaluated
    print(paste("Permuting split at", split.times[ik+1]))

    # Define the lower and upper parts of the split, as well as the ind.interval of time points
    lower1  = split.times[ik]
    T.split = split.times[ik+1]
    upper1  = split.times[ik+2]

    # Define the T.block and permute
    T.block = curr.subj[which(x<=upper1 & x>lower1),]

    # Loop through to find the sum of the left and right sides for each run
    curr.results = c()
    for(i in 1:n.rep){
      perm.block = permute_split(T.block)

      # Fit NMF to the left and right sides
      l.NMF = nmf(perm.block[rownames(perm.block) %in% which(x<=T.split & x>lower1),], rank=n.rank, method=alg.type)
      r.NMF = nmf(perm.block[rownames(perm.block) %in% which(x<=upper1 & x>T.split),], rank=n.rank, method=alg.type)

      curr.results = c(curr.results, sum(l.NMF@residuals + r.NMF@residuals), use.names = FALSE)
    }

    # Compile results into refit.results matrix
    perm.results[[ik]] = curr.results
  }
  return(perm.results)
}
