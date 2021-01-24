#===========================================================================
# Estimates stationary networks using non-negative matrix factorization

#' est.net
#' @description Compliment function to FaBiSearch. Estimates sparse stationary networks using non-negative matrix factorization (NMF).
#'
#' @importFrom NMF nmf
#'
#' @param data Multivariate time series, \eqn{Y}, to be analyzed, should be in a matrix format with time points in rows and variables in columns
#' @param nruns Number of runs to use for NMF function, by default is set to =50
#' @param lambda Either a value >1 to use hierarchical clustering on the consensus matrix, or a between 0 and 1 inclusive to
#' use a cutoff value for the consensus matrix
#' @param rank Specifies whether to find rank using optimization procedure (specify "optimal") or use a predetermined rank (integer value, e.g., 4),
#' by default is set to ="optimal"
#' @param algtype Type of algorithm for NMF function, by default is set to ="brunet"
#'
#' @return Adjacency matrix noting temporal dependencies between variables from \eqn{Y}
#' @export

est.net = function(data, nruns = 50, lambda = 7, rank = "optimal", algtype = "brunet") {

  # Define the data as a matrix
  data = as.matrix(data)

  # If rank has not been specified, then it must be found
  if (rank == "optimal"){
    n.rank = fabisearch:::opt.rank(data, nruns, algtype)
  } else {
    n.rank = rank
    print(paste("User defined rank:", n.rank))
  }

  # Run NMF on the data
  nmf.output = nmf(data, method = algtype, rank = n.rank, nrun = nruns)

  # Save the consensus matrix
  cons.matrix = nmf.output@consensus

  # Branch out different calculations based on method type
  if (lambda > 1){
    # Run hclust and cut the tree at the prespecified n.rank
    hc.out = hclust(dist(cons.matrix))
    ct.out = cutree(hc.out, lambda)

    # Convert ct.out into a matrix factor
    matr.fact = matrix(nrow = lambda, ncol = nrow(cons.matrix))
    for (i in 1:lambda){
       matr.fact[i,] = ct.out == i
    }
    matr.fact = matr.fact*1

    # Construct the adjacency matrix by multiplying the matr.fact and it's
    adj.matrix = t(matr.fact) %*% matr.fact
    diag(adj.matrix) = 0

  } else if (lambda <= 1 && lambda >=0){
    # Any values above lambda assigned 1, equal to or less than lambda assigned 0
    cons.matrix[cons.matrix > lambda] = 1
    cons.matrix[cons.matrix <= lambda] = 0

    # Save the final adjacency matrix
    adj.matrix = cons.matrix
    diag(adj.matrix) = 0
  }

  return(adj.matrix)
}
