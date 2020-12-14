#===========================================================================
# Estimates stationary networks using non-negative matrix factorization

#' estimate_network
#' @description Compliment function to FaBiSearch. Estimates sparse stationary networks using non-negative matrix factorization (NMF).
#'
#' @importFrom NMF nmf
#'
#' @param data Multivariate time series, \eqn{T}, to be analyzed, should be in a matrix format with time points in rows and variables in columns
#' @param alg.type Type of algorithm for NMF function -> check ?nmf for details, under "method"
#' @param n.runs Number of runs to use for NMF function
#' @param n.rank Rank value to be used in NMF calculations
#' @param method Either "cluster" to use hierarchical clustering on the consensus matrix, or "cutoff" to use a cutoff value for the consensus matrix
#' @param lambda Threshold parameter, relates to number of clusters if method = "cluster", or a value between 0 and 1 (inclusive) if method = "cutoff"
#'
#' @return Adjacency matrix noting temporal dependencies between variables from \eqn{T}
#' @export
#'
#' @examples

estimate_network = function(data, alg.type, n.rank, n.runs, method, lambda) {

  # data     = Multivariate time series, \eqn{T}, to be analyzed, should be in a matrix format with time points in rows and variables in columns
  # alg.type = Algorithm type -> check ?nmf for details, under "method"
  # n.rank   = Rank value to be used in NMF calculations
  # n.runs   = Number of runs to use for NMF function
  # method   = Either "cluster" to use hierarchical clustering on the consensus matrix, or "cutoff" to use a cutoff value for the consensus matrix
  # lambda   = Threshold parameter, relates to number of clusters if method = "cluster", or a value between 0 and 1 (inclusive) if method = "cutoff"

  # Run NMF on the data
  nmf.output = nmf(data, method = alg.type, rank = n.rank, nrun = n.runs)

  # Save the consensus matrix
  cons.matrix = nmf.output@consensus

  # Branch out different calculations based on method type
  if (method == "cluster"){
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

  } else if (method == "cutoff"){
    # Any values above lambda assigned 1, equal to or less than lambda assigned 0
    cons.matrix[cons.matrix > lambda] = 1
    cons.matrix[cons.matrix <= lambda] = 0

    # Save the final adjacency matrix
    adj.matrix = cons.matrix
    diag(adj.matrix) = 0
  }

  return(adj.matrix)
}
