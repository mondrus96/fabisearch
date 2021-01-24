#===========================================================================
# Estimates stationary networks using non-negative matrix factorization

#' est.net
#' @description Compliment function to FaBiSearch. Estimates sparse stationary networks using non-negative matrix factorization (NMF). Each run of NMF
#' produces an adjacency matrix, and by averaging across all adjacency matrices we generate a consensus matrix where each entry is the probability of two
#' variables or nodes being clustered together
#'
#' @importFrom NMF nmf
#'
#' @param data A matrix with time points in rows and variables in columns. This input is the multivariate time series, \eqn{Y}, to be analyzed.
#' @param nruns A positive integer, by default is set to 5. Defines the number of runs to use for NMF function.
#' @param lambda Either a positive integer to use hierarchical clustering on the consensus matrix, or a positive real number to use as a cutoff value for
#' the consensus matrix.
#' @param rank A positive integer denoting the value of rank to use in the algorithm. By default is set to "optimal" so the function finds the optimal rank
#' to use.
#' @param algtype A character string denoting the type of algorithm for NMF function, by default is set to "brunet".
#'
#' @return Adjacency matrix noting temporal dependencies between variables from \eqn{Y}
#' @export
#'
#' @examples
#' ## Estimating the network for a multivariate dataset, "data" using default settings - outputs as an adjacency matrix
#' est.net(data)
#'
#' ## Estimating the network for a multivariate dataset, "data", specifying the number of runs to 100
#' ## and using hierarchical clustering to generate the adjacency matrix with a cutoff value of 7 clusters
#' est.net(data, nruns = 100, lambda = 7)
#'
#' ## Estimating the network for a multivariate dataset, "data", specifying the number of runs to 100
#' ## and using a cutoff value for the adjacency matrix to enforce sparsity, where the cutoff is 0.5
#' est.net(data, nruns = 100, lambda = 0.5)
#'
#' ## Estimating the network for a multivariate dataset, "data", specifying the rank beforehand at 4
#' est.net(data, rank = 4)
#'
#' ## Estimating the network for a multivariate dataset, "data", using the least square NMF method
#' est.net(data, algtype = "ls-nmf")
#'
#' @author Martin Ondrus, \email{mondrus@ualberta.ca}, Ivor Cribben, \email{cribben@ualberta.ca}
#' @references "Factorized Binary Search: a novel technique for change point detection in multivariate high-dimensional time series networks", Ondrus et al
#' (2021), preprint.

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
