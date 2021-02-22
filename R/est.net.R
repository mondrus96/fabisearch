#===========================================================================
# Estimates stationary networks using non-negative matrix factorization

#' Sparse network estimation using non-negative matrix factorization (NMF) for data between change points
#' @description This function estimates sparse networks using non-negative matrix factorization (NMF) for data between change points.
#'
#' @importFrom NMF nmf
#'
#' @param Y A numerical matrix representing the multivariate time series, with the columns representing its components.
#' @param nruns A positive integer with default value equal to 50. It is used to define the number of runs in the NMF function.
#' @param lambda A positive real number, which defines the clustering method and/or the cutoff value when estimating an adjacency matrix from the computed
#' consensus matrix. If lambda = a positive integer value, say 6, complete-linkage, hierarchical clustering is applied to the consensus matrix and the cutoff is at
#' 6 clusters. If lambda = a positive real number, say 0.5, entries in the consensus matrix with a value greater than or equal to 0.5 are labelled 1,
#' while entries less than 0.5 are labelled 0.
#' @param rank A character string or a positive integer, which defines the rank used in the optimization procedure to detect the change points.
#' If rank = "optimal", which is also the default value, then the optimal rank is used. If rank = a positive integer value, say 4, then a predetermined
#' rank is used.
#' @param algtype A character string, which defines the algorithm to be used in the NMF function. By default it is set to "brunet". See the "Algorithms" section of
#' \code{\link[NMF]{nmf}} for more information on the available algorithms.
#'
#' @return A matrix (or more specifically, an adjacency matrix) denoting the network (or clustering) structure between components of \eqn{Y}.
#' @export
#'
#' @examples
#' ## Estimating the network for a multivariate data set, "sim2" using default settings
#' \donttest{est.net(sim2)}
#'
#' ## Estimating the network for a multivariate data set, "sim2", using hierarchical
#' ## clustering to generate the adjacency matrix with a cutoff value of 7 clusters
#' \donttest{est.net(sim2, nruns = 100, lambda = 7)}
#'
#' ## Estimating the network for a multivariate data set, "sim2", and using a cutoff
#' ## value for the adjacency matrix to enforce sparsity, where the cutoff is 0.5
#' \donttest{est.net(sim2, nruns = 100, lambda = 0.5)}
#'
#' ## Estimating the network for a multivariate data set, "sim2"
#' \donttest{est.net(sim2, rank = 4)}
#'
#' ## Estimating the network for a multivariate data set, "sim2", using the "snmf/l"
#' ## algorithm for NMF
#' \donttest{est.net(sim2, algtype = "snmf/l")}
#'
#' @author Martin Ondrus, \email{mondrus@ualberta.ca}, Ivor Cribben, \email{cribben@ualberta.ca}
#' @references "Factorized Binary Search: a novel technique for change point detection in multivariate high-dimensional time series networks", Ondrus et al.
#' (2021), preprint.

est.net = function(Y, nruns = 50, lambda = 7, rank = "optimal", algtype = "brunet") {

  # Define the Y as a matrix
  Y = as.matrix(Y)

  # If rank has not been specified, then it must be found
  if (rank == "optimal"){
    n.rank = opt.rank(Y, nruns, algtype)
  } else {
    n.rank = rank
    print(paste("User defined rank:", n.rank))
  }

  # Run NMF on the Y
  nmf.output = nmf(Y, method = algtype, rank = n.rank, nrun = nruns)

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
