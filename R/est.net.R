#===========================================================================
# Estimates stationary networks using non-negative matrix factorization

#' Sparse network estimation using non-negative matrix factorization (NMF) for data between change points
#' @description This function estimates sparse networks using non-negative matrix factorization (NMF) for data between change points.
#'
#' @importFrom NMF nmf
#'
#' @param Y An input multivariate time series in matrix format, with variables organized in columns and time points in rows. All entries in Y must be positive.
#' @param lambda A positive real number, which defines the clustering method and/or the cutoff value when estimating an adjacency matrix from the computed
#' consensus matrix. If lambda = a positive integer value, say 6, complete-linkage, hierarchical clustering is applied to the consensus matrix and the cutoff is at
#' 6 clusters. If lambda is a vector of positive integer values, say c(4, 5, 6), the same clustering method is applied for each value sequentially. If lambda
#' = a positive real number, say 0.5, entries in the consensus matrix with a value greater than or equal to 0.5 are labeled 1, while entries less than 0.5 are
#' labeled 0. Similarly, if lambda is a vector of positive real numbers, say c(0.1, 0.3, 0.8), the same thresholding method is applied for each value
#' sequentially.
#' @param nruns A positive integer with default value equal to 50. It is used to define the number of runs in the NMF function.
#' @param rank A character string or a positive integer, which defines the rank used in the optimization procedure to detect the change points.
#' If rank = "optimal", which is also the default value, then the optimal rank is used. If rank = a positive integer value, say 4, then a predetermined
#' rank is used.
#' @param algtype A character string, which defines the algorithm to be used in the NMF function. By default it is set to "brunet". See the "Algorithms" section of
#' \code{\link[NMF]{nmf}} for more information on the available algorithms.
#' @param changepoints A vector of positive integers with default value equal to \code{NULL}. It is used to specify whether change points exist in the input \eqn{Y}, and
#' thus whether \eqn{Y} should be split into multiple stationary segments and networks estimated separately for each segment. If change points, say c(100, 200) are specified,
#' \eqn{Y} is split at the 100th and 200th row to correspond to 3 stationary segments. Each stationary segment is then estimated sequentially, and
#' a list is returned where each component corresponds to a stationary segment.
#'
#' @return A matrix (or more specifically, an adjacency matrix) denoting the network (or clustering) structure between components of \eqn{Y}. If lambda is a
#' vector, a list of adjacency matrices is returned, where each element of the list corresponds to an element in lambda.
#' @export
#'
#' @examples
#' \donttest{
#' ## Estimating the network for a multivariate data set, "sim2" with the settings:
#' ## nruns = 10 and lambda = 0.5 where the latter specifies the cutoff based method
#' est.net(sim2, lambda = 0.5, nruns = 4)
#' }
#'
#' @author Martin Ondrus, \email{mondrus@ualberta.ca}, Ivor Cribben, \email{cribben@ualberta.ca}
#' @references "Factorized Binary Search: a novel technique for change point detection in multivariate high-dimensional time series networks", Ondrus et al.
#' (2021), <arXiv:2103.06347>.

est.net = function(Y, lambda, nruns = 50, rank = "optimal", algtype = "brunet", changepoints = NULL) {

  # Define the Y as a matrix
  Y = as.matrix(Y)

  # If rank has not been specified, then it must be found
  if (rank == "optimal"){
    n.rank = opt.rank(Y, nruns, algtype)
  } else {
    n.rank = rank
    print(paste("User defined rank:", n.rank))
  }

  # Convert dataset to a list
  if(!is.null(changepoints)){
    changepoints = c(0, changepoints, nrow(Y))
    Y.split = list()
    for(i in 1:(length(changepoints) - 1)){
      # Define the time indices for the current block
      indices = (changepoints[i] + 1):changepoints[i + 1]
      # Save the current block in the list
      Y.split[[i]] = Y[indices, ]
    }

    # Save Y.split list as Y
    Y = Y.split
  } else if(is.null(changepoints)){
    Y = list(Y)
  }

  # Loop through Y and estimate networks
  output = list()
  for(j in 1:length(Y)){
    # Print progress
    print(paste("Estimating stationary block", j))

    # Run NMF on the Y
    nmf.output = nmf(Y[[j]], method = algtype, rank = n.rank, nrun = nruns)

    # Save the consensus matrix
    cons.matrix = nmf.output@consensus

    if (length(lambda) == 1){
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

        # Construct the adjacency matrix by multiplying the matr.fact and it's transpose
        adj.matrix = t(matr.fact) %*% matr.fact
        diag(adj.matrix) = 0

        # Add labels back to the matrix using the consensus matrix labels
        colnames(adj.matrix) = rownames(adj.matrix) = colnames(cons.matrix)

      } else if (lambda <= 1 && lambda >=0){
        # Any values above lambda assigned 1, equal to or less than lambda assigned 0
        cons.matrix[cons.matrix > lambda] = 1
        cons.matrix[cons.matrix <= lambda] = 0

        # Save the final adjacency matrix
        adj.matrix = cons.matrix
        diag(adj.matrix) = 0
      }
    } else if (length(lambda) > 1){
      # Branch out different calculations based on method type
      if (all(lambda > 1)){
        # Run hclust and cut the tree at the prespecified n.rank
        hc.out = hclust(dist(cons.matrix))

        # Loop through vector
        adj.matrix = list()
        for (i in 1:length(lambda)){
          ct.out = cutree(hc.out, lambda[i])

          # Convert ct.out into a matrix factor
          matr.fact = matrix(nrow = lambda[i], ncol = nrow(cons.matrix))
          for (j in 1:lambda[i]){
            matr.fact[j,] = ct.out == j
          }
          matr.fact = matr.fact*1

          # Construct the adjacency matrix by multiplying the matr.fact and it's transpose
          curr.adj.matrix = t(matr.fact) %*% matr.fact
          diag(curr.adj.matrix) = 0

          # Add labels back to the matrix using the consensus matrix labels
          colnames(curr.adj.matrix) = rownames(curr.adj.matrix) = colnames(cons.matrix)

          # Add to the adj.matrix object
          adj.matrix[[i]] = curr.adj.matrix
        }
      } else if (all(lambda <= 1 && lambda >= 0)){
        adj.matrix = list()
        for (i in 1:length(lambda)){
          # Any values above lambda assigned 1, equal to or less than lambda assigned 0
          curr.adj.matrix = cons.matrix
          curr.adj.matrix[curr.adj.matrix > lambda[i]] = 1
          curr.adj.matrix[curr.adj.matrix <= lambda[i]] = 0

          # Save the final adjacency matrix
          diag(curr.adj.matrix) = 0

          # Add to the adj.matrix object
          adj.matrix[[i]] = curr.adj.matrix
        }
      }
    }
    # Add the results to the output
    if (length(Y) == 1){
      output = adj.matrix
    } else if (length(Y) > 1){
      output[[j]] = adj.matrix
    }
  }

  return(output)
}
