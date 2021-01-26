#===========================================================================
# Helper function to find the optimal rank of a dataset

#' opt.rank
#' @description This function finds the optimal rank to use on a dataset using non-negative matrix factorization.
#'
#' @importFrom NMF nmf
#' @importFrom Rdpack reprompt
#'
#' @param Y A numerical matrix representing the multivariate time series, with the columns representing its components.
#' @param nruns A positive integer, by default is set to 50. Defines the number of runs to use for NMF function.
#' @param algtype A character string, which defines the algorithm to be used in the NMF function. By default is set to “brunet” - please see the "Algorithms"
#' section for more information on the available algorithms.
#'
#' @return An integer denoting the optimal rank found.
#' @export
#'
#' @section Algorithms:
#' All algorithms available are presented below, please note the "fabisearch" package builds upon the algorithms available in the "NMF" package
#' \insertCite{Gaujoux2010}{fabisearch}:\cr
#'
#' \code{"brunet"} - This algorithm is based on Kullback-Leibler divergence, from \insertCite{Brunet2004a}{fabisearch}. It uses multiplicative updates from
#' \insertCite{NIPS2000_1861}{fabisearch} with some small enhancements.\cr
#'
#' \code{"lee"} - This algorithm is based on Euclidian distances from \insertCite{NIPS2000_1861}{fabisearch}, and uses simple multiplicative updates.\cr
#'
#' \code{"ls-nmf"} - This is the least-squares NMF method from \insertCite{Wang2006}{fabisearch}. This algorithm uses an altered version of the Euclidian
#' distance based, multiplicative updates from \insertCite{NIPS2000_1861}{fabisearch}. It incorporates weights on each entry of the target matrix.\cr
#'
#' \code{"nsNMF"} - This is the nonsmooth NMF method from \insertCite{Pascual-Montano2006}{fabisearch}. This algorithm uses an altered version of the
#' Kullback-Leibler based, multiplicative updates from \insertCite{NIPS2000_1861}{fabisearch}. It includes an intermediate "smoothing" matrix, which is
#' intended to produce sparser factors.\cr
#'
#' \code{"offset"} - This is the offset NMF method from \insertCite{Badea2008}{fabisearch}. This algorithm uses an altered version of the Euclidian
#' distance based, multiplicative updates from \insertCite{NIPS2000_1861}{fabisearch}. It incorporates an intercept which is intended to reflect a common
#' pattern or baseline amongst components.\cr
#'
#' \code{"pe-nmf"} - This is the pattern-expression NMF from \insertCite{Zhang2008}{fabisearch}. This algorithm utilizes multiplicative updates to minimize
#' a Euclidian distance based objective function. It is further regularized such that the basis vectors effectively express patterns.\cr
#'
#' \code{"snmf/r","snmf/l"} - This is the alternating least-squares (ALS) approach from \insertCite{10.1093/bioinformatics/btm134}{fabisearch}. It uses the
#' non-negative, least-squares algorithm from \insertCite{VanBenthem2004}{fabisearch} to alternatingly estimate the basis and coefficent matrices. It
#' utilizes an Euclidian distance based objective function, and is regularized to promote either sparse basis ("snmf/l") or coefficent ("snmf/r") matrices \cr
#'
#' @examples
#' ## Finding the optimal rank for an input dataset "data" with the default settings
#' opt.rank(data)
#' [1] 5
#'
#' ## Finding the optimal rank for an input dataset "data" with nruns of 100 and the
#' ## default "brunet" algorithm
#' opt.rank(data, nruns = 100)
#' [1] 4
#'
#' ## Finding the optimal rank for an input dataset "data" using the least square
#' ## NMF method and the default nruns
#' opt.rank(data, algtype = "ls-nmf")
#' [1] 8
#'
#' @author Martin Ondrus, \email{mondrus@ualberta.ca}, Ivor Cribben, \email{cribben@ualberta.ca}
#' @references "Factorized Binary Search: a novel technique for change point detection in multivariate high-dimensional time series networks", Ondrus et al
#' (2021), preprint.
#'
#' \insertAllCited{}

opt.rank = function(Y, nruns = 50, algtype = "brunet"){

  print("Finding optimal rank")
  Y = as.matrix(Y)

  # Create a permuted dataset which will be compared with the original Y
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
