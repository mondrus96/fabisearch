#===========================================================================
# Permutes input data by whole rows, that is, permute the time course of all ROIs together

#' permute_split
#' @description This function is embedded in the main FaBiSearch function. It serves to permute the input data by whole rows, that is, permute the time course of all ROIs
#' together
#'
#' @param data Multivariate time series, \eqn{Y}, to be analyzed, should be in a matrix format with time points in rows and variables in columns
#'
#' @return Permuted version of the original multivariate time series, \eqn{Y}
#' @export
#'
#' @examples

permute_split = function(data){

  # data = data to be read into the function

  # Save the row names/labels of the original dataset
  row.names = rownames(data)

  # Create random index
  random.index = sample(seq(1,dim(data)[1],1))

  # Randomize data and change row names to original
  data = data[random.index, ]
  rownames(data) = row.names

  return(data)
}
