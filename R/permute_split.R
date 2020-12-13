#===========================================================================
# Permutes input data by whole rows, that is, permute the time course of all ROIs together

#' permute_split
#'
#' @param data
#'
#' @return
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
