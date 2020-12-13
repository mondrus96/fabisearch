#===========================================================================
# This function prepares the data for the other functions

#' data_setup
#'
#' @param data
#' @param n.subj
#'
#' @return
#' @export
#'
#' @examples

data_setup = function(data,n.subj,T){

  # data   = data set
  # n.subj = number of subjects in the data set
  # T      = number of time points in the data set

  # Defining the actual output
  output = list()
  for (i in 1:n.subj){
    start = ((i-1)*T) + 1
    end = i*T
    curr.subj = data[start:end, ]
    rownames(curr.subj) = 1:T
    output[[i]] = curr.subj
  }
  return(output)
}
