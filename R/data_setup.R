#===========================================================================
# This function prepares the data for the other functions

#' data_setup
#' @description This function is embedded in the main FaBiSearch function. It serves to take the original input time series, \eqn{Y}, and conver it into a more convenient
#' list format. Note, this assumes that timecourse length for each subject is the same.
#'
#' @param data Multivariate time series, \eqn{Y}, to be analyzed, should be in a matrix format with time points in rows and variables in columns
#' @param n.subj Total number of subjects in the dataset
#' @param T Number of time points in the dataset for each subject
#'
#' @return List of original data, \eqn{Y}, where each element in the list is a subject and their associated multivariate time series
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
