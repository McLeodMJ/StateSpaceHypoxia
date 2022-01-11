####### Integrates the log-likelihood over many data points#######

## From Biol 607 page : https://biol607.github.io/lab/07_likelihood.html#13_log-likelihood 
#' Second, as the density functions don’t take kindly to a vector of data and a vector of parameters,
#'  we’ll use rowwise() to iterate over rows, but ungroup() after for other operations. 
######################################################################################################
# y: data
# y_hat: lambda or mean of data structure
## if we want just the likelihood we would take prod instead of sum and make log = F


pois.likel <- function(y, lambda){
  # Data Generating Process
  y_hat <- lambda
  
  #likelihood function
  sum(dpois(y, lambda = y_hat, log = TRUE))
}