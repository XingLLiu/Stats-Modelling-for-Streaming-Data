gka.quantile <- function(gka.summary = NA, phi = NA, ranking = NA){
  # Input: gka.summary = a data frame returned by the gka function
  #        phi = the percentage quantile of interest
  #              If phi*n is not a integer, the floor function is used,
  #        rank = the rank of interest (alternative to phi)
  # Output: response = the observation corresponding to the specified phi or rank
  
  s <- nrow(gka.summary)
  n <- sum(gka.summary[,2])
  r.max.vec <- matrix(NA,nrow=s, ncol=1)
  
  if (!is.na(phi)){
    r <- floor(n*phi)
  }
  else {
    r <- ranking
  }
  
  # Compute the error
  e <- max(gka.summary[,2] + gka.summary[,3])/2
  
  # Choose the corresponding response from the GKA tuples
  if (r > (n - e)){
    response <- gka.summary[s,1]
  }
  else {
    r.max.vec <- cumsum(gka.summary[,2]) + gka.summary[,3]
    response <- gka.summary[which(r.max.vec>(r + e))[1],1]
  }
  
  return(response)
}