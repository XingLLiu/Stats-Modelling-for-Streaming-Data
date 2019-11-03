gka.p <- function(x=NA, gka.summary=NA){
  # Input: gka.summary = a data frame returned by the gka function
  #        x = a realization
  # Output: the empirical Pr(X <= x)
  
  N <- sum(gka.summary[,2])
  s <- nrow(gka.summary)
  # Find max. index j such that x_j <= x
  j <- which(gka.summary[,1] > x)[1] - 1
  if (is.na(j)){
    j <- s
  }
  
  # Find the approximate index of x_j in the stream
  if (j != 0){
    i_j.max <- sum(gka.summary[1:j,2]) + gka.summary[j,3]
    i_j.hat <- i_j.max
  }
  else {
    i_j.hat <- 0
  }
  
  # estimated probability
  prob.hat <- i_j.hat/N
  
  return(prob.hat)
  
}