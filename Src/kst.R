kst <- function(gka.summary=NA, alpha=0.01, F=NA, ...){
  # Input: gka.summary = a summary of tuples returned by GKA
  #        alpha       = significance value of the test,
  #                      must be one of 0.001, 0.01, 0.05, 0.1
  #        F           = string of the name of the referemce distribution
  #        ...         = argument(s) of the function F
  # Output: KS test
  
  N <- sum(gka.summary[,2])
  s <- nrow(gka.summary)
  D.hat <- 0
  F <- get(F, mode = "function", envir = parent.frame())
  
  for (i in 1:s){
    y <- gka.summary[i,1]
    # Find max. index j such that Y_j <= Y
    j <- which(gka.summary[,1] > y)[1] - 1
    if (is.na(j)){
      j <- s
    }
    # Find the approximate index of Y_j in the stream
    i_j.max <- sum(gka.summary[1:j,2]) + gka.summary[j,3]
    i_j.hat <- gka.quantile(gka.summary, ranking=i_j.max)
    i_j.hat <- i_j.max
    
    y_j <- gka.summary[j,1]
    # if (F == 'norm'){
    #   E.hat <- abs(i_j.hat/N - pnorm(y_j, mean = 0, sd = 1))
    # }
    # else if (F == 'exp'){
    #   E.hat <- abs(i_j.hat/N - pexp(y_j, rate=1))
    # }
    E.hat <- abs(i_j.hat/N - F(y_j, ...))
    D.hat <- max(D.hat, E.hat)
  }
  
  if (alpha == 0.001){
    coeff <- 1.94947
  }
  else if (alpha == 0.01){
    coeff <- 1.62762
  }
  else if (alpha == 0.05){
    coeff <- 1.35810
  }
  else if (alpha == 0.1){
    coeff <- 1.22385
  }
  
  k <- coeff/sqrt(N) 
  # return(paste('Statistic: ', format(D.hat,digits=5), 
  #              ' ', format(alpha,digits=5), 
  #              '-level critical value: ', k))
  return(list('result' = paste('Statistic: ', format(D.hat,digits=5), 
                            ' ', format(alpha,digits=5), 
                            '-level critical value: ', k),
              'statistic' = D.hat
              ))
}