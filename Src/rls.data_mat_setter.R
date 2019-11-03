rls.data_mat_setter <- function(data = NA, M = NA, N = NA){
  # Input: data = column vector consisting of the observed data
  # Output: An NxM prewindowing data matrix, X
  
  data <- matrix(data, ncol=1)
  data.mat <- matrix(rep(0,length(data)), ncol = N, nrow = M)
  
  #Initialize data vector at time i
  data.vec <- matrix(rep(0, M),ncol = 1)
  for (i in 1:N){
    
    if(i <= M){
      data.vec[1:i] <- data[i:1]
    } else {
      data.vec <- data[i:(i-M+1)]
    }
    
    #Store this vector into the prewindowing data matrix
    data.mat[,i] <- data.vec
    
  }
  
  
  return(t(data.mat))
}

