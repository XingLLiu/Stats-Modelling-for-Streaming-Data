change.pt <- function(response = NA, data = NA, phi.initial = NA,
                      niterations = 1000, precision = 0.95){
  # Input: 
  # Output:
  
  lambda <- 1
  delta <- 1
  n <- nrow(data)
  data <- cbind(data, matrix(rep(NA, n*2), ncol = 2))
  p <- ncol(data)
  phi <- phi.initial
  
  # Perform 1st iteration
  x2 <- data[1,2]
  condition <- as.numeric(x2 > phi)
  u <- (x2 - phi)*condition
  v <- -condition
  data[1,(p-1):p] <- c(u,v)
  # rls.fit <- rls(response[1], matrix(data[1,],nrow=1), lambda, delta)
  # beta.hat <- rls.fit$beta_hat.vec
  beta.hat <- t(solve(t(data[1,])%*%data[1,])%*%t(data[1,])*response[1])
  
  phi <- beta.hat[p]/as.numeric(beta.hat[p-1]) + phi
  ##########
  if (as.numeric(beta.hat[p-1]) == 0){
    phi <- phi.initial
  }
  #########
  
  for (i in 2:n){
    # Initial setting
    response.new <- response[i]
    data.new <- data[i,]
    
    changept.summary <- changept.iteration(response = response[1:(i-1)], 
                                           data = data[1:(i-1),],
                                           response.new = response.new, 
                                           data.new = data.new,
                                           phi = phi)
    
    data[1:i,] <- changept.summary$data
    ratio <- changept.summary$ratio
    phi <- changept.summary$phi
    
    # print((i == niterations) || (ratio < precision))
    # Check convergence
    # if ((i == niterations) || (ratio >= precision)){
    if (i == niterations){
      print('The value converged.')
      break
    }
  }
  
  return(phi)
  
}