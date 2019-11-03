changept.iteration <- function(response = NA, data = NA,
    response.new = NA, data.new = NA, phi = NA,
    niterations = 1000, precision = 0.95){
  # Input: data.new = a row vector of new data containing extra covariates,
  #        response.new = new response value containing extra covariates,
  #        phi.initial = initial guess of the change pt value,
  #        
  # Output:
  
  lambda <- 1
  delta <- 0.00001
  # Strictly speaking, p should be this minus 2
  p <- length(data.new)
  
  
  # Compute u and v
  x2 <- data.new[2]
  condition <- as.numeric(x2 > phi)
  u <- (x2 - phi)*condition
  v <- -condition
  
  # Fit the model with additional covariates
  data.new[(p-1):p] <- c(u,v)
  data <- matrix(rbind(data,data.new),ncol=p)
  response <- matrix(rbind(matrix(response,ncol=1),response.new),ncol=1)
  # rls.fit <- rls(response, data, lambda, delta)
  # beta.hat <- rls.fit$beta_hat.vec
  # beta.hat <- solve(t(data)%*%data)%*%t(data)%*%response
  N <- nrow(data)
  lm.fit <- lm(response[1:N]~data[1:N,]+0)
  beta.hat <- matrix(lm.fit$coefficients,ncol=1)
  print(beta.hat)
  
  # Improve the change-pt estimate 
  phi.new <- beta.hat[p]/as.numeric(beta.hat[p-1]) + phi
  #######
  if (as.numeric(beta.hat[p-1]) == 0){
    phi.new <- phi
  }
  #######
  
  # Check convergence
  ratio <- (phi.new - phi)/phi
  
  changept.summary <- list('response' = response,
                           'data' = data,
                           'phi' = phi.new,
                           'ratio' = ratio
                          )
  return(changept.summary)
}