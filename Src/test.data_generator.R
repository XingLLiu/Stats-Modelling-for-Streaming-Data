rls.test <- function(M = NA, N = NA, lambda = NA, delta = NA){
  # Input: N and M as in the RLS algorithm
  # Output: A list consisting of the true vector beta and beta^hat
  
  # Generate the true beta vector randomly
  beta_true.vector <- matrix(c(runif(n = M, min=-2, max=2)), ncol=1) 
  sigma_squared <- 1
  # Generate the data matrix X in prewindowing form
  realizations <- matrix(runif(n = N, min = -10, max = 10), ncol = 1)
  x.data <- rls.data_mat_setter(realizations, M, N)
  # Calculate Y
  response.data <- matrix(x.mat, nrow = N) %*% matrix(beta_true.vec, ncol = 1)
  #
  # response.data <- test.data_generator(beta_true.vector, x.data)
  # response.data <- matrix(beta.true * realizations + rnorm(n = N, mean = 0, sd = 1), ncol = 1)
  #
  
  # Calculate beta^hat
  rls.fit <- rls(response.data, x.data , lambda, delta)

  beta.comp <- list('beta_true' = beta_true.vector, 'beta_hat' = rls.fit)
  
  return(beta.comp)
  
}