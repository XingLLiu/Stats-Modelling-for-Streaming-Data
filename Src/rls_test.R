rls.test <- function(M = NA, N = NA, lambda = NA, delta = NA){
  # Input: N and M as in the RLS algorithm
  # Output: A list consisting of the true vector beta and beta^hat
  
  # Generate the true beta vector randomly
  beta_true.vec <- matrix(c(runif(n = M, min=-2, max=2)), ncol=1) 
  sigma_squared <- 1
  # Generate the data matrix X in prewindowing form
  realizations <- matrix(runif(n = N, min = -10, max = 10), ncol = 1)
  x.mat <- rls.data_mat_setter(realizations, M, N)
  # Calculate Y
  response.data <- matrix(x.mat, nrow = N) %*% matrix(beta_true.vec, ncol = 1)
  #
  # response.data <- test.data_generator(beta_true.vec, x.mat)
  # response.data <- matrix(beta.true * realizations + rnorm(n = N, mean = 0, sd = 1), ncol = 1)
  #
  
  # Calculate beta^hat
  # Start timing
  ptm <- proc.time()
  rls.fit <- rls(response.data, x.mat , lambda, delta)
  # End timing
  time <- proc.time() - ptm
  rls.summary <- list('beta_true' = beta_true.vec, 
                    'beta_hat' = rls.fit$beta_hat.vec, 
                    'y_hat.vec' = rls.fit$y_hat.vec,
                    'RSS' = rls.fit$RSS,
                    'R2.vec' = rls.fit$R2.vec,
                    'y_hat.vec' = rls.fit$y_hat.vec,
                    'y.vec' = response.data,
                    'time_elapsed' = time
                    )
  
  return(rls.summary)
  
}