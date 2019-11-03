set.seed(12)
N <- 100
M <- 2
beta_true.vector <- matrix(c(runif(n = M, min=-2, max=2)), ncol=1) 
sigma_squared <- 1
realizations <- matrix(runif(n = N, min = -10, max = 10), ncol = 1)
x.data <- rls.data_mat_setter(realizations, M, N)
response.data <- test.data_generator(beta_true.vector, x.data)
#
#response.data <- matrix(beta.true * realizations + rnorm(n = N, mean = 0, sd = 1), ncol = 1)
#

lambda <- 1
sigma <- 1
rls.fit <- rls(response.data, x.data , lambda, sigma)
