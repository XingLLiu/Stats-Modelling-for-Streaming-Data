rls.iterations<-function( x.vec = NA, beta_hat.vec = NA, D = NA, 
                         y.val = NA, lambda = NA){
  # Input: n = the nth iteration, x.vec = the nth data vector x.vec(n)
  #        beta_hat.vec = beta.hat.vec(n-1), y.val = data vector at time n, D = D(n-1)
  
  gamma.vec <- lambda^(-1) * D %*% x.vec/as.numeric((1 + lambda^(-1)*t(x.vec) %*% D %*% x.vec))
  zeta <- as.numeric(y.val - t(beta_hat.vec) %*% x.vec)
  beta_hat.vec <- beta_hat.vec + gamma.vec * zeta
  D <- lambda^(-1)*D - lambda^(-1) * gamma.vec %*% t(x.vec) %*% D
  
  # Compute estimated Y value
  y_hat.val <- as.numeric(t(beta_hat.vec) %*% x.vec)
  
  output <- list('beta.hat' = beta_hat.vec, 
                 'D' = D, 
                 'y_hat.val' = y_hat.val,
                 'zeta' = zeta
                 )
  return(output)
}
  

