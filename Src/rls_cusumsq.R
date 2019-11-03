rls.cusumsq <- function(p=NA, i=NA, cusumsq.vec=NA, cusumsq.sum=NA, alpha=NA){
  # Input: 
  # Output: The CUSUMSQ statistic and decision rule.
  
  # Write up for multiple p-values
  if (alpha == 1){
    a1 <- 1.0729830
    a2 <- -0.6698868
    a3 <- -0.5816458
    significance.level <- '10%'
  }
  else if (alpha == 2){
    a1 <- 1.2238734
    a2 <- -0.6700069
    a3 <- -0.7351697
  significance.level <- '5%'
  }
  else if (alpha == 3){
    a1 <- 1.3581015
    a2 <- -0.6701218
    a3 <- -0.8858694
    significance.level <- '2.5%'
  }
  else if (alpha == 4){
    a1 <- 1.5174271
    a2 <- -0.6702672
    a3 <- -1.0847745
    significance.level <- '1%'
  }
  else {
    a1 <- 1.6276236
    a2 <- -0.6703724
    a3 <- -1.2365861
    significance.level <- '0.5%'
  }
  
  # 0 = H_0 not rejected, 1 = H_0 rejected
  switch <- 0
  
  # Compute the critical value
  c <- a1/sqrt(i) + a2/i + a3/i^{3/2}
  
  mean <- (c((p+1):i)-p)/(i-p)
  cusumsq.stat <- max(abs(cumsum(cusumsq.vec[(p+1):i])/cusumsq.sum - mean))
  
  if (cusumsq.stat > c){
    switch <- 1
  }
  
  summary <- list( 'switch' = switch,
                   'cusumsq.stat' = cusumsq.stat,
                   'critical.value' = c,
                   'significance.level' = significance.level 
                  )
 }