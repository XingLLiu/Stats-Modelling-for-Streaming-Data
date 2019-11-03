rls.cp <- function(i=NA, switch = NA, count=NA, count.cts=NA, n.ref = NA, n.tes = NA, log.ratio = NA, ratio.bar = NA,
                   ratio.var = NA, y=NA, data=NA, lambda=NA, 
                   delta=NA,beta_hat.vec=NA, D=NA, tao.hat=NA, cusumsq.vec=NA, cusumsq.sum=NA, alpha=NA){
  # Input: h = open reference point (h>0)
  # Output:
  
  # switch = 2 if a cp was detected previously, 1 if a cp was detected at this iteration, 0 if not.
  switch.new <- 0
  p <- ncol(data)
  # Let p=1 if the input data is a vector of ones
  if (length(p) == 0){
    p <- 1
  }
  beta_hat.vec.beforecp <- numeric()
  data <- matrix(data, ncol=p)
  
  rls.fit.ref <- rls.fast(y[(i-n.tes-n.ref+1):(i-n.tes)], data[(i-n.tes-n.ref+1):(i-n.tes),], lambda, delta)
  rss.ref <- rls.fit.ref$RSS

  # beta.hat.ref <- rls.fit.ref$beta_hat.vec
  # res.ref <- y[(t-n.tes+1):t]-data[(t-n.tes+1):t,] %*% beta.hat.ref
  # rss.ref <- t(res.ref)%*%res.ref
  
  rls.fit.tes <- rls.fast(y[(i-n.tes+1):i], data[(i-n.tes+1):i,], lambda, delta)
  rss.tes <- rls.fit.tes$RSS
  
  sigma.sq.hat <- as.numeric(rss.ref+rss.tes)/(n.tes+n.ref-p)
  log.ratio[i] <- ((n.ref-n.tes)/2)*log(2*pi) + ((n.ref-n.tes)/2)*log(sigma.sq.hat) + (1/2)*rss.tes + (1/2)*rss.ref
  
  ###########
  # Determine the threshold using MA
  sample.size <- min(i-n.ref-n.tes+1,n.ref+n.tes )
  # Ref MA
  ratio.bar[i,1] <- mean(log.ratio[(i-n.tes-n.ref+1):(i-n.tes)])
  ratio.var[i,1] <- ( 1/max((sample.size-1),1)*(sum(log.ratio[(i-n.tes-n.ref+1):(i-n.tes)]^2) + 
                                                  sample.size*ratio.bar[i,1]^2 - 
                                                  2*ratio.bar[i,1]*sum(log.ratio[(i-n.tes-n.ref+1):(i-n.tes)])) )
  # Tes MA
  ratio.bar[i,2] <- mean(log.ratio[(i-n.tes+1):i]) 
  # + sample.size*ratio.bar[i,2]^2 - 2*ratio.bar[i,2]*sum(log.ratio[(i-n.tes+1):i]) 
  
  # Compute the threshold value
  ratio.ref.ub <- ratio.bar[i,1]+sqrt(ratio.var[i,1])
  if (i >= (n.ref*3+n.tes) && ratio.bar[i,2] > ratio.ref.ub){
    # Start counting
    count <- count + 1
  }
  else {
    count <- 0
  }
  
  #######
  # Raise an alarm if a change point is confirmed. 
  if (count == (n.tes )){
    
    # Double check with the CUSUMSQ test
    cusumsq.results <- rls.cusumsq(p, i, cusumsq.vec, cusumsq.sum, alpha)
    cusumsq.switch <- as.numeric(cusumsq.results$switch)
    cusumsq.stat <- cusumsq.results$cusumsq.stat
    critical.value <- cusumsq.results$critical.value
    significance.level <- cusumsq.results$significance.level
    
    # Confirm the change point if it passes the two tests
    if (cusumsq.switch == 1){
      cp <- (i-n.tes+1)
      print(paste('Change point detected: i = ', cp))
      print(paste('CUSUMSQ statistic: ', format(cusumsq.stat, digits = 6),
                '(Critical value: )', format(critical.value, digits = 6),
                'Rejected at significance level: ', significance.level))
      
      tao.hat <- cbind(tao.hat,cp)
      # Reset count for the next cp
      count <- 0
      # Reset initialization
      beta_hat.vec.beforecp <- beta_hat.vec
      beta_hat.vec <- as.matrix(rep(0,p), ncol=1)
      D <- delta^(-1)*diag(p)
    }
  }
  
  ########
  # Detect continuous cp
  if (i >= (n.ref*2+n.tes) && ratio.bar[i,2] >= ratio.bar[(i-1),2]){
    count.cts <- count.cts + 1
  }
  else{
    count.cts <- 0
  }

  if (count.cts == n.tes){
    # print(c('Continuous change point detected: i = ', i))
    count.cts <- 0
  }
  
  cp.summary <- list('switch' = switch.new,
                     'count' = count,
                     'count.cts' = count.cts,
                     'log.ratio' = log.ratio,
                     'ratio.bar' = ratio.bar,
                     'ratio.var' = ratio.var,
                     'ratio.ref.ub'= ratio.ref.ub,
                     'beta_hat.vec' = beta_hat.vec,
                     'D' = D,
                     'tao.hat' = tao.hat,
                     'beta_hat.vec.beforecp' = beta_hat.vec.beforecp
                     )
  
  return(cp.summary)
  
}