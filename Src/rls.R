rls<-function(response = NA, data = NA, lambda = NA, delta = NA, intercept = FALSE, 
              cp = FALSE, n.tes = 10, n.ref = 10, alpha = 2){
  # Input: a dataframe with the top row being the most recent data set
  #                    with one column of data and the most recent one at the end
  # Output: the least square vector of predictors
  
  
  # Set no. of predictors
  p <- ncol(data)
  if (length(p) == 0){
    p <- 1
  }
  data <- matrix(cbind(response, data), ncol = p+1)
  # Set no. of iterations
  n <- nrow(data)
  # Initialize a matrix for storing the LSE
  beta_hat.mat <- matrix(rep(NA, p*n), nrow = p)
  # Initialize a row vector for storing the RSS
  RSS.vec <- matrix(rep(NA, n), nrow = 1)
  # Initialize a row vector for storing R^2
  R_squared.vec <- matrix(rep(NA, n), nrow = 1)
  # Initialize a column vector for storing the estimated y(i) values
  y_hat.vec <- matrix(rep(NA, n), ncol = 1)
  
  # Initialize a matrix for storing the leverages
  # Column i of leverage.mat = leverages at time i
  leverage.mat <- matrix(rep(NA, n*n), ncol = n)
  # Initialize a matrix to locate any outliers
  outlier.mat <- matrix(rep(NA, n*n), ncol = n)
  # Initialize a matrix to store the studentized residuals
  stud_res.mat <- matrix(rep(NA, n*n), ncol = n)
  # Initialize a matrix to store the Cook's distance
  cook.mat <- matrix(rep(NA, n*n), ncol = n)
  
  # Initialize quantities for change-pt detection
  # Quandt's method
  count <- 0
  count.cts <- 0
  log.ratio <- matrix(rep(0,n),ncol = 1)
  ratio.bar <- matrix(rep(0,2*n),ncol=2)
  ratio.var <- matrix(rep(0,2*n),ncol=2)
  tao.hat <- numeric()
  # CUSUMSQ method
  cusumsq.vec <- matrix(rep(0,n),ncol = 1)
  cusumsq.max <- -Inf
  cusumsq.sum <- 0
  # Multiple models
  beta.vec.many <- numeric()
  
  
  ###############
  # Include the intercept term if specified
  if (intercept == TRUE){
    data['Intercept'] <- rep(1,n)
    #Find the intercept column
    l = which(colnames(data) == 'Intercept')
    #Rearrange data: response, intercept, ...
    data <- data[, c(1, l, c(1:ncol(data))[c(-1,-l)])]
    p <- ncol(data)-1
  }
  
  ###############
  #Initialize the 0th iteration
  beta_hat.vec <- as.matrix(rep(0,p), ncol=1)
  D <- delta^(-1)*diag(p)
  
  ###############
  # Initialize quantities for computing R^2
  RSS <- 0
  y.var <- 0
  y.bar <- 0
  ###############
  
  # Begin iterations
  for (i in (p+1):n ){

    ###############
    #Extract the response vector at time i
    y <- response[i]
    #Extract the data vector at time i
    x.vec <- data[i,-1]
    ###############
    # Detect hange-point by default
    # if (i >= (n.ref+n.tes) && cp == TRUE){
    #   cp.summary <- rls.cp(i, switch, count n.ref, n.tes, log.ratio,
    #                        response, data[,-1], lambda, delta)
    #   switch <- cp.summary$switch
    #   log.ratio[i] <- cp.summary$log.ratio[i]
    #   if (switch == 1){
    #     print(c('Change point detected: i = ', i))
    #     #Restart initialization
    #     beta_hat.vec <- as.matrix(rep(0,p), ncol=1)
    #     D <- delta^(-1)*diag(p)
    #     # # Start from the current iteration
    #     # results <- rls.iterations(x.vec, beta_hat.vec, D, y, lambda)
    #   }
    # }
    if (i >= (n.ref+n.tes) && cp == TRUE){
      cp.summary <- rls.cp(i, switch, count,count.cts, n.ref, n.tes, log.ratio, ratio.bar,
                            ratio.var, response, as.matrix(data[,-1]), lambda, 
                           delta,beta_hat.vec, D, tao.hat, cusumsq.vec, cusumsq.sum, alpha)
      # Set the next iteration
      count <- cp.summary$count
      count.cts <- cp.summary$count.cts
      log.ratio[i] <- cp.summary$log.ratio[i]
      ratio.bar[i,] <- cp.summary$ratio.bar[i,]
      ratio.var[i,] <- cp.summary$ratio.var[i,]
      beta_hat.vec <- cp.summary$beta_hat.vec
      D <- cp.summary$D
      tao.hat <- cp.summary$tao.hat
      
      beta_hat.vec.beforecp <- cp.summary$beta_hat.vec.beforecp
      if (length(beta_hat.vec.beforecp)>0){
        beta.vec.many <- cbind(beta.vec.many,beta_hat.vec.beforecp)
      }
    }
    ###############
    #Perform the current iteration
    results <- rls.iterations(x.vec, beta_hat.vec, D, y, lambda)
    ###############
    #Set the next iteration
    beta_hat.vec <- results$beta.hat
    D <- results$D
    ###############
    # Compute the leverage of the observations
    for (j in 1:i){
      leverage.mat[j,i] <- data[j,-1] %*% D %*% data[j,-1]
    }
    # Compute the studentized residuals
    e.vec <- (response[1:i] - data[1:i,-1]%*%beta_hat.vec)
    weighted.e.vec <- e.vec * sqrt(c(lambda^((i-1):0)))
    sigma.hat <- as.numeric(t(weighted.e.vec)%*%weighted.e.vec)/(n-p)
    leverage.vec <- leverage.mat[1:i,i]
    stud_res.vec <- e.vec/sqrt(abs(1-leverage.vec)*sigma.hat)
    stud_res.mat[1:i,i] <- stud_res.vec
    # Compute the Cook's distance
    cook.mat[1:i,i] <- stud_res.vec^2 * leverage.vec/((1-leverage.vec)*p)
    
    # # Only store the latest leverage
    # y.hat <- results$y_hat.val
    # leverage <- x.vec %*% D %*% x.vec
    # leverage.mat[i,N] <- leverage
    # e.vec <- response[1:i] - data[1:i,-1]%*%beta_hat.vec
    # sigma.hat <- as.numeric(t(e.vec)%*%(e.vec))/(i-p)
    # stand_res.vec <- e.vec/sqrt(1-leverage)
    # stud_res.mat[i,N] <- stand_res.vec/sqrt(sigma.hat)
    # cook.mat[i,N] <- stand_res.vec^2 * leverage/((1-leverage)*p)
    ###############
    # Compute R^2
    # y.bar <- y.bar + 1/i * (y - y.bar)
    # RSS <- lambda*RSS + 1/i * ((y - x.vec %*% beta_hat.vec)^2 - lambda*RSS)
    # y.var <- lambda*y.var + 1/i * ((y - y.bar)^2 - lambda*y.var)
    # R_squared <- (y.var - RSS)/y.var
    
    # Compute R^2
    y.bar <- y.bar + 1/i * (y - y.bar)
    zeta <- results$zeta
    RSS <- lambda*RSS + zeta*(y - x.vec %*% beta_hat.vec)
    if (i == 1){
      RSS <- (y - x.vec %*% beta_hat.vec)^2
    }
    y.var <- t(response[1:i] - y.bar)%*%(response[1:i] - y.bar)
    R_squared <- (y.var - RSS)/y.var
    ###############
    # Store statistic for the CUSUMSQ
    zeta <- results$zeta
    cusumsq.i <- (zeta/sqrt(1+leverage.vec[i]))^2
    cusumsq.vec[i] <- cusumsq.i
    cusumsq.sum <- cusumsq.sum + cusumsq.i
    ###############
    # Store the results into corresponding matrices
    beta_hat.mat[,i] <- beta_hat.vec
    RSS.vec[,i] <- RSS
    y_hat.vec[i,] <- results$y_hat.val
    R_squared.vec[,i] <- R_squared
  }
  
  ###############
  # Flag those with a high leverage
  # Remove the empty entries in outlier.mat
  outlier.mat <- which(leverage.mat > 2*p/n, arr.ind = TRUE)
  outlier.mat <- which(leverage.mat[,n] > 2*p/n, arr.ind = TRUE)
  
  # Return the results
    # Cp detection results
  cp.summary <- list('log.ratio' = log.ratio,
                     'ratio.bar' = ratio.bar,
                     'ratio.var' = ratio.var,
                     'tao.hat' = tao.hat
                     )
  gof.summary <- list('RSS.vec' = RSS.vec,
                      'R2.vec' = R_squared.vec,
                      'leverage.mat' = leverage.mat,
                       'outlier.mat' = outlier.mat,
                       'Studentized.residuals' = stud_res.mat,
                       'Cook.distance' = cook.mat
                     )
  
  #############
  beta.vec.many <- cbind(beta.vec.many,beta_hat.vec)
  #############
  rls.summary <-  list('beta_hat.vec' = as.matrix(beta_hat.mat[,n],ncol = 1), 
                       'beta_hat.mat' = beta_hat.mat,
                       'y_hat.vec' = y_hat.vec,
                       'RSS' = RSS.vec[,n], 
                       'RSS.vec' = RSS.vec,
                       'R2.vec' = R_squared.vec,
                       'leverage.vec' = leverage.mat[,n],
                       'leverage.mat' = leverage.mat,
                       'outlier.mat' = outlier.mat,
                       'Studentized.residuals' = stud_res.mat,
                       'Cook.distance' = cook.mat,
                       'gof.summary' = gof.summary,
                       'cp.summary' = cp.summary,
                       'beta.vec.many' = beta.vec.many
                       )

  
  return(rls.summary)
}
      
    
  