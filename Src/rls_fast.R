rls.fast<-function(response = NA, data = NA, lambda = NA, delta = NA, intercept = FALSE){
  # Input: a dataframe with the top row being the most recent data set
  #                    with one column of data and the most recent one at the end
  # Output: the least square vector of predictors
  
  
  # Set no. of predictors
  p <- ncol(data)
  if (length(p) == 0){
    p <- 1 
  }
  response <- matrix(response, ncol=1)
  data <- matrix(cbind(response, data), ncol = p+1)
  p <- ncol(data)-1
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
  for (i in 1:n ){

    ###############
    #Extract the response vector at time i
    y <- response[i]
    #Extract the data vector at time i
    x.vec <- data[i,-1]
    ###############
    #Perform the current iteration
    results <- rls.iterations(x.vec, beta_hat.vec, D, y, lambda)
    ###############
    #Set the next iteration
    beta_hat.vec <- results$beta.hat
    D <- results$D
    ###############
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
    # Store the results into corresponding matrices
    beta_hat.mat[,i] <- beta_hat.vec
    RSS.vec[,i] <- RSS
    y_hat.vec[i,] <- results$y_hat.val
    R_squared.vec[,i] <- R_squared
  }
  
  ###############
  # Return the results
  rls.summary <-  list('beta_hat.vec' = as.matrix(beta_hat.mat[,n],ncol = 1), 
                       'beta_hat.mat' = beta_hat.mat,
                       'RSS' = RSS.vec[,n], 
                       'RSS.vec' = RSS.vec,
                       'R2.vec' = R_squared.vec,
                       'y_hat.vec' = y_hat.vec
  )
  
  
  return(rls.summary)
}


