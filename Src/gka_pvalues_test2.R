alpha <- 0.05
F <- 'qnorm'
h.u <- qchisq(0.99, 2)
h.l <- qchisq(0.95, 2)
learn.time <- 200
count.limit <- 5
phase.no <- 2
new.phase <- 0

epsilon <- 2/N
gka.summary <- data.frame(matrix(NA, ncol=3, nrow=1))
colnames(gka.summary) <- c('v','g','delta')
gka.summary.temp <- data.frame(matrix(NA, ncol=3, nrow=1))
colnames(gka.summary.temp) <- c('v','g','delta')
gka.summary.total <- list()
for (i in 1:phase.no){
  gka.summary.total[[i]] <- gka.summary
}
s <- 0
n.mat <- matrix(0,ncol=phase.no,nrow=1)
p.mat <- matrix(0,ncol=1,nrow=N)
gka.switch <- 0
gka.count <- matrix(0,ncol=phase.no,nrow=1)
sum.index <- 1

for (n in 0:(length(response)-1)){
  if (n >learn.time){
    # Compute p-values
    y <-response[n+1]
    gka.summary <- gka.summary.total[[1]]
    s <- nrow(gka.summary)
    # Find max. index j such that Y_j <= Y
    j <- which(gka.summary[,1] > y)[1] - 1
    if (is.na(j)){
      j <- s
    }
    
    # Find the approximate index of Y_j in the stream
    if (j != 0){
      i_j.max <- sum(gka.summary[1:j,2]) + gka.summary[j,3]
      i_j.hat <- i_j.max
    }
    else {
      i_j.hat <- 0
    }
    
    # p-value
    # n.base = sample size based on the baseline distribution
    n.base <- sum(gka.summary[,2])
    p.val <- min(1 - i_j.hat/n.base, i_j.hat/n.base)
    p.mat[(n+1),] <- -2*log(p.val)#p.val
    
    # Identify if a change in state has occurred
    if (p.mat[(n+1),] > h.u && gka.switch == 0){
      gka.count[1] <- gka.count[1] + 1
      gka.count[-1] <- 0
      if (gka.count[1] == count.limit){
        start.time <- n+1-count.limit
        gka.switch <- 1
        gka.count[1] <- 0
        print(paste('Started:', start.time))
        sum.index <- 2
        new.phase <- 1
      }
    }
    else if (p.mat[(n+1),] < h.l && gka.switch == 1){
      gka.count[2] <- gka.count[2] + 1
      gka.count[-2] <- 0
      if (gka.count[2] == count.limit){
        end.time <- n+1-count.limit
        gka.switch <- 0
        gka.count[2] <- 0
        print(paste('Ended:', end.time))
        sum.index <- 1
        new.phase <- 0
        # Reset contemporary summary 
        gka.summary.temp[,] <- NA
      }
    }
    else {
      gka.count[,] <- 0
    }
  }
  
  # if (gka.switch == 0){
  if (TRUE){
    v <- response[n+1]
    gka.summary <- gka.summary.total[[sum.index]]
    # Use the temporary summary if a new phase is detected
    if (new.phase == 1){
      gka.summary <- gka.summary.temp
    }
    # Fill in the first 1 iteration
    if (n.mat[sum.index] == 0){
      gka.summary[1,] <- c(v,1,0)
      # Update size of summary
      s <- nrow(gka.summary)
    }
    else{
      s <- nrow(gka.summary)
      if (n.mat[sum.index] %% (1/(2*epsilon)) == 0){
        
        i = s - 1
        while (i >= 2){
          
          j = i-1
          delta.i <- gka.summary[i,3]
          g.sum <- sum(gka.summary[j:i,2])
          v <- gka.summary[i,1]
          
          while (j >= 2 && ((g.sum + delta.i) < 2*epsilon*n)){
            j <- j - 1
            g.sum <- g.sum + gka.summary[j,2]
          }
          
          # Tune one index up
          j <- j + 1
          
          # DELETE phase
          if (j < i){
            # Merge tuples from j to i
            gka.summary <- gka.summary[-((j+1):i),]
            gka.summary[j,] <- data.frame('v'=v, 'g'=g.sum-gka.summary[(j-1),2], 'delta'=delta.i)
          }
          
          # Continue from the largest integer smaller than j
          i <- j - 1
          # Update size of the summary
          s <- nrow(gka.summary)
        }
      }
      
      
      # INSERT phase
      s <- nrow(gka.summary)
      v.0 <- gka.summary[1,1]
      v.s_1 <- gka.summary[s,1]
      
      # Extreme cases
      tuple.new <- data.frame('v'=NA, 'g'=NA, 'delta'=NA)
      if ( v < v.0 ){
        delta <- 0
        new.position <- 0
        gka.summary <- rbind(tuple.new, gka.summary)
      }
      else if ( v > v.s_1 ){
        delta <- 0
        new.position <- s
        gka.summary <- rbind(gka.summary, tuple.new)
      }
      else{
        # Find appropriate index i
        new.position <- which( v < gka.summary[,1] )[1] - 1
        delta <- gka.summary[new.position,2] + gka.summary[new.position,3] - 1
        gka.summary <- rbind(gka.summary, tuple.new)
        gka.summary[(new.position+2):(s+1), ] <- gka.summary[(new.position+1):s, ]
      }
      
      # Insert new tuple
      tuple.new <- data.frame('v'=v, 'g'=1, 'delta'=delta)
      gka.summary[(new.position+1),] <- tuple.new
      # Update size of summary
      s <- nrow(gka.summary)
    }

    # Update the no. of current data
    n.mat[sum.index] <- n.mat[sum.index] + 1
    # Update the current summary
    gka.summary.total[[sum.index]] <- gka.summary
  }
}



