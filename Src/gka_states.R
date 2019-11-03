gka.states(response=NA, h.u=NA,h.l=NA,count.limit=3){
  # Input: response = column vector of observed values
  #        h.u = upper threshold
  #        h.l = lower threshold
  #        count.limit = max. positive integer over which H0 is rejected
  # Output: estimated start and end time of active states.
  # Assumption: only 2 states involved.
  
  phase.no <- 2
  epsilon <- 2/N
  time.hat <- data.frame(matrix(NA, ncol=2, nrow=1))
  colnames(time.hat) <- c('Start','End')
  gka.summary <- data.frame(matrix(NA, ncol=3, nrow=1))
  colnames(gka.summary) <- c('v','g','delta')
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
      prob.hat <- gka.p(y, gka.summary)
      p.val <- min(1 - prob.hat, prob.hat)
      p.mat[(n+1),] <- -2*log(p.val)
      
      # Identify if a change in state has occurred
      if (p.mat[(n+1),] > h.u && gka.switch == 0){
        gka.count[1] <- gka.count[1] + 1
        gka.count[-1] <- 0
        if (gka.count[1] == count.limit){
          start.time <- n+1-count.limit
          gka.switch <- 1
          gka.count[1] <- 0
          time.hat <- rbind(time.hat, data.frame('Start'=start.time,'End'=NA))
          sum.index <- 2
        }
      }
      else if (p.mat[(n+1),] < h.l && gka.switch == 1){
        gka.count[2] <- gka.count[2] + 1
        gka.count[-2] <- 0
        if (gka.count[2] == count.limit){
          end.time <- n+1-count.limit
          gka.switch <- 0
          gka.count[2] <- 0
          # print(paste('Ended:', end.time))
          time.hat[[nrow(time.hat)]][2] <- end.time
          sum.index <- 1
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
        gka.summary <- gka.insert(v, gka.summary)
        # Update size of summary
        s <- nrow(gka.summary)
      }
      
      # Update the no. of current data
      n.mat[sum.index] <- n.mat[sum.index] + 1
      # Update the current summary
      gka.summary.total[[sum.index]] <- gka.summary
    }
  }
  
  
  
}