########################
# Base and active states
# Initial setting
alpha <- 0.05
F <- 'qnorm'
h.u <- qchisq(0.99, 2)
h.l <- qchisq(0.95, 2)
learn.time <- 200
count.limit <- 4

# Generate data
N <- 800
N.test <- N - learn.time
response <- matrix(NA, nrow=N, ncol=1)
response <- rnorm(N, 0, 1)
response[1:learn.time] <- rnorm(learn.time, 0, 1)
response[(learn.time+1):(learn.time+N.test/3)] <- rnorm(N.test/3, 0, 1)
response[(learn.time+N.test/3+1):(learn.time+2*N.test/3)] <- rnorm(N.test/3, 5, 1)
response[(learn.time+2*N.test/3+1):(N)] <- rnorm(N.test/3, 0, 1)
# Generate data
N <- 1200
state.no <- 5
N.test <- N - learn.time
response <- matrix(NA, nrow=N, ncol=1)
response <- rnorm(N, 0, 1)
response[1:learn.time] <- rnorm(learn.time, 0, 1)
response[(learn.time+1):(learn.time+N.test/state.no)] <- rnorm(N.test/state.no, 0, 1)
response[(learn.time+N.test/state.no+1):(learn.time+2*N.test/state.no)] <- rnorm(N.test/state.no, 5, 1)
response[(learn.time+2*N.test/state.no+1):(learn.time+3*N.test/state.no)] <- rnorm(N.test/state.no, 0, 1)
response[(learn.time+3*N.test/state.no+1):(learn.time+4*N.test/state.no)] <- rnorm(N.test/state.no, 5, 1)
response[(learn.time+4*N.test/state.no+1):(N)] <- rnorm(N.test/state.no, 0, 1)


epsilon <- 2/N
gka.summary <- data.frame(matrix(NA, ncol=3, nrow=1))
colnames(gka.summary) <- c('v','g','delta')
s <- 0
p.mat <- matrix(0,ncol=1,nrow=N)
gka.switch <- 0
gka.count.base <- 0
gka.count.act <- 0

for (n in 0:(length(response)-1)){
  if (n >learn.time){
    # Compute p-values
    y <-response[n+1]
    # Find max. index j such that Y_j <= Y
    j <- which(gka.summary[,1] > y)[1] - 1
    if (is.na(j)){
      j <- s
    }
    # What if y is smaller than all entries in gka.summary?

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
      gka.count.base <- gka.count.base + 1
      gka.count.act <- 0
 
      if (gka.count.base == count.limit){
        start.time <- n+1-count.limit
        gka.switch <- 1
        gka.count.base <- 0
        print(paste('Started:', start.time))
      }
    }
    else if (p.mat[(n+1),] < h.l && gka.switch == 1){
      gka.count.act <- gka.count.act + 1
      gka.count.base <- 0
      if (gka.count.act == count.limit){
        end.time <- n+1-count.limit
        gka.switch <- 0
        gka.count.act <- 0
        print(paste('Ended:', end.time))
      }
    }
    else {
      gka.count.base <- 0
      gka.count.act <- 0
    }
  }
  
  if (gka.switch == 0){
    v <- response[n+1]
    # Fill in the first 1 iteration
    if (n == 0){
      gka.summary[1,] <- c(v,1,0)
      # Update size of summary
      s <- nrow(gka.summary)
    }
    else{
      
      if (n %% (1/(2*epsilon)) == 0){
        
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
      s <- s + 1
    }
  }
}
plot(p.mat)
abline(v=(learn.time+N.test/3), col='red', lty=2)
abline(v=(learn.time+2*N.test/3), col='red', lty=2)


#########################################
#########################################
# Performance analysis
# Time needed: 8mins
L <- 600
time.mat <- matrix(0,nrow=2, ncol=L)
for (l in 1:L){
  # Generate data
  learn.time <- 100
  
  # N <- 850
  # N.test <- N - learn.time
  # response <- matrix(NA, nrow=N, ncol=1)
  # response <- rnorm(N, 0, 1)
  # response[1:learn.time] <- rnorm(learn.time, 0, 1)
  # response[(learn.time+1):(learn.time+N.test/3)] <- rnorm(N.test/3, 0, 1)
  # response[(learn.time+N.test/3+1):(learn.time+2*N.test/3)] <- rnorm(N.test/3, 5, 1)
  # response[(learn.time+2*N.test/3+1):(N)] <- rnorm(N.test/3, 0, 1)
  # 
  N <- 850
  learn.time <- 100
  transit.time <- 100
  N.test <- N - learn.time
  response <- matrix(NA, nrow=N, ncol=1)
  response <- rnorm(N, 0, 1)
  response[1:learn.time] <- rnorm(learn.time, 0, 1)
  response[(learn.time+1):(learn.time+N.test/3)] <- rnorm(N.test/3, 0, 1)
  response[(learn.time+N.test/3+1):(learn.time-transit.time+2*N.test/3)] <- rnorm((N.test/3-transit.time), 10, 1)
  for (i in (learn.time-transit.time+2*N.test/3+1):(learn.time+2*N.test/3)){
    response[i] <- rnorm(1, (-0.1*i+60), 1)
  }
  response[(learn.time+2*N.test/3+1):(N)] <- rnorm(N.test/3, 0, 1)
  
  # Initial setting
  h.u <- qchisq(0.99, 2)
  h.l <- qchisq(0.95, 2)
  count.limit <- 3
  phase.no <- 2
  
  epsilon <- 2/N
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
          time.mat[1,l] <- start.time
          gka.switch <- 1
          gka.count[1] <- 0
          # print(paste('Started:', start.time))
          sum.index <- 2
        }
      }
      else if (p.mat[(n+1),] < h.l && gka.switch == 1){
        gka.count[2] <- gka.count[2] + 1
        gka.count[-2] <- 0
        if (gka.count[2] == count.limit){
          end.time <- n+1-count.limit
          time.mat[2,l] <- end.time
          gka.switch <- 0
          gka.count[2] <- 0
          # print(paste('Ended:', end.time))
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


###############
h2.duration<-hist(time.mat2[2,]-time.mat2[1,])
h2.duration$density <- h2.duration$counts/N
h1.duration<-hist(time.mat1[2,]-time.mat1[1,])
h1.duration$density <- h1.duration$counts/N
h3.duration<-hist(time.mat3[2,]-time.mat3[1,])
h3.duration$density <- h3.duration$counts/N
h2<-hist(time.mat2)
h1<-hist(time.mat1)
h3<-hist(time.mat3)

par(mfrow=c(2,3)) 
plot(h2.duration,freq=FALSE,
     main=expression(paste('Density of estimated duration of active state. ', Delta,tau,'=250')),
     xlab='Time length',
     ylab='Density',xlim=c(245,270)
     )

plot(h1.duration,freq=FALSE,
     main=expression(paste('Density of estimated duration of active state. ', Delta,tau,'=250')),
     xlab='Time length',
     ylab='Density'#,xlim=c(245,265)
     )

plot(h3.duration,freq=FALSE,
     main=expression(paste('Density of estimated duration of active state. ', Delta,tau,'=250.','Continuous change.')),
     xlab='Time length',
     ylab='Density'
     )

plot(h2,freq=FALSE,
     main=expression(paste('Estimated change points. ', mu[0],'=0',mu[1],'=5')),
     xlab='Time',
     ylab='Density')

plot(h1,freq=FALSE,
     main=expression(paste('Estimated change points. ', mu[0],'=0',mu[1],'=10')),
     xlab='Time',
     ylab='Density')

plot(h3,freq=FALSE,
     main=expression(paste('Estimated change points. ', mu[0],'=0',mu[1],'=10.', 'Continuous change.')),
     xlab='Time',
     ylab='Density')




