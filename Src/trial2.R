# CUSUM Robustness test
# Time needed: 40 mins
start.time <- Sys.time()
tao.hat <- matrix(rep(0, 600*4), nrow=4)
for (k in 1:600){
 
  N <- 301
  p <- 3
  intercept <- matrix(rep(1,N), ncol = 1)
  data <- cbind(intercept, matrix(c(runif(n=N*(p-1), min=-10,max=10)),ncol=(p-1)))
  beta <- matrix(rep(NA, N*(p-1)),ncol = (p-1))
  beta <- t(cbind(rep(5,N),beta))
  y <- matrix(rep(NA,N),ncol = 1)
  # #########
  # beta[2,]<- -8
  # beta[3,]<-3
  beta[2,]<- -8
  beta[3,1:149]<-1.5
  beta[3,150:N]<- -0.3
  
  if (k > 300){
    beta[3,] <- 1.5
  }
  #########
  # beta[2,]<- -8
  # beta[3,1:399]<-1.5
  # beta[3,400:699]<-0.4
  # beta[3,700:1000]<--0.3
  #########
  for (i in 1:N){
    y[i] <- data[i,]%*%beta[,i]
  }
  y <- y + matrix(c(rnorm(n=N, mean=0, sd=1)), ncol=1)
  
  ###########
  # CUSUM
  lambda<- 1
  delta <-  0.00001
  res <- matrix(rep(0, N*N), ncol=N)
  cusum.stat <- matrix(rep(0, N), ncol=1)
  cusum.count <- 0
  a<- 0.9479
  ref.val <- matrix(rep(0,N),ncol=1)
  ref.val <- 2*a/sqrt(N-p) * (1:N) + (a*N-3*a*p)/sqrt(N-p)
  ###########
  # Quandts
  n.ref <- 20
  n.tes <- 20
  log.ratio <- matrix(rep(0,N),ncol=1)
  ratio.bar <- matrix(rep(0,2*N),ncol=2)
  ratio.var <- matrix(rep(0,2*N),ncol=2)
  quandts.count <- 0
  quandts.stop <- 0
  ###########
  # Quandts improved
  quandts.imp.stud.res <- matrix(rep(0, N),ncol=1)
  quandts.imp.count <- 0
  ###########
  # CUSUMSQ
  rls.fit.sq <- rls.fast(response = y, data = data, lambda = 1,delta = 0.00001) 
  beta.hat.mat <- rls.fit.sq$beta_hat.mat
  stud.res <- matrix(rep(0, N),ncol=1)
  ###########
  # KS Test
  

  # start.time <- Sys.time()
  for (t in (p+1):N){
      #############
      #CUSUM
      rls.fit <- rls.fast(y[1:t], data[1:t,], lambda, delta)
      # res[1:t,t] <- matrix(diag(rls.fit$Studentized.residuals), ncol=1)
      # cusum.stat[t] <- sum(res[(p+1):t,t])
      cusum.res <- y[t] - (rls.fit$y_hat.vec)[t]
      # cusum.res <- y[t] - data[t,]%*%rls.fit$beta_hat.vec
      cusum.lev <- data[t,]%*%solve(t(data[1:t,]) %*% data[1:t,])%*%data[t,]
      res[1,t] <- cusum.res/sqrt(1-cusum.lev)
      cusum.sigma.sq.hat <- (rls.fit$RSS)/(t-p)
      # cusum.stat[t] <- sum(res[1,])/sqrt(cusum.sigma.sq.hat)
      cusum.stat[t] <- sum(res[1,])/sqrt(cusum.sigma.sq.hat)
     
      if (abs(cusum.stat[t]) > ref.val[t] && cusum.count == 0){
        tao.hat[1, k] <- t
        cusum.count <- 1
      }

      # ###########
      if (t > (n.ref+n.tes)){
        # Quandts
        rls.fit.ref <- rls.fast(y[(t-n.tes-n.ref+1):(t-n.tes)], data[(t-n.tes-n.ref+1):(t-n.tes),], lambda, delta)
        beta.hat.ref <- rls.fit.ref$beta_hat.vec
        res.ref <- y[(t-n.tes+1):t]-data[(t-n.tes+1):t,] %*% beta.hat.ref
        rss.ref <- t(res.ref)%*%res.ref
  
        rls.fit.tes <- rls.fast(y[(t-n.tes+1):t], data[(t-n.tes+1):t,], lambda, delta)
        rss.tes <- rls.fit.tes$RSS
  
        sigma.sq.hat <- as.numeric(rss.ref+rss.tes)/(n.tes+n.ref-p)
        log.ratio[t] <- ((n.ref-n.tes)/2)*log(2*pi) + ((n.ref-n.tes)/2)*log(sigma.sq.hat) + (1/2)*rss.tes + (1/2)*rss.ref
  
        sample.size <- min(t-n.ref-n.tes+1,n.ref+n.tes )
        # Ref
        ratio.bar[t,1] <- mean(log.ratio[(t-n.tes-n.ref+1):(t-n.tes)])
        ratio.var[t,1] <- ( 1/max((sample.size-1),1)*(sum(log.ratio[(t-n.tes-n.ref+1):(t-n.tes)]^2)
                                                      + sample.size*ratio.bar[t,1]^2
                                                      - 2*ratio.bar[t,1]*sum(log.ratio[(t-n.tes-n.ref+1):(t-n.tes)])) )
        # Tes
        ratio.bar[t,2] <- mean(log.ratio[(t-n.tes+1):t])
        ratio.var[t,2] <- ( 1/max((sample.size-1),1)*(sum(log.ratio[(t-n.tes+1):t]^2)
                                                      + sample.size*ratio.bar[t,2]^2
                                                      - 2*ratio.bar[t,2]*sum(log.ratio[(t-n.tes+1):t])) )
  
        ratio.ref.ub <- ratio.bar[t,1]+sqrt(ratio.var[t,1])
        if (t >= (n.ref*2+n.tes) && ratio.bar[t,2] > ratio.ref.ub && quandts.stop == 0){
          # Start counting
          quandts.count <- quandts.count + 1
        }
        else {
          quandts.count <- 0
        }
  
        if (quandts.count == n.tes){
          quandts.count <- 0
          tao.hat[2,k] <- (t-n.tes+1)
          quandts.stop <- 1
          quandts.imp.count <- n.tes
        }
        # ###########
        # Quandts improved
        if (t >= (n.ref*2+n.tes) && ratio.bar[t,2] > ratio.ref.ub){
          # Start counting
          quandts.imp.count <- quandts.count + 1
        }
        else {
          quandts.imp.count <- 0
        }
  
        if (quandts.imp.count == n.tes){
          for (l in ((p+1):t)){
          quandts.imp.beta.hat <- beta.hat.mat[,(t-1)]
          quandts.imp.stud.res[t] <- (y[t]-data[t,]%*%quandts.imp.beta.hat)/sqrt((1 + data[t,]%*%solve(t(data[1:t,]) %*% data[1:t,])%*%data[t,]))
          }
  
          quandts.imp.stud.res.sq <- quandts.imp.stud.res^2
          quandts.imp.S <- sum((quandts.imp.stud.res.sq[(p+1):t]))
          quandts.imp.s.vec <- cumsum(quandts.imp.stud.res.sq)/quandts.imp.S
          quandts.imp.mean <- (c((p+1):t)-p)/(t-p)
          quandts.imp.stat <- abs(quandts.imp.s.vec[(p+1):t] - quandts.imp.mean)
          quandts.imp.stat.max <- max(quandts.imp.stat)
          quandts.imp.tao.hat <- as.numeric((which(quandts.imp.stat == quandts.imp.stat.max))[1])
          t.prime <- 0.5*(t-p)-1
          quandts.imp.c <- 1.3581015/sqrt(t.prime) -0.6701218/(t.prime) -0.8858694/(t.prime)^{3/2}
  
          if (quandts.imp.stat.max > quandts.imp.c){
            quandts.imp.count <- 0
            tao.hat[4,k] <- t - n.tes + 1
          }
  
          quandts.imp.count <- 0
        }
        
      }

      ###########
      # CUSUMSQ
      beta.hat <- beta.hat.mat[,(t-1)]
      # cusumsq.res <- y[1:i]-data[1:i,]%*%beta.hat
      stud.res[t] <- (y[t]-data[t,]%*%beta.hat)/sqrt((1 + data[t,]%*%solve(t(data[1:t,]) %*% data[1:t,])%*%data[t,]))


      ###########
  }
  # print(Sys.time()-start.time)
  # ###########
  # CUSUM
  # Keep commented out
  # if (cusum.count == 0){
  #   tao.hat[1, k] <- 0
  # }
  # Quandts
  if (quandts.stop == 0){
    tao.hat[2, k] <- 0
  }
  # Quandts improved
  # Keep commented out
  # if (quandts.imp.count == 0){
  #   tao.hat[4, k] <- 0
  # }
  # CUSUMSQ
  stud.res.sq <- stud.res^2
  S <- sum((stud.res.sq[(p+1):N]))
  s.vec <- cumsum(stud.res.sq)/S
  mean <- (c((p+1):N)-p)/(N-p)
  cusumsq.stat <- abs(s.vec[(p+1):N] - mean)
  cusumsq.stat.max <- max(cusumsq.stat)
  cusumsq.tao.hat <- as.numeric((which(cusumsq.stat == cusumsq.stat.max))[1])
  N.prime <- 0.5*(N-p)-1
  c <- 1.3581015/sqrt(N.prime) -0.6701218/(N.prime) -0.8858694/(N.prime)^{3/2}
  if (cusumsq.stat.max > c){
    tao.hat[3, k] <- cusumsq.tao.hat
  }
  if (cusumsq.stat.max < c || cusumsq.stat.max == c) {
    tao.hat[3, k] <- 0
  }
  ################
  # plot((p+1):N,cusum.stat[(p+1):N], ylim=c(-60,60))
  # abline((-3*a*p+a*N)/sqrt(N-p),2*a/sqrt(N-p))
  # abline((3*a*p-a*N)/sqrt(N-p),-2*a/sqrt(N-p))

  # c <- 1.3581015/sqrt(148) -0.6701218/(148) -0.8858694/(148)^{3/2}
  # plot((n.ref+n.tes):N, cusumsq.stat[(n.ref+n.tes):N], ylim=c(-1,1))
  # abline(c,0)
  # abline(-c,0)
  # # what the hell are the thresholds??
  # # The original ones for cusum are correct
  # # The ones for cusumsq are incorrect 
  plot((n.ref+n.tes):N, log.ratio[(n.ref+n.tes):N])
  points((n.ref+n.tes):N, ratio.bar[(n.ref+n.tes):N,1], col='red')
  lines((n.ref+n.tes):N, ratio.bar[(n.ref+n.tes):N,1]+sqrt(ratio.var[(n.ref+n.tes):N,1]), col='green')
  points((n.ref+n.tes):N, ratio.bar[(n.ref+n.tes):N,2], col='blue')
  
}

h <- hist(tao.hat, main = 'Chang point = 100, 200. Total no.=300', 
          xlab = 'Estimated cp')
print(Sys.time() - start.time)

######################################################
cusum.det <- sum(as.numeric(abs(cusum.stat) > ref.val))
print(paste(length(tao.hat),' change point detected. ', tao.hat))

plot(1:N,cusum.stat, ylim=c(-60,60))
abline((a*N-3*a*p)/sqrt(N-p), 2*a/sqrt(N-p))
abline((3*a*p-a*N)/sqrt(N-p), -2*a/sqrt(N-p))
print(Sys.time() - start.time)

######################################################
lb<-140
ub<-160
tao.hat0 <- matrix(rep(NA,300*4),nrow=4)
tao.hat1 <- matrix(rep(NA,300*4),nrow=4)
for (k in 1:4){
  tao.hat0[k,] <- tao.hat[k,1:300]
  tao.hat1[k,] <- tao.hat[k,(301):600]
  print(paste(length(subset(tao.hat0[k,], tao.hat0[k,]>=lb & tao.hat0[k,]<=ub)),
              length(subset(tao.hat1[k,], tao.hat1[k,]>=lb & tao.hat1[k,]<=ub))))
  print(paste(length(subset(tao.hat0[k,], tao.hat0[k,]!=0)),
              length(subset(tao.hat1[k,], tao.hat1[k,]!=0))))
}

for (j in 1:4){
weight.sum <- 0
sse <- 0
for (i in 1:300){
  tao.hat.val <- tao.hat0[j,i]
  weight.sum <- weight.sum + 1 + 1*as.numeric(tao.hat.val>150)
  sse <- sse + (1 + 1*as.numeric(tao.hat.val>150)) * (tao.hat.val - 150)^2
}
print(sse)
}
######################################################
# cusum.tao.hat1 <- tao.hat[1,1:N/2]
# cusum.tao.hat0 <- tao.hat[1,(N/2+1):N]
# quandt.tao.hat1 <- tao.hat[2,]
# quandt.tao.hat1 <- tao.hat[2,]
# cusumsq.tao.hat <- tao.hat[3,]
# quandt.imp.tao.hat <- tao.hat[4,]
# length((cusum.tao.hat1[(lb<=cusum.tao.hat1)])[cusum.tao.hat1[(lb<=cusum.tao.hat1)]<=ub])
# length((quandt.tao.hat[(lb<=quandt.tao.hat)])[quandt.tao.hat[(lb<=quandt.tao.hat)]<=ub])
# length((cusumsq.tao.hat[(lb<=cusumsq.tao.hat)])[cusumsq.tao.hat[(lb<=quandt.tao.hat)]<=ub])
# length((quandt.imp.tao.hat[(lb<=quandt.imp.tao.hat)])[quandt.imp.tao.hat[(lb<=quandt.imp.tao.hat)]<=ub])

########################
rls.fit <- rls.fast(response = y, data = data, lambda = 1,delta = 0.00001) 
               
# lev <- rls.fit$leverage.mat

for (i in (p+1):N){
  beta.hat <- beta.hat.mat[,(i-1)]
  res <- y[1:i]-data[1:i,]%*%beta.hat
  # sigma.sq.hat <- t(res)%*%res/(i-p)
  # stud.res[i] <- res[i]/sqrt((1 + lev[i,i])*sigma.sq.hat)
  stud.res[i] <- (y[i]-data[i,]%*%beta.hat)/sqrt((1 + data[i,]%*%solve(t(data[1:i,]) %*% data[1:i,])%*%data[i,]))
}
# cum.sum <- cumsum(stud.res)
# plot((p+1):N,cum.sum[(p+1):N])
# a <- 0.850
# abline(a*sqrt(N-p) - 2*a*p/sqrt(N-p), 2*a/sqrt(N-p))
# abline(-a*sqrt(N-p) - 2*a*p/sqrt(N-p), 2*a/sqrt(N-p))

stud.res.sq <- stud.res^2
S <- sum((stud.res.sq[(p+1):N]))
s.vec <- matrix(rep(NA, N),ncol=1)
c <- 0.10002
r <- 8
ref.line <- matrix(rep(NA, N*2), ncol=2)
for (j in (p+r):N) {
  s.vec[j] <- sum(stud.res.sq[(p+1):j])/S 
  ref.line[j,1] <- c + (j-p)/(N-p)
  ref.line[j,2] <- -c + (j-p)/(N-p)
}

ind<- (c((p+r):N)-p)/(N-p)
c <- 1.2238734/sqrt((p+r):N) - 0.6700069/((p+r):N) - 0.7351697/((p+r):N)^{3/2}
plot((p+r):N,(s.vec[(p+r):N]-ind),ylim=c(-1, 1))
lines((p+r):N, c)
lines((p+r):N, -c)


plot((p+r):N,s.vec[(p+r):N])
lines((p+r):N, ref.line[(p+r):N,1])
lines((p+r):N, ref.line[(p+r):N,2])
plot((p+1):N,(stud.res.sq[(p+1):N]-(c + (c((p+1):N)-p)/(N-p))) )


gof.summary <- (rls.fit$gof.summary)
stud.res <- gof.summary$Studentized.residuals
stud.res.diag <- diag(stud.res)
cum.sum <- cumsum(stud.res.diag[(p+1):N])
plot((p+1):N,cum.sum)

########################
# Plotting empirical cdf
n <- 10

realizations <- rnorm(n=n, mean=0, sd=1)
vec <- seq(from=(-2),to=2,by=0.1)
cdf.emp <- rep(NA, length(vec))
count.sum <- 0

for (i in 1:length(vec)){
  cdf.emp[i] <- sum(realizations<vec[i])
}
# cdf.emp[length(vec)] <- sum(realizations < vec[length(vec)])
cdf.emp <- cdf.emp/n
plot(pnorm(q=vec),cdf.emp)
plot(ecdf(realizations))

########################
# Outliers trial
N <- 20
p <- 2
intercept <- matrix(rep(1,N), ncol = 1)
data <- cbind(intercept, matrix(c(runif(n=N*(p-1), min=-10,max=10)),ncol=(p-1)))
beta <- matrix(rep(NA, N*(p-1)),ncol = (p-1))
beta <- t(cbind(rep(5,N),beta))
y <- matrix(rep(NA,N),ncol = 1)
# #########
beta[2,]<- -8
# beta[3,]<-3
#########
for (i in 1:N){
  y[i] <- data[i,]%*%beta[,i]
}
y <- y + matrix(c(rnorm(n=N, mean=0, sd=6)), ncol=1)

special.vec <- cbind(y[(N/2)],matrix(data[(N/2),],nrow=1))

# i
data[N/2,2] <- 0
y[N/2] <- -50

# ii
data[N/2,]<-cbind(1,max(data[,-1])*2)
y[N/2] <- data[N/2,] %*% beta[,1] * 0.8

# iii
data[N/2,]<-cbind(1,max(data[,-1])*3)
y[N/2] <- data[N/2,] %*% beta[,1] * 0.4

lambda<- 1
delta <-  0.00001
rls.fit <- rls(response = y, data = data, lambda = 1,delta = 0.00001, 
               cp = FALSE,n.tes=6,n.ref=12)
plot(data[,2], y, xlab='x', ylab='y',sub = '(i)', font.sub=2)
abline(rls.fit$beta_hat.vec[1], rls.fit$beta_hat.vec[2])
points(data[N/2,2], y[N/2], col='red')

rls.fit <- rls(response = y[-N/2], data = data[-N/2,], lambda = 1,delta = 0.00001, 
               cp = FALSE,n.tes=6,n.ref=12)
abline(rls.fit$beta_hat.vec[1], rls.fit$beta_hat.vec[2], col='red')

plot(rls.fit$y_hat.vec, rls.fit$leverage.vec, xlab='Fitted Values', ylab='Leverages')
points(rls.fit$y_hat.vec[N/2], rls.fit$leverage.vec[N/2],col='red', pch=16)
abline(2*p/N, 0, lty=2)

plot(rls.fit$y_hat.vec, rls.fit$Studentized.residuals[,N], xlab='Fitted Values', ylab='Studentized Residuals')
points(rls.fit$y_hat.vec[N/2], rls.fit$Studentized.residuals[N/2,N],col='red', pch=16)
abline(0, 0, lty=2)
plot(rls.fit$y_hat.vec, (rls.fit$y_hat.vec - y), xlab='Fitted Values', ylab='Residuals')
points(rls.fit$y_hat.vec[N/2], rls.fit$Studentized.residuals[N/2,N],col='red', pch=16)
abline(0, 0, lty=2)
plot(rls.fit$leverage.vec, rls.fit$Studentized.residuals[,N], xlab='Leverages', ylab='Studentized Residuals')
points(rls.fit$leverage.vec[N/2], rls.fit$Studentized.residuals[N/2,N],col='red', pch=16)
abline(0, 0, lty=2)

plot(rls.fit$y_hat.vec, rls.fit$Cook.distance[,N], xlab='Fitted Values', ylab='Residuals')
points(rls.fit$y_hat.vec[N/2], rls.fit$Cook.distance[N/2,N],col='red', pch=16)
abline(0, 0, lty=2)

data[N/2,] <- special.vec[-1]
y[N/2] <- special.vec[1]


################################################
################################################
# GK Algorithm
N <- 1000
response <- rnorm(N, 0, 1)
start.time <- Sys.time()
summary <- gka(response, 0.001)
gka.quantile(summary,phi = 0.5)
print(Sys.time() - start.time)
sort(response)[which(sort(response) == gka.quantile(summary,phi = 0.5))]

response.quantiles <- matrix(NA, nrow=100, ncol=1)
cdf.sample <- matrix(NA, nrow=N, ncol=1)
for (i in 1:100){
  response.quantiles[i] <- gka.quantile(summary,rank = 10*i)
}
for (i in 1:N){
  cdf.sample[i] <- sum(sort(response) < sort(response)[i])/N
}
plot(response.quantiles, seq(1,N,length.out=100)/N, type='s')
lines(sort(response),cdf.sample, col='red')

########################
# KG Test
N <- 1000
response <- rnorm(N, 0, 1)
gka.summary <- gka(response, 2/N)
kst(gka.summary, 0.1, 'pnorm', 0,1)
ks.test(x = response,y ='pnorm',0,1,alternative = 'two.sided',exact=TRUE)

########################
# CDF Plot
N <- sum(gka.summary[,2])
s <- nrow(gka.summary)
D.hat <- 0
F <- 'norm'
F.hat <- matrix(NA, nrow=s, ncol=1)
F.true <- matrix(NA, nrow=s, ncol=1)

for (i in 1:s){
  y <- gka.summary[i,1]
  # Find max. index j such that Y_j <= Y
  j <- which(gka.summary[,1] > y)[1] - 1
  if (is.na(j)){
    j <- s
  }
  # Find the approximate index of Y_j in the stream
  i_j.max <- sum(gka.summary[1:j,2]) + gka.summary[j,3]
  i_j.hat <- gka.quantile(gka.summary, ranking=i_j.max)
  i_j.hat <- i_j.max
  
  y_j <- gka.summary[j,1]
  if (F == 'norm'){
    F.true[i] <- pnorm(y_j, mean = 0, sd = 1)
  }
  else if (F == 'exp'){
    F.true[i] <- pexp(y_j, rate=1)
  }
  F.hat[i] <- i_j.hat/N
}
plot(gka.summary[,1], F.hat, type='s')
lines(gka.summary[,1], F.true, type='l', col='red')
plot(F.hat,F.true)
abline(0,1, col='red')
#####################################
# Space needed for the GKA plot
N <- 10000
summary.length <- matrix(NA, nrow=N,ncol=1)
response <- rnorm(N, 0, 1)

# Initial setting
summary <- data.frame(matrix(NA, ncol=3, nrow=1))
colnames(summary) <- c('v','g','delta')
s <- 0
# COMPRESS phase
for (n in 0:(length(response)-1)){
  v <- response[n+1]
  # Fill in the first 2 iterations
  if (n == 0){
    summary[1,] <- c(v,1,0)
    # Update size of summary
    s <- nrow(summary)
  }
  else{
    
    if (n %% (1/(2*epsilon)) == 0){
      
      i = s - 1
      while (i >= 2){
        
        j = i-1
        delta.i <- summary[i,3]
        g.sum <- sum(summary[j:i,2])
        v <- summary[i,1]
        
        while (j >= 2 && ((g.sum + delta.i) < 2*epsilon*n)){
          j <- j - 1
          g.sum <- g.sum + summary[j,2]
        }
        
        # Tune one index up
        j <- j + 1
        
        # DELETE phase
        if (j < i){
          # Merge tuples from j to i
          summary <- summary[-((j+1):i),]
          summary[j,] <- data.frame('v'=v, 'g'=g.sum-summary[(j-1),2], 'delta'=delta.i)
        }
        
        # Continue from the largest integer smaller than j
        i <- j - 1
        # Update size of the summary
        s <- nrow(summary)
      }
    }
    
    
    # INSERT phase
    v.0 <- summary[1,1]
    v.s_1 <- summary[s,1]
    
    # Extreme cases
    tuple.new <- data.frame('v'=NA, 'g'=NA, 'delta'=NA)
    if ( v < v.0 ){
      delta <- 0
      new.position <- 0
      summary <- rbind(tuple.new, summary)
    }
    else if ( v > v.s_1 ){
      delta <- 0
      new.position <- s
      summary <- rbind(summary, tuple.new)
    }
    else{
      # Find appropriate index i
      new.position <- which( v < summary[,1] )[1] - 1
      delta <- summary[new.position,2] + summary[new.position,3] - 1
      summary <- rbind(summary, tuple.new)
      summary[(new.position+2):(s+1), ] <- summary[(new.position+1):s, ]
    }
    
    # Insert new tuple
    tuple.new <- data.frame('v'=v, 'g'=1, 'delta'=delta)
    summary[(new.position+1),] <- tuple.new
    # Update size of summary
    s <- s + 1
  }
  summary.length[n+1] <- nrow(summary)
}
plot(1:N,summary.length,type='l')

########################
# Performance Analysis

N <- 60
response <- matrix(NA, nrow=N, ncol=1)
response <- rnorm(N, 0, 1)
response[1:(N/2)] <- rnorm(N/2, 0, 1)
response[(N/2+1):(N)] <- rnorm(N/2, 2, 1)
# response[1:(N/3)] <- rnorm(N/3, 0, 1)
# response[(N/3+1):(2*N/3)] <- rnorm(N/3, 2, 1)
# response[(2*N/3+1):(N)] <- rnorm(N/3, 0, 1)

kst.statistic <- matrix(NA, nrow=N,ncol=1)

# Initial setting
epsilon <- 2/N
gka.summary <- data.frame(matrix(NA, ncol=3, nrow=1))
colnames(gka.summary) <- c('v','g','delta')
s <- 0
# COMPRESS phase
kst.switch <- 0
for (n in 0:(length(response)-1)){
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


  kst.summary <- kst(gka.summary, 0.1, 'pnorm', 0,1)
  kst.statistic[n+1] <- kst.summary$statistic
  k <- 1.35810/sqrt(n) 
  if (kst.statistic[n+1] > k && kst.switch == 0){
    print(paste('cp detected at ', n))
    kst.switch <- 1
  }
}

plot(1:N, kst.statistic)
lines(1:N, 1.35810/sqrt(1:N), col='red')

########################
# Base and active states
N <- 60
response <- matrix(NA, nrow=N, ncol=1)
response <- rnorm(N, 0, 1)
response[1:(N/2)] <- rnorm(N/2, 0, 1)
response[(N/2+1):(N)] <- rnorm(N/2, 2, 1)
# response[1:(N/3)] <- rnorm(N/3, 0, 1)
# response[(N/3+1):(2*N/3)] <- rnorm(N/3, 2, 1)
# response[(2*N/3+1):(N)] <- rnorm(N/3, 0, 1)

kst.statistic <- matrix(NA, nrow=N,ncol=1)

# Initial setting
epsilon <- 2/N
gka.summary <- data.frame(matrix(NA, ncol=3, nrow=1))
colnames(gka.summary) <- c('v','g','delta')
s <- 0
# COMPRESS phase
kst.switch <- 0
for (n in 0:(length(response)-1)){
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
  
  
  kst.summary <- kst(gka.summary, 0.1, 'pnorm', 0,1)
  kst.statistic[n+1] <- kst.summary$statistic
  k <- 1.35810/sqrt(n) 
  if (kst.statistic[n+1] > k && kst.switch == 0){
    print(paste('cp detected at ', n))
    kst.switch <- 1
  }
}





#############################################
#############################################
# Real case analysis
# Outliers: c(163, 168, 230, 266, 294, 318, 583, 650, 728)
y <- matrix((read.csv('Chrom_13_GBM31.csv'))[,ncol(
            read.csv('Chrom_13_GBM31.csv'))],ncol=1)
data <- matrix(1, ncol=1, nrow=nrow(y))

plot(y[-c( 163, 168, 318, 728)])
rls.fit <- rls(response = y[-c( 163, 168, 318, 728)],
               data = data[-c( 163, 168, 318, 728)],lambda = 1,delta = 0.0001,
               cp = TRUE,n.tes = 12,n.ref = 12)
points(rls.fit$y_hat.vec,col='red')

###
y <- matrix((read.csv('Chrom_13_GBM31.csv')),ncol=7)
y <- matrix(unlist(y),ncol=7)
data <- matrix(1, ncol=1, nrow=nrow(y))
data <- cbind(data, y[,5])
y <- matrix(y[,7], ncol=1)
###

plot(y[-c( 163, 168, 318, 728)])
rls.fit <- rls(response = y[-c( 163, 168, 318, 728)],
               data = data[-c( 163, 168, 318, 728),],lambda = 1,delta = 0.0001,
               cp = TRUE,n.tes = 10,n.ref = 10)
points(rls.fit$y_hat.vec,col='red')

###
# data <- matrix(as.numeric(unlist((read.csv('Chrom_13_GBM31.csv'))[,c(5:7)])),ncol=3)
# data <- cbind(matrix(1,ncol=1,nrow=797),data)
rls.fit <- rls(response = y,
               data = data,lambda = 1,delta = 0.0001,
               cp = TRUE,n.tes = 12,n.ref = 12)

y.complete <- y
data.complete <- data
y <- matrix(y.complete[-c(  318)],ncol=1)
data <- matrix(data.complete[-c( 318),],ncol=ncol(data))

###########
N <- nrow(y)
p <- ncol(data)
# CUSUM
lambda<- 1
delta <-  0.00001
res <- matrix(rep(0, N*N), ncol=N)
cusum.stat <- matrix(rep(0, N), ncol=1)
cusum.count <- 0
a<- 0.9479
ref.val <- matrix(rep(0,N),ncol=1)
ref.val <- 2*a/sqrt(N-p) * (1:N) + (a*N-3*a*p)/sqrt(N-p)
###########
# Quandts
n.ref <- 10
n.tes <- 10
log.ratio <- matrix(rep(0,N),ncol=1)
ratio.bar <- matrix(rep(0,2*N),ncol=2)
ratio.var <- matrix(rep(0,2*N),ncol=2)
quandts.count <- 0
quandts.stop <- 0
###########
# Quandts improved
quandts.imp.stud.res <- matrix(rep(0, N),ncol=1)
quandts.imp.count <- 0
###########
# CUSUMSQ
rls.fit.sq <- rls.fast(response = y, data = data, lambda = 1,delta = 0.00001) 
beta.hat.mat <- rls.fit.sq$beta_hat.mat
stud.res <- matrix(rep(0, N),ncol=1)
###########
# KS Test


# start.time <- Sys.time()
for (t in (p+1):N){
  #############
  #CUSUM
  rls.fit <- rls.fast(y[1:t], data[1:t,], lambda, delta)
  cusum.res <- y[t] - (rls.fit$y_hat.vec)[t]
  cusum.res <- y[t] - data[t,]%*%rls.fit$beta_hat.vec
  cusum.lev <- data[t,]%*%solve(t(data[1:t,]) %*% data[1:t,])%*%data[t,]
  res[1,t] <- cusum.res/sqrt(1-cusum.lev)
  cusum.sigma.sq.hat <- (rls.fit$RSS)/(t-p)
  cusum.stat[t] <- sum(res[1,])/sqrt(cusum.sigma.sq.hat)
  # if (t >585){
  #   rls.fit <- rls.fast(y[582:t], data[582:t,], lambda, delta)
  #   cusum.res <- y[t] - (rls.fit$y_hat.vec)[t]
  #   cusum.res <- y[t] - data[t,]%*%rls.fit$beta_hat.vec
  #   cusum.lev <- data[t,]%*%solve(t(data[1:t,]) %*% data[1:t,])%*%data[t,]
  #   res[1,t] <- cusum.res/sqrt(1-cusum.lev)
  #   cusum.sigma.sq.hat <- (rls.fit$RSS)/(t-p)
  #   cusum.stat[t] <- sum(res[1,(582:t)])/sqrt(cusum.sigma.sq.hat)
  # }
  
  if (abs(cusum.stat[t]) > ref.val[t] && cusum.count == 0){
    tao.hat[1, k] <- t
    cusum.count <- 1
  }
  
  # ###########
  if (t > (n.ref+n.tes)){
    # Quandts
    rls.fit.ref <- rls.fast(y[(t-n.tes-n.ref+1):(t-n.tes)], data[(t-n.tes-n.ref+1):(t-n.tes),], lambda, delta)
    beta.hat.ref <- rls.fit.ref$beta_hat.vec
    res.ref <- y[(t-n.tes+1):t]-data[(t-n.tes+1):t,] %*% beta.hat.ref
    rss.ref <- t(res.ref)%*%res.ref
    
    rls.fit.tes <- rls.fast(y[(t-n.tes+1):t], data[(t-n.tes+1):t,], lambda, delta)
    rss.tes <- rls.fit.tes$RSS
    
    sigma.sq.hat <- as.numeric(rss.ref+rss.tes)/(n.tes+n.ref-p)
    sigma.sq.hat0 <- as.numeric(rss.ref)/(n.ref - p)
    sigma.sq.hat1 <- as.numeric(rss.tes)/(n.tes - p)
    
    # log.ratio[t] <- ((n.ref-n.tes)/2)*log(2*pi) + ((n.ref-n.tes)/2)*log(sigma.sq.hat) + (1/2)*rss.tes + (1/2)*rss.ref
    log.ratio[t] <- ((n.ref-n.tes)/2)*log(2*pi) + ((n.ref)/2)*log(sigma.sq.hat0) - 
                    ((n.tes)/2)*log(sigma.sq.hat1) + (1/2)*rss.tes + (1/2)*rss.ref
    
    sample.size <- min(t-n.ref-n.tes+1,n.ref+n.tes )
    # Ref
    ratio.bar[t,1] <- mean(log.ratio[(t-n.tes-n.ref+1):(t-n.tes)])
    ratio.var[t,1] <- ( 1/max((sample.size-1),1)*(sum(log.ratio[(t-n.tes-n.ref+1):(t-n.tes)]^2)
                                                  + sample.size*ratio.bar[t,1]^2
                                                  - 2*ratio.bar[t,1]*sum(log.ratio[(t-n.tes-n.ref+1):(t-n.tes)])) )
    # Tes
    ratio.bar[t,2] <- mean(log.ratio[(t-n.tes+1):t])
    ratio.var[t,2] <- ( 1/max((sample.size-1),1)*(sum(log.ratio[(t-n.tes+1):t]^2)
                                                  + sample.size*ratio.bar[t,2]^2
                                                  - 2*ratio.bar[t,2]*sum(log.ratio[(t-n.tes+1):t])) )
    
    ratio.ref.ub <- ratio.bar[t,1]+sqrt(ratio.var[t,1])
    if (t >= (n.ref*2+n.tes) && ratio.bar[t,2] > ratio.ref.ub && quandts.stop == 0 && t>50){
      # Start counting
      quandts.count <- quandts.count + 1
    }
    else {
      quandts.count <- 0
    }
    
    if (quandts.count == n.tes){
      quandts.count <- 0
      tao.hat[2,k] <- (t-n.tes+1)
      quandts.stop <- 1
      # quandts.imp.count <- n.tes
    }
    # ###########
    # Quandts improved
    if (t >= (n.ref*2+n.tes) && ratio.bar[t,2] > ratio.ref.ub){
      # Start counting
      quandts.imp.count <- quandts.imp.count + 1
    }
    else {
      quandts.imp.count <- 0
    }
    
    if (quandts.imp.count == n.tes){
      for (l in ((p+1):t)){
        quandts.imp.beta.hat <- beta.hat.mat[,(t-1)]
        quandts.imp.stud.res[t] <- (y[t]-data[t,]%*%quandts.imp.beta.hat)/sqrt((1 + data[t,]%*%solve(t(data[1:t,]) %*% data[1:t,])%*%data[t,]))
      }
      
      quandts.imp.stud.res.sq <- quandts.imp.stud.res^2
      quandts.imp.S <- sum((quandts.imp.stud.res.sq[(p+1):t]))
      quandts.imp.s.vec <- cumsum(quandts.imp.stud.res.sq)/quandts.imp.S
      quandts.imp.mean <- (c((p+1):t)-p)/(t-p)
      quandts.imp.stat <- abs(quandts.imp.s.vec[(p+1):t] - quandts.imp.mean)
      quandts.imp.stat.max <- max(quandts.imp.stat)
      quandts.imp.tao.hat <- as.numeric((which(quandts.imp.stat == quandts.imp.stat.max))[1])
      t.prime <- 0.5*(t-p)-1
      quandts.imp.c <- 1.3581015/sqrt(t.prime) -0.6701218/(t.prime) -0.8858694/(t.prime)^{3/2}
      
      if (quandts.imp.stat.max > quandts.imp.c){
        quandts.imp.count <- 0
        tao.hat[4,k] <- t - n.tes + 1
      }
      
      quandts.imp.count <- 0
    }
    
  }
  
  ###########
  # CUSUMSQ
  beta.hat <- beta.hat.mat[,(t-1)]
  # cusumsq.res <- y[1:i]-data[1:i,]%*%beta.hat
  stud.res[t] <- (y[t]-data[t,]%*%beta.hat)/sqrt((1 + data[t,]%*%solve(t(data[1:t,]) %*% data[1:t,])%*%data[t,]))
  
  
  ###########
}
if (quandts.stop == 0){
  tao.hat[2, k] <- 0
}
stud.res.sq <- stud.res^2
S <- sum((stud.res.sq[(p+1):N]))
s.vec <- cumsum(stud.res.sq)/S
mean <- (c((p+1):N)-p)/(N-p)
cusumsq.stat <- abs(s.vec[(p+1):N] - mean)
cusumsq.stat.max <- max(cusumsq.stat)
cusumsq.tao.hat <- as.numeric((which(cusumsq.stat == cusumsq.stat.max))[1])
N.prime <- 0.5*(N-p)-1
c <- 1.3581015/sqrt(N.prime) -0.6701218/(N.prime) -0.8858694/(N.prime)^{3/2}
if (cusumsq.stat.max > c){
  tao.hat[3, k] <- cusumsq.tao.hat
}
if (cusumsq.stat.max < c || cusumsq.stat.max == c) {
  tao.hat[3, k] <- 0
}
################
plot((p+1):N,cusum.stat[(p+1):N], ylim=c(-80,100), col='grey',xaxt='n',
     main='CUSUM', xlab='Time',yla='CUSUM statistic')
abline((-3*a*p+a*N)/sqrt(N-p),2*a/sqrt(N-p), col='red')
abline((3*a*p-a*N)/sqrt(N-p),-2*a/sqrt(N-p), col='red')
abline(v=768, col='blue',lty=2)
abline(v=527, col='blue',lty=2)
# cusum.stat2 <- cusum.stat[582:N]
points(588:N, cusum.stat2[588:N - 583],col='black')
points(cusum.stat[1:588],col='black')
axis(side=1,at=c(0,200,400,527,600,768))

c <- 1.3581015/sqrt(N.prime) -0.6701218/(N.prime) -0.8858694/(N.prime)^{3/2}
plot((n.ref+n.tes):N, cusumsq.stat[(n.ref+n.tes):N], ylim=c(-0.4,0.4),
     main='CUSUMSQ', xlab='Time', ylab='CUSUMSQ statistic')
abline(c,0,col='red',lty=2)
abline(-c,0,col='red',lty=2)
abline(v=524, col='blue', lty=2)
axis(side=1,at=c(0,200,400,524,600,800))

plot((n.ref+n.tes):N, log.ratio[(n.ref+n.tes):N], col='gray',
      main='Quandt\'s', xlab='Time', ylab='loglikelihood ratio',xaxt='n')
#points((n.ref+n.tes):N, ratio.bar[(n.ref+n.tes):N,1], col='dark green')
lines((n.ref+n.tes):N, ratio.bar[(n.ref+n.tes):N,1]+sqrt(ratio.var[(n.ref+n.tes):N,1]), col='red',lty=4)
points((n.ref+n.tes):N, ratio.bar[(n.ref+n.tes):N,2], col='blue')
abline(v='582', lty=2, col='blue')
axis(side=1,at=c(0,200,400,582,800))

print(tao.hat[,1])

y <- y.complete
data <- data.complete
# Results with outliers removed: c(527,580,522,580)

# Plots
plot(y, main='Chromosome 13 in GBM31 against time', xlab='Time',
     ylab='Chromosome 13', col='grey', pch=19, xaxt='n',
     cex.lab=1.5,cex.axis=1.5,cex.main=2)
lines(1:797,rbind(matrix(-0.25, nrow=length(1:539)),
                matrix(0.01357442, nrow=length(540:797))),col='red',lwd=2)
abline(v=540, col='blue', lty=2)
abline(0,0,lty=2)
axis(side=1,at=c(0,200,400,540,600,800),cex.axis=1.5)

# CUSUM, CUSUMSQ and Combined
# 527 580 522 580
tao.hat.val <- tao.hat[4,1]
plot(y, main='Combined', xlab='Time',
    ylab='Chromosome 13', col='grey', pch=19,xaxt='n',
    cex.lab=1.75,cex.axis=2,cex.main=3)
lines(1:797,rbind(matrix(-0.25, nrow=length(1:(tao.hat.val-1))),
                  matrix(0.01357442, nrow=length(tao.hat.val:797))),col='red',lwd=2)
abline(v=tao.hat.val, col='blue', lty=2)
abline(0,0,lty=2)
axis(side=1,at=c(0,200,400,tao.hat.val,600,800),cex.axis=2) 

# Quandt's
tao.hat.val <- tao.hat[2,1]
plot(y, main='Quandt\'s', xlab='Time',
     ylab='Chromosome 13', col='grey', pch=19, xaxt='n',
     cex.lab=1.75,cex.axis=2,cex.main=3)
lines(1:797,rbind(matrix(-0.223605, nrow=length(1:(159-1))),
                  matrix(-0.25, nrow=length(159:(tao.hat.val-1))),
                  matrix(0.01357442, nrow=length(tao.hat.val:797))),col='red',lwd=2)
abline(v=tao.hat.val, col='blue', lty=2)
abline(v=159, col='blue', lty=2)
abline(0,0,lty=2)
axis(side=1,at=c(0,159,200,400,tao.hat.val,800),cex.axis=2) 
