N <- 100
p <- 5
intercept <- matrix(rep(1,N), ncol = 1)
data <- cbind(intercept, matrix(c(runif(n=N*(p-1), min=-3,max=3)),ncol=p-1))
# y <- data %*% matrix(c(2),ncol = 1) + c(rnorm(n = N, mean=0, sd=1))
# data <- matrix(c(runif(n=N*p, min=-3,max=3)), ncol=p)
# beta.vec <- matrix(seq(1,2*N,2), ncol=1)^2
beta.vec <- matrix(c(1:p), ncol=1)
# beta.vec <- matrix(cbind(rep(2,N/2),rep(5,N/2)),ncol=1)
y <- data %*% beta.vec + matrix(c(rnorm(n=N, mean=0, sd=1)), ncol=1)

##################
rls.fit <- rls(response = y, data = data, lambda = 1,delta = 0.00001)
rls.fit$R2.vec[length(y)]
rls.fit$beta_hat.vec

system.time(
  rls.fit <- rls(response = y, data = data, lambda = 1,delta = 0.00001)
)
lev <- rls.fit$leverage.mat

plot(1:N, beta.vec)
points(1:N, rls.fit$beta_hat.mat[1,], col='red')
plot(data, y)
points(data, as.numeric(rls.fit$beta_hat.vec)*data,col='red')

##################
# Leverages
rls.fit <- rls(response = y, data = data, lambda = 1,delta = 0.00001)
lm.fit <- lm(y~data+0)

lev <- rls.fit$leverage.mat
tail(lev[,N])
tail(hatvalues(lm.fit))

##################
# Studentized residuals
rls.fit <- rls(response = y, data = data, lambda = 1,delta = 0.00001)
lm.fit <- lm(y~data+0)

res <- rls.fit$Studentized.residuals
tail(res[,N])
tail(lm.fit$residuals/sqrt((1-hatvalues(lm.fit))*as.numeric(t(lm.fit$residuals)%*%(lm.fit$residuals)/(N-p))))

##################
# Cook's distance
rls.fit <- rls(response = y, data = data, lambda = 1,delta = 0.00001)
lm.fit <- lm(y~data+0)

cooks <- rls.fit$Cook.distance
tail(cooks[,N])
tail(cooks.distance(lm.fit))

plot(1:N,cooks.distance(lm.fit))
points(1:N,cooks[,N],col='red')

##################
# Varying beta, discrete
N <- 100
p <- 2
intercept <- matrix(rep(1,N), ncol = 1)
data <- cbind(intercept, matrix(c(runif(n=N*(p-1), min=-3,max=3)),ncol=p-1))
beta <- matrix(rep(NA, N),ncol = 1)
beta <- t(cbind(rep(5,N),beta))
y <- matrix(rep(NA,N),ncol = 1)
k=2
for (i in 1:k){
  lb <- (N/k)*(i-1)+1
  ub <- (N/k)*i
  beta[2,lb:ub] <- rep(i,N/k)
}
for (i in 1:N){
  y[i] <- data[i,]%*%beta[,i]
}
y <- y + matrix(c(rnorm(n=N, mean=0, sd=1)), ncol=1)

par(mfrow=c(2,2))
lambda <- c(1,0.9,0.6,0.2)
# for (i in lambda){
#   rls.fit <- rls(response = y, data = data, lambda = lambda[i],delta = 0.00001)
#   beta.hat <- rls.fit$beta_hat.mat
#   plot(1:N,beta[2,])
#   points(1:N, beta.hat[2,],col='red')
# }
rls.fit <- rls(response = y, data = data, lambda = 1,delta = 0.00001)
beta.hat <- rls.fit$beta_hat.mat
plot(1:N,beta[2,])
points(1:N, beta.hat[2,],col='red')
rls.fit <- rls(response = y, data = data, lambda = 0.9,delta = 0.00001)
beta.hat <- rls.fit$beta_hat.mat
plot(1:N,beta[2,])
points(1:N, beta.hat[2,],col='red')
rls.fit <- rls(response = y, data = data, lambda = 0.6,delta = 0.00001)
beta.hat <- rls.fit$beta_hat.mat
plot(1:N,beta[2,])
points(1:N, beta.hat[2,],col='red')
rls.fit <- rls(response = y, data = data, lambda = 0.2,delta = 0.00001)
beta.hat <- rls.fit$beta_hat.mat
plot(1:N,beta[2,])
points(1:N, beta.hat[2,],col='red')

##################
# Change point: beta changes depending on covariates
N <- 100
p <- 2
intercept <- matrix(rep(1,N), ncol = 1)
data <- cbind(intercept, matrix(c(runif(n=N*(p-1), min=-10,max=10)),ncol=(p-1)))
beta <- matrix(rep(3, N),ncol = 1)
beta <- t(cbind(rep(2,N),beta))
cond <- as.numeric(data[,2]>1)
beta[2,] <- cond*t(beta[2,])
y <- matrix(rep(NA,N),ncol = 1)
for (i in 1:N){
  y[i] <- data[i,]%*%beta[,i]
}

phi.initial <- 1
phi <- change.pt(response = y, data = data, phi.initial = 1)

##########
N <- 100
p <- 2
intercept <- matrix(rep(1,N), ncol = 1)
data <- cbind(intercept, matrix(c(runif(n=N*(p-1), min=-10,max=10)),ncol=(p-1)))

beta <- matrix(rep(3, N),ncol = 1)
beta <- t(cbind(rep(2,N),beta))
cond <- as.numeric(data[,2]>(0))
beta[2,] <- cond*t(beta[2,])
y <- matrix(rep(NA,N),ncol = 1) 
for (i in 1:N){
  y[i] <- data[i,]%*%beta[,i]
}
y <- y + matrix(c(rnorm(n=N, mean=0, sd=1)), ncol=1)

phi.initial <- 8
data.large <- matrix(cbind(data, matrix(rep(NA, N*2), ncol=2)), ncol=(p+2))
phi <- phi.initial
for (k in 1:100){
  condition <- as.numeric(data.large[,2] > phi)
  u <- (data.large[,2] - phi)*condition
  v <- (-condition)
  data.large[,(p+1):(p+2)] <- matrix(c(u,v),ncol=2)
  beta.hat <- lm(y~(data.large[,-1]))$coefficients
  phi <-(as.numeric(beta.hat[p+2])/as.numeric(beta.hat[p+1]) + phi)
  print(phi)
}

##########
N <- 100
p <- 2
intercept <- matrix(rep(1,N), ncol = 1)
data <- cbind(intercept, matrix(c(runif(n=N*(p-1), min=-3,max=3)),ncol=p-1))
beta <- matrix(rep(NA, N),ncol = 1)
beta <- t(cbind(rep(5,N),beta))
y <- matrix(rep(NA,N),ncol = 1)
k=2
for (i in 1:k){
  lb <- (N/k)*(i-1)+1
  ub <- (N/k)*i
  beta[2,lb:ub] <- rep(i,N/k)
}
for (i in 1:N){
  y[i] <- data[i,]%*%beta[,i]
}

# y <- matrix(diag(data[1:N,]%*%beta[,1:N]),ncol=1)
y <- y + matrix(c(rnorm(n=N, mean=0, sd=1)), ncol=1)


likelihood <- -Inf
rss <- Inf
tao.hat <- p+1
tao.hat2<- p+1
lambda <- 1
delta <- 0.00001
for (t in (p+1):(N-p-1)){
  
  rls.fit1 <- rls.fast(y[1:t], data[1:t,], lambda, delta)
  rss1 <- rls.fit1$RSS
  
  rls.fit2 <- rls.fast(y[(t+1):N], data[(t+1):N,], lambda, delta)
  rss2 <- rls.fit2$RSS
  
  sigma.sq.hat <- as.numeric((rss1+rss2)/N)
  
  likelihood.new <- -(N/2)*log(2*pi) - (N/2)*log(sigma.sq.hat) - N/2
  # Equivalent as minimizing RSS
  
  if (likelihood.new > likelihood){
    tao.hat <- t
    likelihood <- likelihood.new
  }
}
print(tao.hat)
print(tao.hat2)

##########
# Change point detection
start.time <- Sys.time()
n.ref <- 10
n.tes <- 10
log.ratio <- matrix(rep(0,N),ncol=1)

lambda<- 1
delta <-  0.00001
ratio.bar <- matrix(rep(0,2*N),ncol=2)
ratio.var <- matrix(rep(0,2*N),ncol=2)
# ratio.sum <- 0
# Indicates when to start the next detection procedure
# 0 = Initial setting, detection confirmed when count = n.tes
count <- 0
for (t in (p+1):N){
  
  if (t >= (n.ref+n.tes)){
    rls.fit.ref <- rls.fast(y[(t-n.tes-n.ref+1):(t-n.tes)], data[(t-n.tes-n.ref+1):(t-n.tes),], lambda, delta)
    # rss.ref <- rls.fit.ref$RSS
    # could update the beta estimates directly
    beta.hat.ref <- rls.fit.ref$beta_hat.vec
    res.ref <- y[(t-n.tes+1):t]-data[(t-n.tes+1):t,] %*% beta.hat.ref
    rss.ref <- t(res.ref)%*%res.ref
    
    rls.fit.tes <- rls.fast(y[(t-n.tes+1):t], data[(t-n.tes+1):t,], lambda, delta)
    rss.tes <- rls.fit.tes$RSS
    
    sigma.sq.hat <- as.numeric(rss.ref+rss.tes)/(n.tes+n.ref-p)
    log.ratio[t] <- ((n.ref-n.tes)/2)*log(2*pi) + ((n.ref-n.tes)/2)*log(sigma.sq.hat) + (1/2)*rss.tes + (1/2)*rss.ref
  
    # MA and MV
    #######
    # sample.size <- min(t-n.ref-n.tes+1,n.ref+n.tes )
    # if (t == (n.ref+n.tes)){
    #   ratio.sum <- ratio.sum + log.ratio[t]
    # }
    # else{
    #   ratio.sum <- ratio.sum + log.ratio[t] - log.ratio[t-n.ref-n.tes]
    # }
    # ratio.bar[t] <- ratio.sum/sample.size
    # if ((t-n.ref-n.tes-1) > 0){
    #   ratio.var[t] <- 1/(sample.size-1)* ((sample.size-2)*ratio.var[t-1] 
    #                                       + (log.ratio[t] - ratio.bar[t-1])*(log.ratio[t] - ratio.bar[t]))
    #   ratio.var[t] <- ratio.bar[t] - 1/(sample.size-1) * (log.ratio[t-n.ref-n.tes+1] - ratio.bar[t])^2
    # }
    #######
    # ratio.bar[t] <- ratio.bar[t-1] + 1/(sample.size) * (log.ratio[t]-ratio.bar[t-1])
    # ratio.bar[t] <- ratio.bar[t] - (1/(sample.size)) * log.ratio[t-n.ref-n.tes+1] 
    # if ((t-n.ref-n.tes-1) > 0){
    # ratio.var[t] <- 1/(sample.size-1)* ((sample.size-2)*ratio.var[t-1] 
    #                 + (log.ratio[t] - ratio.bar[t-1])*(log.ratio[t] - ratio.bar[t]))
    # ratio.var[t] <- ratio.bar[t] - 1/(sample.size-1) * (log.ratio[t-n.ref-n.tes+1] - ratio.bar[t])^2
    # }
    #######
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
    if (t >= (n.ref*2+n.tes) && ratio.bar[t,2] > ratio.ref.ub){
      # Start counting
      count <- count + 1
    }
    else {
      count <- 0
    }
    
    if (count == n.tes-p){
      # Change point confirmed. Raise an alarm.
      print(c('Cp detected: i =',t))
      count <- 0
      # Reset initialization, D = delta^{-1}*I, ...
    }
  }
}
plot(1:N, log.ratio)
plot((n.ref+n.tes):N, log.ratio[(n.ref+n.tes):N])
points((n.ref+n.tes):N, ratio.bar[(n.ref+n.tes):N,1], col='red')
lines((n.ref+n.tes):N, ratio.bar[(n.ref+n.tes):N,1]+sqrt(ratio.var[(n.ref+n.tes):N,1]), col='green')
points((n.ref+n.tes):N, ratio.bar[(n.ref+n.tes):N,2], col='blue')
print(Sys.time() - start.time)
##################
# Change-pt detection (piecewise constant)
N <- 1000
p <- 3
intercept <- matrix(rep(1,N), ncol = 1)
data <- cbind(intercept, matrix(c(runif(n=N*(p-1), min=-10,max=10)),ncol=p-1))
beta <- matrix(rep(NA, N*(p-1)),ncol = (p-1))
beta <- t(cbind(rep(5,N),beta))
# Write up the case for p > 2
# beta[3,] <- rep(5,N)
y <- matrix(rep(NA,N),ncol = 1)
k=2
for (i in 1:k){
  lb <- (N/k)*(i-1)+1
  ub <- (N/k)*i
  beta[2,lb:ub] <- rep(i,N/k)
}
for (i in 1:N){
  y[i] <- data[i,]%*%beta[,i]
}
y <- y + matrix(c(rnorm(n=N, mean=0, sd=1)), ncol=1)


rls.fit <- rls(response = y, data = data, lambda = 1,delta = 0.00001, 
               cp = TRUE,n.tes=10,n.ref=10, alpha=3)
beta.hat <- rls.fit$beta_hat.mat
plot(1:N,beta[2,])
points(1:N, beta.hat[2,],col='red')

cp.summary <- rls.fit$cp.summary
log.ratio <- cp.summary$log.ratio
ratio.bar <- cp.summary$ratio.bar
ratio.var <- cp.summary$ratio.var
plot(1:N, log.ratio)
points(1:N,ratio.bar[,1], col='blue')
lines(1:N,ratio.bar[,1]+sqrt(ratio.var[,1]), col='purple',lty = 4)
points(1:N,ratio.bar[,2], col='red')

##################
# Change point detection (continuous)
N <- 300
p <- 2
intercept <- matrix(rep(1,N), ncol = 1)
data <- cbind(intercept, matrix(c(runif(n=N*(p-1), min=-10,max=10)),ncol=p-1))
beta <- matrix(rep(NA, N*(p-1)),ncol = (p-1))
beta <- t(cbind(rep(5,N),beta))
# Write up the case for p > 2
# beta[3,] <- rep(5,N)
y <- matrix(rep(NA,N),ncol = 1)
k=2
for (i in 1:k){
  lb <- (N/k)*(i-1)+1
  ub <- (N/k)*i
  beta[2,lb:ub] <- rep(i,N/k)
}
beta[2,51:N] <- c(2:51)
for (i in 1:N){
  y[i] <- data[i,]%*%beta[,i]
}
y <- y + matrix(c(rnorm(n=N, mean=0, sd=1)), ncol=1)


rls.fit <- rls(response = y, data = data, lambda = 1,delta = 0.00001, 
               cp = TRUE,n.tes=10,n.ref=10, alpha =5)
beta.hat <- rls.fit$beta_hat.mat
plot(1:N,beta[2,])
points(1:N, beta.hat[2,],col='red')



n.ref <- 10
n.tes <- 10
log.ratio <- matrix(rep(0,N),ncol=1)

lambda<- 1
delta <-  0.00001
ratio.bar <- matrix(rep(0,2*N),ncol=2)
ratio.var <- matrix(rep(0,2*N),ncol=2)
count <- 0
for (t in (p+1):N){
  
  if (t >= (n.ref+n.tes)){
    rls.fit.ref <- rls.fast(y[(t-n.tes-n.ref+1):(t-n.tes)], data[(t-n.tes-n.ref+1):(t-n.tes),], lambda, delta)
    beta.hat.ref <- rls.fit.ref$beta_hat.vec
    
    rls.fit.tes <- rls.fast(y[(t-n.tes+1):t], data[(t-n.tes+1):t,], lambda, delta)
    beta.hat.tes <- rls.fit.tes$beta_hat.vec
   
    tao.hat <- (beta.hat.tes[1]-beta.hat.ref[1])/as.numeric(beta.hat.ref[2]
                                                 - beta.hat.tes[2]) 
    if (tao.hat <= data[t,2] || tao.hat >= data[(t+1),2]){
      
    }
  }
}
plot(1:N, log.ratio)
plot((n.ref+n.tes):N, log.ratio[(n.ref+n.tes):N])
points((n.ref+n.tes):N, ratio.bar[(n.ref+n.tes):N,1], col='red')
lines((n.ref+n.tes):N, ratio.bar[(n.ref+n.tes):N,1]+sqrt(ratio.var[(n.ref+n.tes):N,1]), col='green')
points((n.ref+n.tes):N, ratio.bar[(n.ref+n.tes):N,2], col='blue')
print(Sys.time() - start.time)
##################
# Change-pt detection histgram plot
start.time <- Sys.time()
tao.hat <-as.numeric()
for (k in 1:300){
  N <- 300
  p <- 2
  intercept <- matrix(rep(1,N), ncol = 1)
  data <- cbind(intercept, matrix(c(runif(n=N*(p-1), min=-10,max=10)),ncol=p-1))
  beta <- matrix(rep(NA, N*(p-1)),ncol = (p-1))
  beta <- t(cbind(rep(5,N),beta))
  # beta[3,] <- rep(3,N)
  y <- matrix(rep(NA,N),ncol = 1)
  k=3
  for (i in 1:k){
    lb <- (N/k)*(i-1)+1
    ub <- (N/k)*i
    beta[2,lb:ub] <- rep(i,N/k)
  }
  for (i in 1:N){
    y[i] <- data[i,]%*%beta[,i]
  }
  y <- y + matrix(c(rnorm(n=N, mean=0, sd=1)), ncol=1)
  
  
  rls.fit <- rls(response = y, data = data, lambda = 1,delta = 0.00001, 
                 cp = TRUE,n.tes=6,n.ref=12)
  cp.summary <- rls.fit$cp.summary
  tao.hat <- cbind(tao.hat,as.numeric(cp.summary$tao.hat))
}
h <- hist(tao.hat, main = 'Chang point = 100, 200. Total no.=300', 
          xlab = 'Estimated cp')
print(Sys.time() - start.time)