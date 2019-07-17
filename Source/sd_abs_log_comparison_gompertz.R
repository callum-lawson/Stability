### Is it better to measure N variability on absolute or log scale? ###

# Simulations -------------------------------------------------------------

v <- 1
k <- 10
sig <- log(2) # average 2-fold variation in growth rate
tT <- 10^6
X0 <- 0
eps <- rnorm(tT,0,sig)
Xt_gompertz <- Xt_logistic <- vector("numeric",length=tT)
Xt_gompertz[1] <- Xt_logistic[1] <- X0

dX_dt_gompertz <- function(X,v,k,eps){
  v * (log(k) - X) + eps
}

dX_dt_logistic <- function(X,v,k,eps){
  v * (k - exp(X)) + eps
}

for(t in 2:tT){
  Xt_gompertz[t] <- Xt_gompertz[t-1] + dX_dt_gompertz(Xt_gompertz[t-1],v,k,eps[t])
  Xt_logistic[t] <- Xt_logistic[t-1] + dX_dt_logistic(Xt_logistic[t-1],v,k,eps[t])
}

# Population size ---------------------------------------------------------

### HISTORGRAMS

par(mfrow=c(2,2),mar=c(4.5,4.5,1.5,1.5))

hist(Xt_gompertz,breaks=1000)
abline(v=log(k),lty=3,col="red")
hist(exp(Xt_gompertz),breaks=1000)
abline(v=k,lty=3,col="red")
hist(Xt_logistic,breaks=1000)
abline(v=log(k),lty=3,col="red")
hist(exp(Xt_logistic),breaks=1000)
abline(v=k,lty=3,col="red")

### K ESTIMATION 

# log scale  
exp(mean(Xt_gompertz)) # unbiased
exp(mean(Xt_logistic)) # biased

# absolute scale 
mean(exp(Xt_gompertz)) # biased
mean(exp(Xt_logistic)) # unbiased

### SD ESTIMATION

# log scale 
exp(sd(Xt_gompertz)) 
exp(sd(Xt_logistic))

# absolute scale
sd(exp(Xt_gompertz))
sd(exp(Xt_logistic))

# Population growth -------------------------------------------------------

Rt_gompertz <- diff(Xt_gompertz)
Rt_logistic <- diff(Xt_logistic)

### HISTORGRAMS

par(mfrow=c(2,2))
hist(Rt_gompertz,breaks=1000)
abline(v=0,lty=3,col="red")
hist(exp(Rt_gompertz),breaks=1000)
abline(v=1,lty=3,col="red")
hist(Rt_logistic,breaks=1000)
abline(v=0,lty=3,col="red")
hist(exp(Rt_logistic),breaks=1000)
abline(v=1,lty=3,col="red")

### K ESTIMATION 

# log scale  
exp(mean(Rt_gompertz)) # unbiased
exp(mean(Rt_logistic)) # unbiased

# absolute scale 
mean(exp(Rt_gompertz))
mean(exp(Rt_logistic))
  # biased measure, shouldn't use (Lewontin & Cohen)

### SD ESTIMATION

# log scale 
exp(sd(Rt_gompertz)) 
exp(sd(Rt_logistic))

# absolute scale
sd(exp(Rt_gompertz))
sd(exp(Rt_logistic))
