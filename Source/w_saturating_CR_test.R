### Test the validity of my w parameterisation of saturating-FR CR model ###

require("rootSolve")
require("deSolve")
require("phaseR")

# Comparing functional response outputs -----------------------------------

f_classic <- function(R,k,a,h,mu){
  1/mu * alpha*a*R / (1+a*h*R)
}

f_new <- function(R,k,a,h,mu){
  x <- R/k
  w <- 1/mu * alpha*a*k / (1+a*h*k)
  xhalf <- 1 / (a*h*k)
  w * x * (xhalf+1) / (xhalf+x)
}

parms <- list(
  k <- 2,
  a <- 1,
  h <- 10,
  mu <- 0.1
)

Rseq <- 10^(seq(-3,3,length.out=100))

fout_classic <- with(parms, f_classic(R=Rseq,k,a,h,mu))
fout_new <- with(parms, f_new(R=Rseq,k,a,h,mu))

matplot(log(Rseq),cbind(fout_classic,fout_new),type="l")

# Dynamics ----------------------------------------------------------------

# Something wrong with handling time
# Not to do with suppression calculation, as this is correct

G = function(N,v,k){
  # v*(k/N - 1)
  v*(1 - N/k)
}

g = function(n){
  # 1/n - 1
  1 - n
}

# Feeding

f = function(n,b){
  (b + 1) / (b + n)
}

d = function(w,b){
  b / (w * (b + 1) - 1)
} 

parms1 <- list(
  v = 2,
  k = 2,
  a = 1,
  h = 0.5,
  alpha = 0.25,
  mu = 0.1
)

dx_dt1 <- function(t,y,parms){
  with(parms,{
    N <- exp(y)
    dx1 <- G(N[1],v,k) - a * N[2] / (1 + a * h * N[1])
    dx2 <- alpha * a * N[1] / (1 + a * h * N[1]) - mu
    list(c(dx1=dx1,dx2=dx2))
  })
}

x01 <- c(x1=0,x2=0)
tT <- 10^4
tmax <- 10^2
tseq1 <- seq(0,tmax,length.out=tT)
lnxt1 <- ode(y=x01,times=tseq1,func=dx_dt1,parms=parms1)[,-1] 

par(mfrow=c(1,1))
matplot(tseq1,lnxt1,type="l",lty=2)

Cstar <- runsteady(y=x01,time=c(0,Inf),func=dx_dt1,parms=parms1)$y[2]
lNstar <- c(log(parms1$k),Cstar)
  # only works if we don't have a limit cycle

parms2 <- with(parms1, list(
  w = (1/mu) * alpha * a * k / (1 + a * h * k),
  b = 1 / (a * h * k),
  p = mu / v
))

with(c(parms1,parms2), abline(h=log(k)+log(d(w,b)),col="blue"))
  # suppression (d) is correct - matches observed ndotdot from simulation

dx_dt2 <- function(t,y,parms){
  with(parms,{
    n <- exp(y)
    ndd <- d(w,b)                               
    dx1 <- g(n[1]) - g(ndd) / f(ndd,b) * f(n[1],b) * n[2] # loss term 1/2, should be 2/3
    dx2 <- p * (w * f(n[1],b) * n[1] - 1)       
    list(c(dx1=dx1,dx2=dx2))
  })
}

x02 <- x01 - lNstar
tseq2 <- seq(0,tmax,length.out=tT) * parms1$v # time in tau units 
lnxt2 <- ode(y=x02,times=tseq2,func=dx_dt2,parms=c(parms1,parms2))[,-1] + rep(lNstar,each=tT)

matplot(tseq1,lnxt2,type="l",lty=3,lwd=3,add=T)
# tseq1 in both cases because want to plot on same (original) timescale
