### Testing whether my R0 parameterisation gives same dynamics as traditional C-R model ###

require("rootSolve")
require("deSolve")
require("phaseR")

source("Source/consumer_resource_functions.R")

# Logistic ----------------------------------------------------------------

zparms <- list(
  sigma = c(0,0), # climate sensitivity
  zmu = 0,
  zsig = 0,
  zl = 1000 # wavelength
)

parms1 <- list(
  v = 3,
  k = 2,
  a = 3,
  alpha = 0.5,
  mu = 2
)

parms2 <- with(parms1, list(
  R0 = alpha * a / mu * k,
  p = mu / (k*v)
))

dx_dt1 <- function(t,y,parms){
  with(parms,{
    z <- zt_cyclic(t,parms)
    N <- exp(y)
    dx1 <- v * (k - N[1]) - a * N[2]    # + sigma[1] * z
    dx2 <- alpha * a * N[1] - mu        # + sigma[2] * z
    list(c(dx1=dx1,dx2=dx2))
  })
}

dx_dt2 <- function(t,y,parms){
  with(parms,{
    z <- zt_cyclic(t,parms) # would need to change this for tau
    N <- exp(y)
    dx1 <- 1 - N[1] - (R0-1)/R0 * N[2]  # + sigma[1] * z
    dx2 <- p * (R0 * N[1] - 1)          # + sigma[2] * z
    list(c(dx1=dx1,dx2=dx2))
  })
}

Rstar <- with(parms1, k)
Cstar <- with(parms1, v/a*(k-mu/(alpha*a)))

x01 <- c(x1=0,x2=0)
x02 <- x01 - log(c(Rstar,Cstar))
tT <- 10^4
tmax <- 10^2
tseq1 <- seq(0,tmax,length.out=tT)
tseq2 <- seq(0,tmax,length.out=tT) * parms1$v * parms1$k # time in tau units

lnxt1 <- ode(y=x01,times=tseq1,func=dx_dt1,parms=c(parms1,zparms))[,-1] 
lnxt2 <- ode(y=x02,times=tseq2,func=dx_dt2,parms=c(parms2,parms1,zparms))[,-1] + rep(log(c(Rstar,Cstar)),each=tT)

par(mfrow=c(1,1))
matplot(tseq1,lnxt1,type="l",lty=2)
matplot(tseq1,lnxt2,type="l",lty=3,lwd=3,add=T)
  # tseq1 in both cases because want to plot on same (original) timescale

# Abiotic -----------------------------------------------------------------

zparms <- list(
  sigma = c(0,0), # climate sensitivity
  zmu = 0,
  zsig = 0,
  zl = 1000 # wavelength
)

parms1 <- list(
  v = 2,
  k = 2,
  a = 2,
  b = 1.5,
  alpha = 0.5,
  beta = 0.4,
  mu = 0.05,
  m = 0.1
)

Rstar <- with(parms1, k)
Cstar <- with(parms1, v * (k*alpha/mu - 1/a))
Cdd <- with(parms1, m/(beta*b))
Rddd <- with(parms1, k/(a/v*(Cdd)+1))
Pstar <- with(parms1, (alpha*a*Rddd - mu)/b)
logstar <- log(c(Rstar,Cstar,Pstar))

parms2 <- with(parms1, list(
  R0 = alpha * a / mu * Rstar,
  S0 = beta * b / m * Cstar,
  p = mu / v,
  q = m / mu
))

log(c(Rddd,Cdd,1)) - logstar
with(parms2,log(c(
  1/(1+1/S0*(R0-1)),
  1/S0,
  1
)))
# check on analytical calculations of equilibria

dx_dt1 <- function(t,y,parms){
  with(parms,{
    z <- zt_cyclic(t,parms)
    N <- exp(y)
    dx1 <- v * (k/N[1] - 1) - a * N[2]  # + sigma[1] * z
    dx2 <- alpha * a * N[1] - mu - b*N[3]       # + sigma[2] * z
    dx3 <- beta * b * N[2] - m
    list(c(dx1=dx1,dx2=dx2,dx3=dx3))
  })
}

dx_dt2 <- function(t,y,parms){
  with(parms,{
    z <- zt_cyclic(t,parms) # would need to change this for tau
    N <- exp(y)
    dx1 <- 1/N[1] - 1 - (R0-1) * N[2]  # + sigma[1] * z
    dx2 <- p * (R0 * N[1] - 1 - (R0/(1+(R0-1)/S0) - 1) * N[3]) # + sigma[2] * z
    dx3 <- p * q * (S0 * N[2] - 1)
    list(c(dx1=dx1,dx2=dx2,dx3=dx3))
  })
}

x01 <- c(x1=0,x2=0,x3=0)
x02 <- x01 - logstar
tT <- 10^4
tmax <- 10^2
tseq1 <- seq(0,tmax,length.out=tT)
tseq2 <- seq(0,tmax,length.out=tT) * parms1$v # time in tau units

lnxt1 <- ode(y=x01,times=tseq1,func=dx_dt1,parms=c(parms1,zparms))[,-1] 
lnxt2 <- ode(y=x02,times=tseq2,func=dx_dt2,parms=c(parms2,zparms))[,-1] + rep(logstar,each=tT)

par(mfrow=c(1,1))
matplot(tseq1,lnxt1,type="l",lty=2)
matplot(tseq1,lnxt2,type="l",lty=3,lwd=3,add=T)
# tseq1 in both cases because want to plot on same (original) timescale

# C-2R --------------------------------------------------------------------

### Logistic

zparms <- list(
  sigma = c(0,0), # climate sensitivity
  zmu = 0,
  zsig = 0,
  zl = 1000 # wavelength
)

parms1 <- list(
  v1 = 1,
  v2 = 1,
  k1 = 1,
  k2 = 1,
  a1 = 0.3,
  a2 = 0.25,
  alpha1 = 0.5,
  alpha2 = 0.5,
  mu = 0.1
)

# dx_dt1 <- function(t,y,parms){
#   with(parms,{
#     z <- zt_cyclic(t,parms)
#     N <- exp(y)
#     dx1 <- v1 * (k1/N[1] - 1) - a1 * N[3]  # + sigma[1] * z
#     dx2 <- v2 * (k2/N[2] - 1) - a2 * N[3]  # + sigma[1] * z
#     dx3 <- alpha1 * a1 * N[1] + alpha2 * a2 * N[2] - mu
#     list(c(dx1=dx1,dx2=dx2,dx3=dx3))
#   })
# }

dx_dt1 <- function(t,y,parms){
  with(parms,{
    z <- zt_cyclic(t,parms)
    N <- exp(y)
    dx1 <- v1 * (1 - N[1]/k1) - a1 * N[3]  # + sigma[1] * z
    dx2 <- v2 * (1 - N[2]/k2) - a2 * N[3]  # + sigma[1] * z
    dx3 <- alpha1 * a1 * N[1] + alpha2 * a2 * N[2] - mu
    list(c(dx1=dx1,dx2=dx2,dx3=dx3))
  })
}

x01 <- c(x1=0,x2=0,x3=0)
tT <- 10^4
tmax <- 10^3
tseq1 <- seq(0,tmax,length.out=tT)

lnxt1 <- ode(y=x01,times=tseq1,func=dx_dt1,parms=c(parms1,zparms))[,-1]

logstar <- lnxt1[tT,]
  # take final simulation value as N* estimate

Rstar1 <- with(parms1, k1)
Rstar2 <- with(parms1, k2)

parms2 <- with(parms1, list(
  w1 = alpha1 * a1 / mu * Rstar1,
  w2 = alpha2 * a2 / mu * Rstar2,
  q = v2 / v1, # sets v1 as base timescale?
  p = mu / v1 # p = mu / v1
))

Rdd_dash <- function(w1,w2){
  (w2/w1 + 1/w2 - 1) / (w2/w1 + w1/w2)
}

with(parms2, log(Rdd_dash(w1,w2)))
with(parms2, log(Rdd_dash(w2,w1)))
  # unlogged values give percentage of k

# dx_dt2 <- function(t,y,parms){
#   with(parms,{
#     z <- zt_cyclic(t,parms) # would need to change this for tau
#     N <- exp(y)
#     dx1 <- 1/N[1] - 1 - lossf(w1,w2) * N[3] 
#     dx2 <- q * (1/N[2] - 1 - lossf(w2,w1) * N[3])
#       # does q include the loss term?
#     dx3 <- p * (w1 * N[1] + w2 * N[2] - 1)
#     list(c(dx1=dx1,dx2=dx2,dx3=dx3))
#   })
# }

dx_dt2 <- function(t,y,parms){
  with(parms,{
    z <- zt_cyclic(t,parms) # would need to change this for tau
    N <- exp(y)
    dx1 <- 1 - N[1] - (1 - Rdd_dash(w1,w2)) * N[3] 
    dx2 <- q * ( 1 - N[2] - (1 - Rdd_dash(w2,w1)) * N[3] )
      # does q include the loss term?
    dx3 <- p * (w1 * N[1] + w2 * N[2] - 1)
      # may need a v2 term in here too?
    list(c(dx1=dx1,dx2=dx2,dx3=dx3))
  })
}

logdot <- with(parms1,c(log(k1),log(k2),logstar[3]))
x02 <- x01 - logdot
tseq2 <- seq(0,tmax,length.out=tT) * parms1$v1 # time in tau units

lnxt2 <- ode(y=x02,times=tseq2,func=dx_dt2,parms=c(parms2,zparms))[,-1] + rep(logdot,each=tT)

par(mfrow=c(1,1))
matplot(tseq1,lnxt1,type="l",lty=2)
matplot(tseq1,lnxt2,type="l",lty=3,lwd=3,add=T)
# tseq1 in both cases because want to plot on same (original) timescale

# Varying perturbation fluctuation speed ----------------------------------
# (from previous code)
lzlmin <- -1
lzlmax <- 3
nzl <- 20
zlseq <- 10^seq(lzlmin,lzlmax,length.out=nzl)
tburn <- table(tseq<(tmax/10))[1] # first 10% of timepoints dropped

parms_temp <- parms
lnx_sd <- matrix(nr=nzl,nc=2) # two species

for(i in 1:nzl){
  parms_temp$zl <- zlseq[i]
  lnxt_temp <- ode(y=x0,times=tseq,func=dx_dt,parms=parms_temp) 
  lnx_sd[i,] <- apply(lnxt_temp[-(1:tburn),-1],2,sd)
}

matplot(log10(zlseq),lnx_sd,type="l")
matplot(log10(zlseq),apply(lnx_sd,1,diff),type="l")

