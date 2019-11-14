### Integrate dynamics for LV predator-prey model with Sine wave perturbations ###

require("rootSolve")
require("deSolve")
require("phaseR")

source("Source/consumer_resource_functions.R")

# Setup -------------------------------------------------------------------

parms <- list(
  r = 2,
  p = 0.5,
  alpha = c(1,0),
  zmu = 0,
  zsig = 0.1,
  zl = 1000 # wavelength
)

dx_dt <- function(t,y,parms){
  with(parms,{
    z <- zt_cyclic(t,parms)
    N <- exp(y)
    dx1 <- 1 - N[1] - r * N[2] + alpha[1] * z
    dx2 <- p * (r * N[1] - 1)  + alpha[2] * z
    list(c(dx1=dx1,dx2=dx2))
  })
}

x0 <- c(dx1=0,dx2=0)
tT <- 10^4
tmax <- 10^3
tseq <- seq(0,tmax,length.out=tT)

par(mfrow=c(1,2))
lnxt <- ode(y=x0,times=tseq,func=dx_dt,parms=parms) 
matplot(tseq,lnxt[,-1],type="l")
plot(lnxt[,2],lnxt[,3])

# Varying perturbation fluctuation speed ----------------------------------

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

