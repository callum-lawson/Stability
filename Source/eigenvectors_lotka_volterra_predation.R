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
tseq <- seq(0,10^3,length.out=tT)

lnxt <- ode(y=x0,times=tseq,func=dx_dt,parms=parms) 
matplot(lnxt[,1],lnxt[,-1],type="l")