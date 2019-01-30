### Exploring closed-system steady states ###

require("deSolve")
require("rootSolve")
require("phaseR")

recycle <- function(t,y,parms){
  with(parms,{
    R <- y[1]
    C <- y[2]
    dR <- 1 - C - (1 + theta + beta*theta*C)*R
    dC <- (alpha*beta*theta*R - omega*theta)*C
    list(v*c(dR,dC))
  })
}

y0 <- rep(1/3,2) # rep(1/3,3)
tseq <- seq(0,10,length.out=1000)

RCstarana <- function(parms){
  with(parms, {
    R <- omega/(alpha*beta)
    C <- (1-R*(1+theta))/(1+theta*beta*R)
    c(R,C)
  })
}

parms <- list(
  v=100,
  theta=1000,
  omega=0.001,
  alpha=1,
  beta=10
  )

# parms <- list(
#   v=1,
#   theta=1000,
#   omega=1/4.307801,
#   alpha=0.621764,
#   beta=1000*0.05489633
# )

# trial <- ode(y=y0,times=tseq,func=recycle,parms=parms)
# matplot(trial[,-1],type="l",ylim=c(0,max(ystar)))
# ( ystar <- steady(y=y0,parms=parms,fun=recycle,times=c(0,Inf),method="runsteady")$y )
( ystar <- RCstarana(parms) )

cxlim <- c(0,2*ystar[1])
cylim <- c(0,2*ystar[2])

par(mfrow=c(1,1))
flowField(recycle,xlim=cxlim,ylim=cylim,parameters=parms,points=30,add=FALSE)
nullclines(recycle,xlim=cxlim,ylim=cylim,parameters=parms,points=100,add.legend=FALSE)

r_recycle <- function(C,parms){
  with(parms, {
    alpha * beta * (1-C) / (1 + 1/theta + beta * C) - theta * omega
  })
}

curve(r_recycle(x,parms))
abline(h=0,col="red",lty=3)
curve(r_recycle(10^(x),parms),xlim=c(-3,0))
abline(h=0,col="red",lty=3)

library(RColorBrewer)
ntheta <- 50
thetacols <- heat.colors(ntheta) # brewer.pal(ntheta,"PuOr")
thetaseq <- 10^(seq(0,3,length.out=ntheta))
newparms <- list()
for(i in 1:ntheta){
  newparms[[i]] <- parms
  newparms[[i]]$theta <- thetaseq[i]
  curve(r_recycle(10^(x),newparms[[i]]),add=T,col=thetacols[i])
}
    
# Explicit B --------------------------------------------------------------

explicit <- function(t,y,parms){
  with(parms,{
    R <- y[1]
    B <- y[2]
    dR <- B - theta*R
    dB <- theta*R * B
    list(v*c(dR,dB))
  })
}

newparms <- list(
  v = 1,
  theta = 1
)

flowField(explicit,xlim=c(0,1),ylim=c(0,1),parameters=newparms,points=30,add=FALSE)
nullclines(explicit,xlim=c(0,1),ylim=c(0,1),parameters=newparms,points=100,add.legend=FALSE)
abline(1,-1,col="red")
