### Effects of refuges on attack rates and functional response ###

ess <- function(R,r,k,xi,c){
  p <- 1/(1+c)*(c-xi*k/(2*r*R))
  p[p<=0] <- 0
  p[p>=1] <- 1
  return(p)
}

drdtheta <- function(R,r,k,xi,c,theta){
  2*r*R*(theta*(1+c)-c)/k - xi
}

Rstar <- function(r,k,xi,c){
  xi*k/(2*r*c)
}

dR_Rdt <- function(R,r,k,xi,c){
  r*(1-c*R/k) + xi
}

inflow <- function(r,k,xi,c){
  R_A <- Rstar(r,k,xi,c)
  dR_Rdt(R=R_A,r,k,xi,c)
}

thetastar <- function(c){
  c/(1+c)
}

curve(ess(R=x,r=5,k=100,xi=1,c=1),xlim=c(0,100),ylim=c(0,1))
curve(ess(R=x,r=5,k=100,xi=2,c=2),add=T,lty=2)
curve(ess(R=x,r=5,k=100,xi=4,c=4),add=T,lty=3)
curve(ess(R=x,r=5,k=100,xi=10,c=10),add=T,lty=4)
curve(ess(R=x,r=5,k=100,xi=100,c=100),add=T,lty=4)
curve((x-10)/x,add=T,col="red")

Rstar(r=5,k=100,xi=1,c=1)
curve(inflow(r=5,k=100,xi=1*x,c=1*x),xlim=c(0,100))
thetastar(c=100)

f <- function(R,C,r,k,xi=1,c=1){
  theta <- ess(R,r,k,xi*C,c)
  theta * R
}

curve(f(R=x,C=1,r=5,xi=1,k=100),xlim=c(0,100))
curve(f(R=x,C=2,r=5,xi=1,k=100),add=T,lty=2)
curve(f(R=x,C=5,r=5,xi=1,k=100),add=T,lty=3)

dR_Rdt_bar <- function(R,r,k,xi,c){
  theta <- ess(R,r,k,xi,c)
  (1-theta) * dR_Rdt((1-theta)*R,r,k,xi,c) + theta * dR_Rdt(theta*R,r,k,xi=0,c=1)
}

curve(dR_Rdt_bar(R=10^x,r=5,k=100,xi=1,c=5),xlim=c(-2,2),ylim=c(1,8),n=10^3)
curve(dR_Rdt_bar(R=10^x,r=5,k=100,xi=2,c=1),add=T,lty=2,n=10^3)
curve(dR_Rdt_bar(R=10^x,r=5,k=100,xi=2,c=5),add=T,lty=3,n=10^3)
curve(dR_Rdt_bar(R=10^x,r=5,k=100,xi=0,c=1),add=T,lty=4,n=10^3)

nAstar <- function(nB,r,k,xi,c){
  1/c*(xi*k/r+nB)
}
curve(nAstar(nB=x,r=5,k=100,xi=1,c=1),xlim=c(0,100),ylim=c(0,100),n=10^3)

irate <- function(nB,r,k,xi,c){
  nA <- nAstar(nB,r,k,xi,c)
  dR_Rdt(nA,r,k,xi,c) * nA
}
curve(irate(nB=x,r=5,k=100,xi=1,c=1),xlim=c(0,100),n=10^3)
curve(irate(nB=x,r=0.75*5,k=0.75*100,xi=1,c=1),add=T,lty=2)
curve(irate(nB=x,r=0.5*5,k=0.5*100,xi=1,c=1),add=T,lty=3)
curve(irate(nB=x,r=5,k=100,xi=0,c=1),add=T,col="red")
  # but this assumes that all individuals produced by A will move to B
  # how to siphon-off the individuals that stay in A, 
  # (as well as those that move from B->A)?
