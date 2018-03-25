### Effects of refuges on attack rates and functional response ###

# Functions ---------------------------------------------------------------

drdtheta <- function(R,r,k,xi,c,theta){
  2*r/k*(c(1-theta)-theta)*R - xi
}

theta <- function(R,r,k,xi,c){
  p <-  c/(c+1) * (1-xi*k/(2*r*c*R)) 
  p[p<=0] <- 0
  p[p>=1] <- 1
  return(p)
}

theta2 <- function(R,r,k,xi,c){
  p <-  c/(c+1) * (1-xi*k/(r*c*R)) 
  p[p<=0] <- 0
  p[p>=1] <- 1
  return(p)
}

Rstar <- function(r,k,xi,c){
  xi*k/(2*r*c)
}

dR_Rdt <- function(R,r,k,xi,c){
  r*(1-c*R/k) - xi
}

inflow <- function(r,k,xi,c){
  R_A <- Rstar(r,k,xi,c)
  dR_Rdt(R=R_A,r,k,xi,c)
}

thetastar <- function(c){
  c/(1+c)
}

f <- function(R,C,r,k,xi=1,c=1){
  theta <- theta(R,r,k,xi*C,c)
  theta * R
}

dR_Rdt_bar <- function(R,r,k,xi,c){
  theta <- theta(R,r,k,xi,c)
  (1-theta) * dR_Rdt((1-theta)*R,r,k,xi=0,c) + theta * dR_Rdt(theta*R,r,k,xi,c=1)
}

dR_Rdt_bar2 <- function(R,r,k,xi,c){
  theta <- theta2(R,r,k,xi,c)
  (1-theta) * dR_Rdt((1-theta)*R,r,k,xi=0,c) + theta * dR_Rdt(theta*R,r,k,xi,c=1)
}

nAstar <- function(nB,r,k,xi,c){
  1/c*(xi*k/r+nB)
}

irate <- function(nB,r,k,xi,c){
  nA <- nAstar(nB,r,k,xi,c)
  theta <- theta(nA+nB,r,k,xi,c)
  theta * dR_Rdt(nA,r,k,xi=0,c) * nA - (1-theta) * dR_Rdt(nB,r,k,xi,c=1) * nB
}

# gain from protected - loss of increase to protected
# = A production * allocation determined by theta
# = mean total growth rate of pop - mean growth rate of B portion

# dR_dt_plus <- function(nB,r,k,xi,c){
#   dR_Rdt(nB,r,k,xi,c)*nB + irate(nB,r,k,xi,c)
# }

dR_dt_plus <- function(nB,r,k,xi,c,xi2=0){
  nA <- nAstar(nB,r,k,xi,c)
  theta <- theta(nA+nB,r,k,xi,c)
   ( dR_Rdt(nA,r,k,xi=0,c)*nA + (dR_Rdt(nB,r,k,xi,c=1)-xi2)*nB )
}

# Optimal habitat selection -----------------------------------------------

curve(theta(R=x,r=5,k=100,xi=1,c=1),xlim=c(0,100),ylim=c(0,1))
curve(theta(R=x,r=5,k=100,xi=2,c=2),add=T,lty=2)
curve(theta(R=x,r=5,k=100,xi=4,c=4),add=T,lty=3)
curve(theta(R=x,r=5,k=100,xi=10,c=10),add=T,lty=4)
curve(theta(R=x,r=5,k=100,xi=100,c=100),add=T,lty=4)
curve((x-10)/x,add=T,col="red")

curve(theta(R=x,r=5,k=100,xi=2,c=5),xlim=c(0,100),ylim=c(0,1))

Rstar(r=5,k=100,xi=1,c=1)
curve(inflow(r=5,k=100,xi=1*x,c=1*x),xlim=c(0,100))
thetastar(c=100)

curve(f(R=x,C=1,r=5,xi=1,k=100),xlim=c(0,100))
curve(f(R=x,C=2,r=5,xi=1,k=100),add=T,lty=2)
curve(f(R=x,C=5,r=5,xi=1,k=100),add=T,lty=3)

# Mean growth rate --------------------------------------------------------

curve(dR_Rdt_bar(R=10^x,r=5,k=100,xi=1,c=5),xlim=c(-2,2),ylim=c(1,8),n=10^3)
curve(dR_Rdt_bar(R=10^x,r=5,k=100,xi=2,c=1),add=T,lty=2,n=10^3,col="red")
curve(dR_Rdt_bar(R=10^x,r=5,k=100,xi=2,c=5),add=T,lty=3,n=10^3,col="orange")
curve(dR_Rdt_bar(R=10^x,r=5,k=100,xi=0,c=1),add=T,lty=4,n=10^3,col="blue")

par(mar=c(4,4,1,1))
curve(dR_Rdt_bar(R=10^x,r=5,k=100,xi=1,c=5),xlim=c(0,2),n=10^3,xlab="ln N",ylab=expression(r[t]))
abline(h=0,lty=3,col="gray")
abline(v=log10(Rstar(r=5,k=100,xi=1,c=5)),col="gray",lty=2)
curve(dR_Rdt(R=10^x,r=5,k=100,xi=0,c=5),add=T,col="blue",lty=2)
curve(dR_Rdt(R=10^x,r=5,k=100,xi=1,c=1),add=T,col="red",lty=3)
curve(dR_Rdt_bar2(R=10^x,r=5,k=100,xi=1,c=5),col="green",add=T)

curve(theta(R=10^x,r=5,k=100,xi=1,c=5),xlim=c(0,2),ylim=c(0,1))
curve(theta2(R=10^x,r=5,k=100,xi=1,c=5),add=T,col="green")
abline(v=log10(Rstar(r=5,k=100,xi=1,c=5)),col="gray",lty=2)

curve(nAstar(nB=x,r=5,k=100,xi=1,c=1),xlim=c(0,100),ylim=c(0,100),n=10^3)

# Immigration rates -------------------------------------------------------

### varying intercept (r and k in the same way)

curve(irate(nB=x,r=5,k=100,xi=1,c=1),xlim=c(0,125),n=10^3,
      xlab=expression(N[B]),ylab="i")
curve(irate(nB=x,r=0.8*5,k=0.8*100,xi=1,c=1),add=T,lty=2)
curve(irate(nB=x,r=0.6*5,k=0.6*100,xi=1,c=1),add=T,lty=3)
curve(irate(nB=x,r=5,k=100,xi=0,c=1),add=T,col="red")
  # B is independent of A when have identical growth functions
abline(h=0,lty=3,col="gray")

curve(irate(nB=x,r=5,k=100,xi=1,c=1)/x,xlim=c(2,125),n=10^3,
      xlab=expression(N[B]),ylab=expression(i/N[B]))
curve(irate(nB=x,r=0.8*5,k=0.8*100,xi=1,c=1)/x,add=T,lty=2)
curve(irate(nB=x,r=0.6*5,k=0.6*100,xi=1,c=1)/x,add=T,lty=3)
curve(irate(nB=x,r=5,k=100,xi=0,c=1)/x,add=T,col="red")
abline(h=0,lty=3,col="gray")

  # but this assumes that all individuals produced by A will move to B
  # how to siphon-off the individuals that stay in A, 
  # (as well as those that move from B->A)?

### varying k

curve(irate(nB=x,r=5,k=100,xi=1,c=1),xlim=c(0,125),n=10^3,
      ylab="absolute immigration")
curve(irate(nB=x,r=5,k=0.8*100,xi=1,c=1),add=T,lty=2)
curve(irate(nB=x,r=5,k=0.6*100,xi=1,c=1),add=T,lty=3)
  # NB: this is changing the k of *both* protected and exposed by the same fraction
abline(h=0,lty=3,col="gray")

### varying relative intercept (consumption)

curve(irate(nB=x,r=5,k=100,xi=2,c=1),xlim=c(0,125),n=10^3,
      ylab="absolute immigration")
curve(irate(nB=x,r=5,k=100,xi=1,c=1),add=T,lty=2)
curve(irate(nB=x,r=5,k=100,xi=0.5,c=1),add=T,lty=3)
  # higher consumption rate ->
  #   higher immigration at low B densities (because more sheltering in A)
  #   lower immigration at higher B densities (because more leave for A)
abline(h=0,lty=3,col="gray")

### varying relative k 

curve(irate(nB=x,r=5,k=100,xi=1,c=1),xlim=c(0,125),n=10^3,
      ylab="absolute immigration")
curve(irate(nB=x,r=5,k=100,xi=1,c=2),add=T,lty=2)
curve(irate(nB=x,r=5,k=100,xi=1,c=5),add=T,lty=3)
  # smaller k -> protected acts more like constant immigration
abline(h=0,lty=3,col="gray")

# Overall growth rate -----------------------------------------------------

curve(dR_dt_plus(nB=x,r=5,k=100,xi=0.5,c=1),xlim=c(0,125),n=10^3,
      xlab=expression(N[B]),ylab=expression(delta*N[B]/delta*t)
      )
curve(dR_dt_plus(nB=x,r=5,k=100,xi=1,c=1),n=10^3,add=T,lty=2)
curve(dR_dt_plus(nB=x,r=5,k=100,xi=2,c=1),n=10^3,add=T,lty=3)

curve(dR_Rdt(R=x,r=5,k=100,xi=0.5,c=1)*x,n=10^3,add=T,lty=1,col="red")
curve(dR_Rdt(R=x,r=5,k=100,xi=1,c=1)*x,n=10^3,add=T,lty=2,col="red")
curve(dR_Rdt(R=x,r=5,k=100,xi=2,c=1)*x,n=10^3,add=T,lty=3,col="red")

abline(h=0,lty=3,col="gray")
  # higher consumer pressure -> 
  #   higher immigration at low densities
  #   but lower growth at higher densities 
  #   (because of higher death rates and higher immmigration rates)
  # habitat selection -> stabilising:
  #   higher growth below K
  #   but lower growth above K

# Proportional refuge -----------------------------------------------------

curve(dR_dt_plus(nB=x,r=5,k=100,xi=0,c=1,xi2=1),xlim=c(0,125),n=10^3,
      xlab=expression(N[B]),ylab=expression(delta*N[B]/delta*t)
)
curve(dR_Rdt(x,r=5,k=100,xi=1,c=1)*x,add=T,col="red")


# Trial growth functions --------------------------------------------------

dR_dt_plus2 <- function(nB,r,k,xi,c){
  nA <- nAstar(nB,r,k,xi,c)
  theta <- theta(nA+nB,r,k,xi,c) # theta calculated using NEW variables
  theta * ( nA*(r*(1-c*nA/k)) +  nB*(r*(1-nB/k) - xi) )
}

dR_dt_plus3 <- function(nB,r,k,xi,c){
  ((-k * xi + k * r - nB * r) * (2* c *nB *r + k *xi + 2 *nB *r))/(2 *(c + 1)* k *r)
}

curve(dR_dt_plus2(nB=x,r=5,k=100,xi=0.5,c=1),xlim=c(0,125),n=10^3,
      xlab=expression(N[B]),ylab=expression(delta*N[B]/delta*t)
)
curve(dR_dt_plus2(nB=x,r=5,k=100,xi=1,c=1),n=10^3,add=T,lty=2)
curve(dR_dt_plus2(nB=x,r=5,k=100,xi=2,c=1),n=10^3,add=T,lty=3)

curve(dR_dt_plus3(nB=x,r=5,k=100,xi=0.5,c=1),xlim=c(0,125),n=10^3,
      xlab=expression(N[B]),ylab=expression(delta*N[B]/delta*t)
)
curve(dR_dt_plus3(nB=x,r=5,k=100,xi=1,c=1),n=10^3,add=T,lty=2)
curve(dR_dt_plus3(nB=x,r=5,k=100,xi=2,c=1),n=10^3,add=T,lty=3)

dR_dt_plus4 <- function(nB,r,k,xi,c){
 1/c*((1+c)*dR_Rdt(nB,r,k,xi,c=1) + xi*k/nB*(1-xi/r))
}
  
curve(dR_dt_plus4(nB=x,r=5,k=100,xi=0.5,c=1),n=10^3,add=T,col="blue")

dR_Rdt_bar2 <- function(nB,r,k,xi,c){
  r - 1/k*( nB+xi*k/r*( 1 - c*nB/( (c+1)*nB+xi*k/r) ) )
}

curve(dR_Rdt_bar(R=10^x,r=5,k=100,xi=1,c=5),xlim=c(-2,2),ylim=c(1,8),n=10^3)
curve(dR_Rdt_bar2(nB=10^x,r=5,k=100,xi=1,c=5),xlim=c(-2,2),ylim=c(1,8),n=10^3)


# Dynamic changes with consumer density -----------------------------------

# Is this similar to consumer interference?
#   Yes: lower feeding rates at higher consumer densities
#   No: affects DD of resource population too
