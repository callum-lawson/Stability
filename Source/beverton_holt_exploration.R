####################################################################
# Explore how random variation in Beverton-Holt mortality function #
# affects the distribution of population sizes                     #
####################################################################

# Functions ---------------------------------------------------------------

BH <- function(n_in,m0,m1){
  n_out <- n_in * exp(-m0) / ( 1 + (m1/m0)*(1-exp(-m0))*n_in )
  return(n_out)
}
# assuming T = 1

# Exploration -------------------------------------------------------------

ln0_min <- -5
ln0_max <- 0
nseq <- 10^3
neps <- 10^3
alpha_m <- 0.1
beta_m <- 1

# n0_seq <- exp(seq(ln0_min,ln0_max,length.out=nseq))
# eps_m <- seq(-1,1,length.out=neps)
n0_seq <- exp(rnorm(nseq,0,1.5))
eps_m <- rnorm(neps,0,2)
  
m0 <- m1 <- n0 <- n1 <- matrix(NA,nr=nseq,nc=neps)

m0[] <- rep(exp(outer(alpha_m,eps_m,"+")),each=nseq)
m1[] <- rep(exp(outer(beta_m,eps_m,"+")),each=nseq)
n0[] <- rep(n0_seq,times=neps)
n1 <- BH(n0,m0,m1)

par(mfrow=c(2,1))
hist(n0,breaks=100,prob=T)
hist(n1,breaks=100,prob=T)

hist(n1[1,],breaks=100)
hist(n1[2,],breaks=100)
hist(n1[3,],breaks=100)

matplot(n0[,1],n1,type="l")
abline(0,1,col="purple",lty=3)

(K <- BH(10^10,m0[1,],m1[1,]))

# Plant example -----------------------------------------------------------
# From Yodzis and work on desert annuals

### Without seed dormancy

curve(log(1-1/exp(x)),col="red",xlim=c(0,0.5))
curve(log(0.9-1/exp(x)),col="blue",add=T)
# drop-off occurs earlier at lower K

lnbh <- function(y,m0,m1,T3=0.3){
  c1 <- exp(-m0*T3)
  c2 <- (1-exp(-m0)*T3)*m1/m0
  log((y*c1-1)/(c2*y))
}

m0 <- 1.1
m1 <- 0.1
curve(lnbh(x,m0,m1),col="blue",xlim=c(1,3))
m1 <- 0.01
curve(lnbh(x,m0,m1),col="red",add=T)
# changing c2 only changes equilibrium N by constant factor (shifts intercept)

m0 <- 1
m1 <- 0.01
c1 <- exp(-m0)
c2 <- (1-exp(-m0))*m1/m0
curve(log((exp(x)*c1-1)/(c2*exp(x))),col="blue",xlim=c(0,3))
m0 <- 0.5
c1 <- exp(-m0)
c2 <- (1-exp(-m0))*m1/m0
curve(log((exp(x)*c1-1)/(c2*exp(x))),col="red",add=T)
# drop-off occurs earlier at *lower* K

# Y is just scaling up and down, changing c2 does the same thing?

### With seed dormancy

lnbhdor <- function(y,m0,m1,g,T3=0.3){
  s0 <- exp(-m0)
  c1 <- exp(-m0*T3)
  c2 <- (1-exp(-m0)*T3)*m1/m0
  log((c1*g*y-g*s0+s0-1) / (c2*g*y*((g-1)*s0+1)))
}

g <- 0.4
m0 <- 1
m1 <- 0.01
curve(lnbhdor(x,m0,m1,g),col="blue",xlim=c(1,3))
m0 <- 0.5
curve(lnbhdor(x,m0,m1,g),col="red",add=T)
  # Same result with G - drop-off earlier for lower K

# How does Y affect BH curve? ---------------------------------------------

BH <- function(n_in,m0,m1){
  n_out <- n_in * exp(-m0) / ( 1 + (m1/m0)*(1-exp(-m0))*n_in )
  return(n_out)
}


y <- 3 # log scale
m0 <- 1
curve(log(BH(exp(x+y),m0,m1)),col="blue",xlim=c(-1,5))
curve(log(BH(exp(x+y-0.5),m0,m1)),col="red",add=T)
m0 <- 2
curve(log(BH(exp(x+y),m0,m1)),col="blue",add=T,lty=2)
curve(log(BH(exp(x+y-0.5),m0,m1)),col="red",add=T,lty=2)
abline(0,1,lty=3)
  # -0.5 -> exp(0.5) decrease in y
  # decrease in y reduces equilibrium N *faster* when max popsize (K) is lower

