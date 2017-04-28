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
curve(lnbh(x,m0,m1)-2.3,col="red",add=T)
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
# changing c1 -> drop-off occurs *earlier* at lower K

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
  # -0.5 -> reduce y by a freaction of exp(0.5)=1.6
  # decrease in y reduces equilibrium N by *more* when max popsize (K) is lower

# Y is just scaling up and down, changing c2 does the same thing?

# How does G affect BH curve? ---------------------------------------------

m0 <- 0.1
m1 <- 0.01
curve(lnbhdor(x,m0,m1,g=0.75),col="blue",xlim=c(0,5),ylim=c(-0.5,5))
curve(lnbhdor(x,m0,m1,g=0.5),col="red",add=T)
curve(lnbhdor(x,m0,m1,g=0.25),col="green",add=T)
  # above a certain Y, get more pop with high G, but as Y decreases, better to have low G

m0 <- 0.1
m1 <- 0.01
curve(lnbhdor(exp(x),m0,m1,g=0.75),col="blue",xlim=c(0,1.5),ylim=c(-0.5,5))
curve(lnbhdor(exp(x),m0,m1,g=0.25),col="red",add=T)
m0 <- 0.2
m1 <- 0.01
curve(lnbhdor(exp(x),m0,m1,g=0.75),col="blue",lty=2,add=T)
curve(lnbhdor(exp(x),m0,m1,g=0.25),col="red",lty=2,add=T)
  # When G is lower, reducing K has bigger effect on drop-off threshold for Y
  # (but will still always reduce it, so doesn't explain positive correlations)

m0 <- 0.1
m1 <- 0.01
curve(lnbhdor(x,m0,m1,g=0.75)-lnbhdor(x,m0,m1,g=0.25),col="blue",xlim=c(0,10))
m0 <- 0.2
m1 <- 0.01
curve(lnbhdor(x,m0,m1,g=0.75)-lnbhdor(x,m0,m1,g=0.25),col="red",add=T)
m0 <- 0.3
m1 <- 0.01
curve(lnbhdor(x,m0,m1,g=0.75)-lnbhdor(x,m0,m1,g=0.25),col="green",add=T)
abline(h=0,lty=3)
  # m1 -> lower germination favoured for given Y 

m0 <- 0.1
m1 <- 0.01
curve(lnbhdor(x,m0,m1,g=0.75)-lnbhdor(x,m0,m1,g=0.25),col="blue",xlim=c(0,10))
m0 <- 0.1
m1 <- 0.015
curve(lnbhdor(x,m0,m1,g=0.75)-lnbhdor(x,m0,m1,g=0.25),col="red",add=T)
abline(h=0,lty=3)
  # m2 -> proportional effect of G on dK is the same

# How does DI mortality relate to K? --------------------------------------

Kcalc <- function(m0,m1){
  exp(-m0)/( (1-exp(-m0)) * (m1/m0) )
} 
  # assuming T=1

curve(Kcalc(exp(x),1),xlim=c(-5,5))
curve(Kcalc(exp(x),2),col="red",add=T)
  # higher m0 -> lower K

# Optimal G - constant environment ----------------------------------------

# maximise population size (ignoring extinction)

np <- 10
yarr <- m0arr <- m1arr <- garr <- array(dim=rep(np,4))
yarr[] <- exp(seq(0,3,length.out=np))
m0arr[] <- exp(seq(-3,0,length.out=np))
m1arr[] <- exp(seq(-3,0,length.out=np))
garr[] <- plogis(seq(-5,5,length.out=np))

m0arr <- aperm(m0arr,c(2,1,3,4))
m1arr <- aperm(m1arr,c(3,2,1,4))
garr <- aperm(garr,c(4,2,3,1))

lKarr <- lnbhdor(yarr,m0arr,m1arr,garr,T3=0.3)
plot(lKarr[np,np,np,]~garr[np,np,np,],type="l")
require(fields)
par(mfrow=c(2,2))
matplot(garr[1,1,1,],t(lKarr[,5,1,]),type="l",lty=1,col=tim.colors(np))
matplot(garr[1,1,1,],t(lKarr[,15,1,]),type="l",lty=1,col=tim.colors(np))
matplot(garr[1,1,1,],t(lKarr[,35,1,]),type="l",lty=1,col=tim.colors(np))
matplot(garr[1,1,1,],t(lKarr[,50,1,]),type="l",lty=1,col=tim.colors(np))
  # m1 has no impact on optimal phenotype (only affects value of K)
  # G=0 is never optimal
  # S0/Y ratio high -> low G is optimal; So/Y ratio low -> G=1 is optimal

# Optimal G - variable environment ----------------------------------------

sgbh <- function(n_in,y,m0,m1,G,T3=0.3){
  So <- exp(-m0)
  Sn <- exp(-m0*T3)
  n_in*((1-G)*So + G*y*Sn / (1 + (1-Sn)*(m1/m0)*n_in))
}

nb <- 100
nt <- 1000 + nb
ns <- 5 # n y sig
narr <- yarr2 <- array(dim=c(nt,dim(garr),ns))
ysig <- exp(seq(-1,1,length.out=np))

yarr2[] <- exp(rnorm(
  n=nt*np^4*ns,
  mean=rep(rep(log(yarr),each=nt),times=ns),
  sd=rep(ysig,each=nt*np^4)
  ))

m0arr2 <- m1arr2 <- garr2 <- array(dim=c(dim(garr),ns))
m0arr2[] <- m0arr
m1arr2[] <- m1arr
garr2[] <- garr

N0 <- 1000
narr[1,,,,,] <- N0

for(t in 2:nt){
  narr[t,,,,,] <- sgbh(narr[t-1,,,,,],yarr2[t,,,,,],m0arr2,m1arr2,garr2)
}
  
lnarr <- log(narr)
muarr <- apply(lnarr[-(1:nb),,,,,],2:6,mean)
sdarr <- apply(lnarr[-(1:nb),,,,,],2:6,sd)

par(mfrow=c(2,2))
matplot(garr2[1,1,1,,1],muarr[5,2,5,,],type="l",lty=1,col=tim.colors(ns),
  ylim=c(2,5))
matplot(garr2[1,1,1,,1],muarr[5,4,5,,],type="l",lty=1,col=tim.colors(ns),
  ylim=c(2,5))
matplot(garr2[1,1,1,,1],muarr[5,6,5,,],type="l",lty=1,col=tim.colors(ns),
  ylim=c(2,5))
matplot(garr2[1,1,1,,1],muarr[5,8,5,,],type="l",lty=1,col=tim.colors(ns),
  ylim=c(2,5))
  # fluctuations favour G<1

ylim <- c(2,5)
par(mfrow=c(2,2))
matplot(garr2[1,1,1,,1],t(muarr[5,,5,,1]),type="l",lty=1,col=tim.colors(np),
  ylim=ylim)
matplot(garr2[1,1,1,,1],t(muarr[5,,5,,2]),type="l",lty=1,col=tim.colors(np),
  ylim=ylim)
matplot(garr2[1,1,1,,1],t(muarr[5,,5,,3]),type="l",lty=1,col=tim.colors(np),
  ylim=ylim)
matplot(garr2[1,1,1,,1],t(muarr[5,,5,,4]),type="l",lty=1,col=tim.colors(np),
  ylim=ylim)

ylim <- c(2,5)
par(mfrow=c(2,2))
matplot(garr2[1,1,1,,1],t(muarr[5,5,,,1]),type="l",lty=1,col=tim.colors(np),
  ylim=ylim)
matplot(garr2[1,1,1,,1],t(muarr[5,5,,,2]),type="l",lty=1,col=tim.colors(np),
  ylim=ylim)
matplot(garr2[1,1,1,,1],t(muarr[5,5,,,3]),type="l",lty=1,col=tim.colors(np),
  ylim=ylim)
matplot(garr2[1,1,1,,1],t(muarr[5,5,,,4]),type="l",lty=1,col=tim.colors(np),
  ylim=ylim)
  # m2 still shifting intercept but not optimal G?

# Optimal plastic G - variable environment --------------------------------







