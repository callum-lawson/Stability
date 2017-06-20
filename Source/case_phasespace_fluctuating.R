#################################################################################
# Phase-space plots when underlying predator or prey paramaters are fluctuating # 
#################################################################################

# Turchin 2003 p. 98
# N = prey
# P = predator
# r0 = prey intrinsic rate of increase (low density, absence of predator)
# k = prey carrying capacity (absence of predators)
# c = rate at which predator encounters prey
# d = prey density at which predator is feeding at 1/2 max capacity
# chi = number of predators produced for each prey eaten
# mu = consumption rate required for predator to sustain itself

require(rootSolve)

r0 <- 10
k <- 100
# c <- 1000
d <- 1000
chi <- 0.5
mu <- 1

dN_dt_f <- function(N,P,c){
  N*( r0*(1-N/k) - c*P/(d+N) )
} 

dP_dt_f <- function(N,P,c){
  P*( chi*(c*N/(d+N) - mu) )
} 

Nstar_f <- function(P,c){
  e <- try( 
    d <- uniroot(dN_dt_f, P=P, c=c, lower=10^-100, upper=10^100), 
    silent = TRUE 
  ) 
  if(class(e)=="try-error") { 
    return(NA) 
  } 
  else{ 
    return(d$root)	
  } 
}

c1 <- 1000
c2 <- 2000

nP <- 100
Pmin <- 0
Pmax <- 10
Pseq <- seq(Pmin,Pmax,length.out=nP)

Nstar1 <- sapply(Pseq,Nstar_f,c=c1)
Nstar2 <- sapply(Pseq,Nstar_f,c=c2)

Ndash1 <- uniroot(dP_dt_f, P=1, c=c1, lower=10^-100, upper=10^100)$root
Ndash2 <- uniroot(dP_dt_f, P=1, c=c2, lower=10^-100, upper=10^100)$root 
  # P=arbitrary

Pstar1 <- uniroot(dN_dt_f, N=Ndash1, c=c1, lower=10^-100, upper=10^100)$root
Pstar2 <- uniroot(dN_dt_f, N=Ndash1, c=c2, lower=10^-100, upper=10^100)$root

plot(Pseq~Nstar1,type="l",col="red",xlab="N",ylab="P")
lines(Pseq~Nstar2,col="blue")
abline(v=Ndash1,col="red",lty=2)
abline(v=Ndash2,col="blue",lty=2)
abline(h=Pstar1,col="red",lty=3)
abline(h=Pstar2,col="blue",lty=3)

### When have predator DD (e.g. interference):

Nmin <- min(c(Nstar1,Nstar2),na.rm=T)
Nmax <- max(c(Nstar1,Nstar2),na.rm=T)
Nseq <- seq(Nmin,Nmax,length.out=nP)

Pstar_f <- function(N,c){
  e <- try( 
    d <- uniroot(dP_dt_f, N=N, c=c, lower=10^-100, upper=10^100), 
    silent = TRUE 
  ) 
  if(class(e)=="try-error") { 
    return(NA) 
  } 
  else{ 
    return(d$root)	
  } 
}
