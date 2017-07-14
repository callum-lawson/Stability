##################################################################################
# Investigate consumer-resource dynamics with multiple prey species using method #
# from Geritz & Kisdi 2004 J Theor Biol                                          #
##################################################################################

# TODO: sort out uniroot function so that doesn't crash when no root

require(deSolve)

# Without consumer mortality ----------------------------------------------

f <- function(R,C,a=1){
  a*R*C
}

g1 <- function(R,u=0.1,K=10){
  u*(1 - R/K)
}

g2 <- function(R,r=1,K=10){
  R*r*(1 - R/K)
}

dR1 <- function(R,C){
  g1(R) - f(R,C)
}

dR2 <- function(R,C){
  g2(R) - f(R,C)
}

kisdi_int <- function(t,R,C,alpha=0.85,delta=0.1){
  Ei <- alpha * f(R,C)    # inflow rate of eggs
  Eo <- -delta            # outflow rate of eggs
  return( 1/Eo * (exp( Eo * (t + log(Ei)/Eo) ) - Ei) )
  # asymptote at -Eo/Ei
}

### Population growth rates

rC <- function(Cs,nt){
  K1s <- uniroot(dR1,C=Cs,lower=10^-10,upper=10^10)$root
  K2s <- uniroot(dR2,C=Cs,lower=10^-10,upper=10^10)$root
  Es <- kisdi_int(t=nt,R=K1s+K2s,C=Cs)
  rC <- log(Es/Cs)
  return(rC)
}

na <- 5
aseq <- exp(seq(-1,1,length.out=na))
nCs <- 100
lCsseq <- seq(-5,-1,length.out=nCs)
Csseq <- exp(lCsseq)
rCmat <- matrix(nr=nCs,nc=na)
for(i in 1:nCs){
  for(j in 1:na){
    f <- function(R,C,a=aseq[j]){
      a*R*C
    }
    rCmat[i,j] <- rC(Cs=Csseq[i],nt=0.1)
  }
}
  
matplot(lCsseq,rCmat,type="l")
abline(h=0,lty=3)

### Temporal dynamics

C0 <- 0.1
ns <- 20
sseq <- 1:ns
nt <- 0.1

Cs <- Es <- K1s <- K2s <- vector()
for(s in 1:ns){
  Cs[s] <- ifelse(s==1,C0,Es[s-1])
  K1s[s] <- uniroot(dR1,C=Cs[s],lower=10^-10,upper=10^10)$root
  K2s[s] <- uniroot(dR2,C=Cs[s],lower=10^-10,upper=10^10)$root
  Es[s] <- kisdi_int(t=nt,R=K1s[s]+K2s[s],C=Cs[s])
}

matplot(sseq,log(cbind(K1s,K2s,Cs)),type="l",lty=1)

# With consumer mortality -------------------------------------------------

kisdi <- function(y,R1,R2,a,mu,delta){
  
  R1 <- y[1]
  R2 <- y[2]
  E  <- y[3]
  C  <- y[4]
  
  dR1 <- R1*g1(R1) - a*R1*C
  dR2 <- R2*g2(R2) - a*R2*C
  dE  <- alpha*a*C*(R1+R2) - delta*E
  dC  <- -mu*C
  
  list(c(dR1,dR2,dE,dC))
  
}

# ode1 <- ode(y=c(R0=R0,C0=C0),times=tseq,func=romac_dis,parms=NULL)
