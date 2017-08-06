### Investigate consumer-resource dynamics with multiple prey species using method ###
### from Geritz & Kisdi 2004 J Theor Biol                                          ###

# TODO: sort out uniroot function so that doesn't crash when no root

require(deSolve)

# Without consumer mortality ----------------------------------------------

f <- function(R,C,a=1){
  a*R*C
}

g1 <- function(R,u1,K1){
  u1*(1 - R/K1)
}

g2 <- function(R,u2,K2){
  u2*(1 - R/K2)
}

# g2 <- function(R,r=1,K=10){
#   R*r*(1 - R/K)
# }

dR1 <- function(R,C,u1,K1){
  g1(R,u1,K1) - f(R,C)
}

dR2 <- function(R,C,u2,K2){
  g2(R,u2,K2) - f(R,C)
}

kisdi_int <- function(t,R,C,alpha=0.85,delta=0.1){
  Ei <- alpha * f(R,C)    # inflow rate of eggs
  Eo <- -delta            # outflow rate of eggs
  if(is.finite(t)){
    return( 1/Eo * ( exp(Eo*(t+log(Ei)/Eo)) - Ei ) )
  }
  if(!is.finite(t)){
    return( -Ei/Eo ) # asymptote
  }
  # integral of eggs over time interval t
  # in this simple case, can be calculated explicitly at any given t
  # (i.e. no need for numerical integration)
  # (calculated with help from Shaopeng - see )
}

### Breakdown

# curve(kisdi_int(t=x,R=1,C=1),xlim=c(0,100))
# curve(kisdi_int(t=10,R=x,C=1),xlim=c(0,10))
  # absolute growth rate is 
  # - saturating function of t
  # - linear function of resource availability (but this is set by C)

### Population growth rates

rC <- function(Cs,p,nt,u1=0.1,u2=0.05,K1=0.5,K2=1){
  
  try1 <- try(  root1 <- uniroot(dR1,C=Cs,u1=u1,K1=K1,lower=10^-10,upper=10^10) )
  if(class(try1)=="try-error") K1s <- NA
  else K1s <- root1$root
  
  try2 <- try(  root2 <- uniroot(dR2,C=Cs,u2=u2,K2=K2,lower=10^-10,upper=10^10) )
  if(class(try2)=="try-error") K2s <- NA
  else K2s <- root2$root

    # resource Ks adjust instantly to current C density
  Es <- kisdi_int(t=nt,R=p*K1s+(1-p)*K2s,C=Cs)
    # C eggs accumulate over interval
  rC <- log(Es/Cs)
    # Change in C = new eggs - death of all adults
  return(rC)
}

np <- 3
na <- 5
nCs <- 100

pseq <- seq(0,1,length.out=np)
aseq <- exp(seq(-1,1,length.out=na))

### Fluctuating attack rates

lCsseq <- seq(-2,0,length.out=nCs)
Csseq <- exp(lCsseq)

rCmat <- array(dim=c(nCs,na,np))
for(i in 1:nCs){
  for(j in 1:na){
    for(k in 1:np){
      f <- function(R,C,a=aseq[j]){
        a*R*C
      }
      rCmat[i,j,k] <- rC(Cs=Csseq[i],pseq[k],nt=Inf) # nt=0.1
    }

  }
}
  
require(fields)
matplot(lCsseq,rCmat[,,1],type="l",col="red",lty=1)
matplot(lCsseq,rCmat[,,3],type="l",col="blue",lty=2,add=T)
matplot(lCsseq,rCmat[,,2],type="l",col="purple",lty=4,add=T)
abline(h=0,lty=3)

### Resource fluctuations - DI

lCsseq <- seq(-0.5,1,length.out=nCs)
Csseq <- exp(lCsseq)

rCmat <- array(dim=c(nCs,na,np))
for(i in 1:nCs){
  for(j in 1:na){
    for(k in 1:np){
      f <- function(R,C,a=1){
        a*R*C
      }
      rCmat[i,j,k] <- rC(Cs=Csseq[i],pseq[k],nt=Inf,
                         u1=aseq[j]/10+0.1,u2=aseq[j]/10+0.05)
    }
  }
}

require(fields)
matplot(lCsseq,rCmat[,,1],type="l",col="red",lty=1)
matplot(lCsseq,rCmat[,,3],type="l",col="blue",lty=2,add=T)
matplot(lCsseq,rCmat[,,2] ,type="l",col="purple",lty=4,add=T)
abline(h=0,lty=3)

### Temporal dynamics

C0 <- 0.1
ns <- 20
sseq <- 1:ns
nt <- 0.1

Cs <- Es <- K1s <- K2s <- vector()
for(s in 1:ns){
  Cs[s] <- ifelse(s==1,C0,Es[s-1])
  
  try1 <- try(  root1 <- uniroot(dR1,C=Cs[s],lower=10^-10,upper=10^10)$root )
  if(class(try1)=="try-error") K1s[s] <- NA
  else K1s[s] <- root1 
  
  try2 <- try(  root2 <- uniroot(dR2,C=Cs[s],lower=10^-10,upper=10^10)$root )
  if(class(try2)=="try-error") K2s[s] <- NA
  else K2s[s] <- root2
  
  Es[s] <- kisdi_int(t=nt,R=K1s[s]+K2s[s],C=Cs[s])
}
  # still doesn't work

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
