### Investigate consumer-resource dynamics with multiple prey species using method ###
### from Geritz & Kisdi 2004 J Theor Biol                                          ###

# TODO: sort out uniroot function so that doesn't crash when no root

require(deSolve)

# Without consumer mortality ----------------------------------------------

f <- function(R,C,p=0.5,a=1,h=0){
  p*a*R*C / (1 + p*a*h*R)
}

g <- function(R,u,v){
  u - v*R
}

# g1 <- function(R,u1,K1){
#   u1*(1 - R/K1)
# }
#  # v = u/K

# g2 <- function(R,r=1,K=10){
#   R*r*(1 - R/K)
# }

dR <- function(R,C,p,a,h,u,v){
  g(R,u,v) - f(R,C,p,a,h)
}

kisdi_int <- function(nt,R1,R2,C,p,a,alpha=0.85,delta=0.1){
  Ei <- alpha * ( f(R1,C,p,a) +  f(R2,C,(1-p),a) )  # inflow rate of eggs
  Eo <- -delta                # outflow rate of eggs
  if(is.finite(nt)){
    return( 1/Eo * ( exp(Eo*(nt+log(Ei)/Eo)) - Ei ) )
  }
  if(!is.finite(nt)){
    return( -Ei/Eo )
  }
  # integral of eggs over time interval t
  # in this simple case, can be calculated explicitly at any given t
  # (i.e. no need for numerical integration)
  # (calculated with help from Shaopeng)
  # *only works for h=0*
}

kisdi_ode <- function(nt,R1,R2,C,p,a,h,alpha=0.85,delta=0.1){
  kidsi_cont <- function(t,y,parms=NULL){
    E <- y[1]
    dE <- alpha * ( f(R1,C,p,a,h) +  f(R2,C,(1-p),a,h) ) - delta * E  
    list(dE)
  }
  Et <- ode(y=c(E0=0),times=c(0,nt),func=kidsi_cont,parms=NULL)
  return(Et[2,2])
    # [2,2] because have to include t=0 in integration
}

### Breakdown

# curve(kisdi_int(t=x,R=1,C=1),xlim=c(0,100))
# curve(kisdi_int(t=10,R=x,C=1),xlim=c(0,10))
  # absolute growth rate is 
  # - saturating function of t
  # - linear function of resource availability (but this is set by C)

### Population growth rates

rC <- function(Cs,nt,p=0.5,q=0.5,a=1,h=0,u1=1,u2=0.5,v1=5,v2=0.5){
  
  if(q>0){
  try1 <- try(
    root1 <- uniroot(dR,C=Cs,p=p,a=a,h=h,u=q*u1,v=v1,lower=10^-10,upper=10^10) 
    )
  if(class(try1)=="try-error") K1s <- NA
  else K1s <- root1$root
  }
  else K1s <- 0
  
  if(q<1){
    try2 <- try(
      root2 <- uniroot(dR,C=Cs,p=(1-p),a=a,h=h,u=(1-q)*u2,v=v2,lower=10^-10,upper=10^10)
      )
    if(class(try2)=="try-error") K2s <- NA
    else K2s <- root2$root
  }
  else K2s <- 0

    # resource Ks adjust instantly to current C density
  if(h==0) Es <- kisdi_int(nt=nt,R1=K1s,R2=K2s,C=Cs,p=p,a=a)
  if(h>0) Es <- kisdi_ode(nt=nt,R1=K1s,R2=K2s,C=Cs,p=p,a=a,h=h)
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

lCsseq <- seq(0,2,length.out=nCs)
Csseq <- exp(lCsseq)

qseq <- pseq!=0.5
rCmat <- array(dim=c(nCs,na,np))
for(i in 1:nCs){
  for(j in 1:na){
    for(k in 1:np){
      rCmat[i,j,k] <- rC(Cs=Csseq[i],nt=Inf,p=pseq[k],a=aseq[j],
                         u1=1,u2=0.5,q=pseq[k]) # nt=0.1
    }
  }
}
  
require(fields)
matplot(lCsseq,rCmat[,,1],type="l",col="red",lty=1)
matplot(lCsseq,rCmat[,,3],type="l",col="blue",lty=3,add=T)
matplot(lCsseq,rCmat[,,2],type="l",col="purple",lty=2,add=T)
abline(h=0,lty=3)

matplot(lCsseq,rCmat[,5,] - rCmat[,1,], type="l",col=c("red","purple","blue")) 

matplot(lCsseq, rCmat[,,2] - rCmat[,,1], col="black", type="l")
  # difference increasing -> generalist more stable
matplot(lCsseq, rCmat[,,2] - rCmat[,,3], col="black", type="l")
  # difference decreasing -> generalist less stable

### Synchronised resource fluctuations (DI)

lCsseq <- seq(1.5,3.5,length.out=nCs)
Csseq <- exp(lCsseq)

rCmat <- array(dim=c(nCs,na,np))
for(i in 1:nCs){
  for(j in 1:na){
    for(k in 1:np){
      rCmat[i,j,k] <- rC(Cs=Csseq[i],pseq[k],nt=Inf,
                         u1=1+aseq[j],u2=0.5+aseq[j],q=pseq[k])
    }
  }
}

require(fields)
matplot(lCsseq,rCmat[,,1],type="l",col="red",lty=1)
matplot(lCsseq,rCmat[,,3],type="l",col="blue",lty=3,add=T)
matplot(lCsseq,rCmat[,,2] ,type="l",col="purple",lty=2,add=T)
abline(h=0,lty=3)

matplot(lCsseq,rCmat[,5,] - rCmat[,1,], type="l",col=c("red","purple","blue"))  # density-independent 

matplot(lCsseq, rCmat[,,2] - rCmat[,,1], col="black", type="l")
  # difference increasing -> generalist more stable
matplot(lCsseq, rCmat[,,2] - rCmat[,,3], col="black", type="l")
  # difference decreasing -> generalist less stable

### Opposite resource fluctuations

lCsseq <- seq(1.5,3.5,length.out=nCs)
Csseq <- exp(lCsseq)

rCmat <- array(dim=c(nCs,na,np))
for(i in 1:nCs){
  for(j in 1:na){
    for(k in 1:np){
      rCmat[i,j,k] <- rC(Cs=Csseq[i],pseq[k],nt=Inf,
                         u1=1+aseq[j],u2=0.5+rev(aseq)[j],q=pseq[k])
    }
  }
}

require(fields)
matplot(lCsseq,rCmat[,,1],type="l",col="red",lty=1)
matplot(lCsseq,rCmat[,,3],type="l",col="blue",lty=3,add=T)
matplot(lCsseq,rCmat[,,2] ,type="l",col="purple",lty=2,add=T)
abline(h=0,lty=3)

matplot(lCsseq,rCmat[,,2] ,type="l",col="purple",lty=4)

matplot(lCsseq,rCmat[,5,] - rCmat[,1,], type="l",col=c("red","purple","blue"))  # density-independent 

matplot(lCsseq, rCmat[,,2] - rCmat[,,1], col="black", type="l")
# difference increasing -> generalist more stable
matplot(lCsseq, rCmat[,,2] - rCmat[,,3], col="black", type="l")
# difference decreasing (slightly) -> generalist less stable

# Handling time analyses --------------------------------------------------

### Attack rates

lCsseq <- seq(0,3,length.out=nCs)
Csseq <- exp(lCsseq)

rCmat <- array(dim=c(nCs,na,np))
for(i in 1:nCs){
  for(j in 1:na){
    for(k in 1:np){
      rCmat[i,j,k] <- rC(Cs=Csseq[i],nt=100,p=pseq[k],a=aseq[j],h=1,
                         u1=1,u2=0.5,v1=5,v2=0.5,q=pseq[k]
                         )
    }
  }
}

require(fields)
matplot(lCsseq,rCmat[,,1],type="l",col="red",lty=1)
matplot(lCsseq,rCmat[,,3],type="l",col="blue",lty=3,add=T)
matplot(lCsseq,rCmat[,,2],type="l",col="purple",lty=2,add=T)
abline(h=0,lty=3)

matplot(lCsseq,matplot(lCsseq,rCmat[,5,] - rCmat[,1,], type="l",col=c("red","purple","blue")) , type="l",col=c("red","purple","blue")) 
 # density-*dependent*

matplot(lCsseq, rCmat[,,2] - rCmat[,,1], col="black", type="l")
# difference increasing -> generalist more stable
matplot(lCsseq, rCmat[,,2] - rCmat[,,3], col="black", type="l")
# difference decreasing -> generalist less stable

### Resource growth

lCsseq <- seq(1.5,3.5,length.out=nCs)
Csseq <- exp(lCsseq)

rCmat <- array(dim=c(nCs,na,np))
for(i in 1:nCs){
  for(j in 1:na){
    for(k in 1:np){
      rCmat[i,j,k] <- rC(Cs=Csseq[i],nt=100,p=pseq[k],a=1,h=2,
                         u1=1+aseq[j],u2=0.5+aseq[j],v1=5,v2=0.5,q=pseq[k]
      )
    }
  }
}

require(fields)
matplot(lCsseq,rCmat[,,1],type="l",col="red",lty=1)
matplot(lCsseq,rCmat[,,3],type="l",col="blue",lty=3,add=T)
matplot(lCsseq,rCmat[,,2],type="l",col="purple",lty=2,add=T)
abline(h=0,lty=3)

matplot(lCsseq,matplot(lCsseq,rCmat[,5,] - rCmat[,1,],type="l",col=c("red","purple","blue"))) # slightly density-*dependent*

matplot(lCsseq, rCmat[,,3] - rCmat[,,1], col="black", type="l")
  # not *quite* density-independent
matplot(lCsseq, rCmat[,,2] - rCmat[,,1], col="black", type="l")
# difference increasing -> generalist more stable
matplot(lCsseq, rCmat[,,2] - rCmat[,,3], col="black", type="l")
# difference decreasing -> generalist less stable

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
