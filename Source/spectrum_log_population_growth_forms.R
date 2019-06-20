### Deriving different basic forms of density-dependence in population growth ###

r_abiotic <- function(X,v,k){
  N <- exp(X)
  K <- exp(k)
  1/N * v * (K - N)
}

r_gomp <- function(X,v,k){
  v * (k - X)
}

r_unknown <- function(X,v,k){
  N <- exp(X)
  1/N * v * (k - X)
}

r_logistic <- function(X,v,k){
  N <- exp(X)
  K <- exp(k)
  v * (K - N)
}

rflist <- list(r_abiotic, r_gomp, r_logistic, r_unknown)
Xseq <- seq(-1,1,length.out=100)
rmat <- sapply(rflist, mapply, Xseq, v=1, k=0)

par(mfrow=c(1,2),mar=c(4.5,4.5,1,1))
matplot(Xseq,rmat,type="l", xlab="ln N", ylab="dN/Ndt")
matplot(exp(Xseq),rmat,type="l", xlab="N", ylab="dN/Ndt")

# Temporal growth curves --------------------------------------------------

require(deSolve)
dlnNdt <- function(t,y,parms){
  with(parms,{
    X <- y[1]
    dX <- FUN(X,v,k)
    list(c(dX))
  })
}

tseq <- seq(0,10,length.out=100)
bparms <- c(v=1,k=0)
X0 <- c(X=-5)
parms_abiotic <- as.list(c(FUN=r_abiotic,bparms))
parms_logistic <- as.list(c(FUN=r_logistic,bparms))

Xt_abiotic <- ode(y=X0,times=tseq,func=dlnNdt,parms=parms_abiotic)[,"X"]
Xt_logistic <- ode(y=X0,times=tseq,func=dlnNdt,parms=parms_logistic)[,"X"]

matplot(tseq,cbind(Xt_abiotic,Xt_logistic),type="l")
matplot(tseq,exp(cbind(Xt_abiotic,Xt_logistic)),type="l")
