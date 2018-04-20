model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    tlag <- t - Delta
    if (tlag < 0) lags <- c(0,0)   # no initial predation
    else lags <- lagvalue(tlag)    # returns lags of R and N
    dR <- r*R*(1 - R/K) - a*R*N
    dN <- a*lags[1]*lags[2] - d*N    
    return(list(c(dR, dN)))  
  }) 
}

p <- c(r=1,K=1,a=1,c=1,d=0.5,Delta=10)
s <- c(R=1,N=0.1)
run(delay=TRUE)

