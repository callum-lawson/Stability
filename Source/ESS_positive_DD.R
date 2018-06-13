
pgrow <- function(N,C,v=10,k,a=0.1,b=0.1,h=10,mu=0.5){
  v*k/(v+a*N) - b*C/(1+a*h*N) - mu
}

hi5 <- function(p,N,C1,C2,v=10,k1=1.5,k2=1,a=0.1,b=0.1,h=10,mu1=0.5,mu2=0.5){
  v*k1/(v+a*p*N) - b*C1/(1+a*h*p*N) - v*k2/(v+a*(1-p)*N) + b*C2/(1+a*h*(1-p)*N) - p*(mu1-mu2)
}

pstarsearch <- Vectorize(
  function(N,C1,C2){
    require(rootSolve)
    try1 <- try(
      root1 <- uniroot.all(hi5, interval=c(0,1), N=N, C1=C1, C2=C2)
    )
    if(class(try1)=="try-error" | length(try1)==0){
      return(NA)
    }
    else{
      return(root1) # returns p
    } 
  },
  vectorize.args=c("N")
)

pstarpos <- function(N,C1,C2){
  ifelse(pgrow(N=N,k=1.5,C=C1)>=0 & pgrow(N=N,k=1,C=C2)>=0,
    pstarsearch(N=N,C1=C1,C2=C2),
    NA
  )
}

xseq <- 0:100

plot(pgrow(N=xseq,k=1.5,C=0),ylim=c(-0.5,1.5),type="l")
lines(pgrow(N=xseq,k=1,C=0),lty=2)
lines(pgrow(N=xseq,k=1.5,C=20),lty=3)

pstar <- pstarpos(N=xseq,C1=20,C2=0)
pstar1 <- sapply(pstar, function(x) x[1])
pstar2 <- sapply(pstar, function(x) if(length(x)==1) x[1] else x[2])

plot(pstar1~xseq,ylim=c(0,1),type="l")
lines(pstar2~xseq,lty=2)
  # phase transition back to p=0 at low N values

plot(pgrow(N=xseq,k=1.5,C=20),ylim=c(-0.5,1.5),type="l")
lines(pgrow(N=xseq,k=1,C=0),lty=2)
abline(v=c(100*pstar1[100],100*(1-pstar1[100])),lty=3,col="red")
abline(v=c(100*pstar2[100],100*(1-pstar2[100])),lty=3,col="blue")

hi6 <- function(N,v=10,k1=1.5,k2=1,a=0.1){
  (a*k1*N + k1*v - k2*v) / (a*k1*N + a*k2*N)
}
lines(hi6(xseq),xlim=c(0,100),type="l",lty=3)
# http://www.wolframalpha.com/input/?i=solve+v+k+%2F+(v+%2B+a+p+N)+-+v+l+%2F+(v+%2B+a+(1-p)+N)+for+p

# Allee effects in both ---------------------------------------------------

plot(pgrow(N=xseq,k=1.5,C=20),ylim=c(-0.5,1.5),type="l")
lines(pgrow(N=xseq,k=1,C=10),lty=2)

qstar <- pstarpos(N=xseq,C1=20,C2=10)
qstar1 <- sapply(qstar, function(x) x[1])
qstar2 <- sapply(qstar, function(x) if(length(x)==1) x[1] else x[2])

plot(qstar1~xseq,ylim=c(0,1),type="l")
lines(qstar2~xseq,lty=2)
  # still double ESS but pop can't sustain at point of flip

