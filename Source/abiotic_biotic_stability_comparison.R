### Comparing self-regulation around k between biotic and abiotic models ###

# NB: code not yet correct
# predation (increasing aC*) is vertical shift until new k is reached = horizontal perturbation

xn_f <- function(x,parms){
  with(parms, {
    kn <- mu / (beta * a)
    kr <- k / kn
    xn <- x + log(kr) # horizontal shift
    return(xn)
  })
}

abiotic <- function(x,parms){
  with(parms, {
    xn <- xn_f(x,parms)
    dxdt <- v * (-1 + k * exp(-xn))
    return(dxdt)
  })
}

biotic <-  function(x,parms){
  with(parms, {
    xn <- xn_f(x,parms)
    dxdt <- v * (1 - 1/k * exp(xn))
    return(dxdt)
  })
}

parms <- list(
  v = 1,
  k = 1,
  mu = 0.05,
  beta = 1,
  a = 1
)

abcurve <- function(parms,scale){
  lkn <- log(mu / (beta * a))  
  curve(abiotic(x,parms),xlim=c(lkn - scale, lkn + scale))
  curve(biotic(x,parms),col="red",add=T)
  abline(h=0,col="gray",lty=3)
}

par(mfcol=c(1,2),mar=c(3,3,3,3))
abcurve(parms,scale=1)
abcurve(parms,scale=0.001)
  # Currently incorrect
