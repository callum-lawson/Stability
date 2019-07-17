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

# Eigenvalues -------------------------------------------------------------

abiotic_f <- function(t,y,parms){
  with(parms, {
    X <- y[1]
    Y <- y[2]
    dX <- exp(-X) - 1 - R_0 * exp(Y)
    dY <- p * (R_0 * exp(X) - 1)
    list( c( dX=dX, dY=dY ) )
  })
}

trial <- ode(y=c(X=0,Y=0),times=0:100,func=abiotic_f,parms=list(R_0=2,p=0.5)) 
plot(trial)

TS_abiotic_sim_f <- function(R_0,p){
  
  parms_sim <- list(
    R_0 = R_0,
    p = p
  )

  equ <- steady(y = c(X=0,Y=0), times=c(0,Inf), func=abiotic_f, parms = parms_sim, method = "runsteady")$y
  jac <- jacobian.full(y=equ,fun=abiotic_f,parms=parms_sim,time=0)
  
  eig <- eigen(jac)
  vec <- eig$vectors
  val <- eig$values
  TS <- val[1] / val[2]
  
}

R_0_seq <- 10^seq(0.01,1,length.out=100)
p_seq <- seq(0.01,0.99,length.out=100)
TS_abiotic_sim_mat <- matrix(nr=100,nc=100)
TS_abiotic_sim_mat[] <- mapply(TS_abiotic_sim_f,R_0=R_0_seq,p=p_seq)

par(mfrow=c(1,1))
image(R_0_seq,p_seq,log(TS_abiotic_sim_mat))
