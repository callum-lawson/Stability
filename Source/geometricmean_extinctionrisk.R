### Does a population's mean log size correlate with its extinction risk? ###

# Functions ---------------------------------------------------------------

rBH <- function(N,r0,k,eps=0){
  N * k / ( N + (k-N)*exp(-r0) ) * exp(eps) # arbitrary
}
  # https://en.wikipedia.org/wiki/Beverton%E2%80%93Holt_model

# Parameters --------------------------------------------------------------

nr <- 15
nk <- 15
ns <- 15
np <- nr*nk*ns

nt <- 1250
nb <- 250
ni <- 10
nT <- nt*ni*np

dd <- expand.grid(
  r0 = exp(seq(0,5,length.out=nr)),
  k = exp(seq(0,5,length.out=nk)),
  sig = exp(seq(-1,1,length.out=ns))
)

# Simulate dynamics -------------------------------------------------------

eps <- array(dim=c(nt,ni,np))
eps[] <- rnorm(nT,0,rep(dd$sig,each=nt*ni))

N0 <- 1
N <- array(dim=c(nt,ni,np))
N[1,,] <- N0

for(t in 1:(nt-1)){
  N[t+1,,] <- with(dd, rBH(N[t,,], r0, k, eps[t,,]))
}

# Analyse output ----------------------------------------------------------

gm <- apply(N[nb:nt,,],3,function(x) mean(log(x)))
tau <- 10^-0.5
ex <- apply(apply(N[nb:nt,,],c(2,3),function(x) TRUE %in% (x<tau)),2,mean)
exr <- qlogis(ifelse(ex==1,0.99,ifelse(ex==0,0.01,ex)))
mte <- apply(apply(N[nb:nt,,],c(2,3),function(x) which(x<tau)[1]),2,mean)

plot(exr~gm)
lines(supsmu(gm,exr),col="red")
summary(lm(exr~gm+I(gm^2)))
  # roughly proportional to extinction risk?

plot(log(mte)~gm)
lines(supsmu(gm,log(mte)),col="red")
summary(lm(log(mte)~gm+I(gm^2)))

