####################################################################################
# Include temperature fluctuations in pred-prey models from Fussmann et al. (2014) #
####################################################################################

arr <- function(Tt,E0,E1,T0=293.15,k=8.617*10^-5){
  E0 * exp( E1*(Tt-T0) / (k*Tt*T0) )
}
  # arrhenius equation

cr <- function(R,C,r,K,y,B,eps,x,tau=0.1){
  dR_Rdt <- r*(1-R/K) - y*C/(B+R)
  dC_Cdt <- eps*(y*R/(B+R) - x)
  Rnew <- R + R*dR_Rdt*tau
  Cnew <- C + C*dC_Cdt*tau
  return(c(R=Rnew,C=Cnew))
}
  # R = prey
  # C = predator
  # r = prey intrinsic rate of increase (low density, absence of predator)
  # K = prey carrying capacity (absence of predators)
  # c = rate at which predator encounters prey
  # B = prey density at which predator is feeding at 1/2 max capacity
  # eps = number of predators produced for each prey eaten
  # x = consumption rate required for predator to sustain itself
  # tau = step length

rE0 <- 8.715*10^-7
KE0 <- 5.623
yE0 <- 8.408*10^-6
BE0 <- 3.664
xE0 <- 2.689*10^-6
rE1 <- 0.84
KE1 <- -0.772
yE1 <- 0.467
BE1 <- -0.114
xE1 <- 0.639
  # change slopes later
  # r units are per SECOND; pop more than triples every 24h
  
nt <- 10^5 # time in SECONDS
Tmu <- 30
Tt <- rep(Tmu + 293.15,nt)

r <- arr(Tt,rE0,rE1)
K <- arr(Tt,KE0,KE1)
y <- arr(Tt,yE0,yE1)
B <- arr(Tt,BE0,BE1)
x <- arr(Tt,xE0,xE1)
  # will have to hold each of these constants for multiple timesteps
  # (otherwise numerical approximation of continuous time doesn't work)

eps <- 0.85

R0 <- 100
C0 <- 10
N <- matrix(nr=nt,nc=2)
N[1,] <-  c(R0,C0)

for(t in 2:nt){
  N[t,] <- cr(N[t-1,1],N[t-1,2],r[t],K[t],y[t],B[t],eps,x[t],tau=100)
}

matplot(log(N),type="l")






