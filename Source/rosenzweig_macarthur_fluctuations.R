##################################################################################
# Simulated Rosenzweig-MacArthur dynamics with discrete temperature fluctuations # 
##################################################################################

# TODO
# Convert to single-species (lagged) model

romac <- function(N,P,r0,k,c,d,chi,mu,tau=0.01){
  dN_Ndt <- r0*(1-N/k) - c*P/(d+N)
  dP_Pdt <- chi*(c*N/(d+N) - mu)
  Nnew <- N + N*dN_Ndt*tau
  Pnew <- P + P*dP_Pdt*tau
  return(c(N=Nnew,P=Pnew))
}
  # Turchin 2003 p. 98
  # N = prey
  # P = predator
  # r0 = prey intrinsic rate of increase (low density, absence of predator)
  # k = prey carrying capacity (absence of predators)
  # c = rate at which predator encounters prey
  # d = prey density at which predator is feeding at 1/2 max capacity
  # chi = number of predators produced for each prey eaten
  # mu = consumption rate required for predator to sustain itself
  # tau = step length

r0 <- 10
k <- 100
c <- 1000
d <- 1000
chi <- 0.5
mu <- 1

nt <- 10^4
N0 <- 100
P0 <- 100
Nmat <- matrix(nr=nt,nc=2)
Nmat[1,] <-  c(N0,P0)
  
for(i in 2:nt){
  Nmat[i,] <- romac(Nmat[i-1,1],Nmat[i-1,2],r0,k,c,d,chi,mu,tau=0.01)
}

matplot(log(Nmat),type="l")






