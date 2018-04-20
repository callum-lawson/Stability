model <- function(t, state, parms){
  with(as.list(c(state,parms)),{
    R <- state[1:n]
    N <- state[(n+1):(2*n)]
    S <- sum(R)
    dR <- b*R*(1-S) - d1*R - a*R*N
    dN <- a*R*N - d2*N
    return(list(c(dR,dN)))    
  }) 
} 

n <- 3                         # number of species
b <- rnorm(n,mean=1,sd=0.1)    # b is a global parameter
p <- c(d1=0.1,d2=0.2,a=1)      # other parameters
R <- rep(0.1/n,n)              # initial condition of R
names(R) <- paste("R",seq(1,n),sep="")
N <- rep(0.01/n,n)             # initial condition of N
names(N) <- paste("N",seq(1,n),sep="")
s <- c(R,N)                    # combine R and N into s
run()
