##################################################################################
### Exploratory time-series plots of predator-prey dynamics in continuous time ###
##################################################################################

alpha <- 10 # prey growth
beta <- 0.1 # encounter rate
gamma <- 0.1  # prey to predator conversion
delta <- 0.1  # death of predator
  
dx <- function(x,y,alpha,beta,abiotic=F){
  if(abiotic==F) return(alpha*x - beta*x*y)
  if(abiotic==T) return(alpha - beta*x*y)
  }
dy <- function(x,y,beta,gamma,delta){
  beta*gamma*x*y - delta*y
  }

nt <- 10^5
x <- y <- vector("numeric",length=nt)
x0 <- 100
y0 <- 10

x[1] <- x0
y[1] <- y0
steplen <- 10^-3
#i0 <- 0.1
#iota <- i0*rep(rep(c(0,x0),c(90,1)),length.out=nt)
#iota <- i0*rep(1,nt)

for(t in 2:nt){
  x[t] <- x[t-1] + dx(x[t-1],y[t-1],alpha=1,beta=0.1,abiotic=F)*steplen # + iota[t-1]
  y[t] <- y[t-1] + dy(x[t-1],y[t-1],beta=0.1,gamma,delta)*steplen
  }

matplot(1:nt,log(cbind(x,y)),type="l")

#i0 <- 0.1
#iota <- i0*rep(rep(c(0,x0),c(90,1)),length.out=nt)
#iota <- i0*rep(1,nt)

###
# a <- 2
# b <- -1
# abiotic <- function(x) a/x-b
# curve(abiotic(x),xlim=c(0,10))
# curve(x+x*abiotic(x),xlim=c(0,10))
# abline(0,1,col="red",lty=2)
