### Effects of refuges on attack rates and functional response ###

# Prey always find refuges ------------------------------------------------

af <- function(x,a,p){
  ifelse(x<=p,0,a*(x-p)/x)
}
curve(af(x,2,0.01),xlim=c(0,10),xlab="resource density",ylab="attack rate")
curve(af(x,2,0.1),xlim=c(0,10),add=T,col="red")
curve(af(x,2,1),xlim=c(0,10),add=T,col="orange")
curve(af(x,2,5),xlim=c(0,10),add=T,col="blue")

# Prey find refuges with probability p ------------------------------------
# (or spend a fraction of their time p in them)

ff <- function(x,a,h,p){
  ifelse(x<=p,0,a*(x-p)/(1+a*h*(x-p)) )
}
curve(ff(x,2,1,1),xlim=c(0,10))
curve(af(x,2,1),add=T,col="red")
curve(af(x,1,2),add=T,col="blue")

# Hill function -----------------------------------------------------------

real <- function(x,a=2,h=2,theta=1.5){
  a*x^theta
}
  # from Real / Hill functional response 
curve(real(x,theta=1.5),xlim=c(0,10))
curve(real(x,theta=3),xlim=c(0,10))

# Turchin models ----------------------------------------------------------

require(deSolve)
tseq <- 0:100
dAB <- function(t,y,parms){
  A <- y[1]
  B <- y[2]
  with(parms, {
    dA <- b*B - a*A
    dB <- a*A - b*B
    list(c(dA=dA,dB=dB))
  })
}
AB <- ode(y=c(2,1),times=tseq,func=dAB,parms=list(a=1,b=2))

dAB <- function(t,y,parms){
  A <- y[1]
  B <- y[2]
  with(parms, {
    dA <- (c*B + r*A)*(1-A/m)
    dB <- s*A - c*B*(1-A/m) - d*B
    list(c(dA=dA,dB=dB))
  })
}

s <- 1
c <- 1
m <- 1
d <- 1
r <- 1
AB <- ode(y=c(0.1,1),times=tseq,func=dAB,parms=list(s=s,c=c,m=m,d=d,r=r))
matplot(tseq,log(AB[,2:3]),type="l")
legend("bottomright",legend=LETTERS[1:2],col=1:2,lty=1:2,bty="n")

# Migration models --------------------------------------------------------

require(deSolve)
tseq <- 0:100

dAB <- function(t,y,parms){
  A <- y[1]
  B <- y[2]
  with(parms, {
    dA <- b*(A+B) - d*(A+B)^2 + m2*B*(1-A/k2) - m1*A*(1-B/k1)
    dB <- b*(A+B) - d*(A+B)^2 + m1*A*(1-B/k1) - m2*B*(1-A/k2) - x*B
    list(c(dA=dA,dB=dB))
  })
}

mypars <- list(
  b = 10,
  d = 1,
  m1 = 1000,
  m2 = 1000,
  k1 = 0.001,
  k2 = 0.001,
  x = 10
)

AB <- ode(y=c(1,1),times=tseq,func=dAB,parms=mypars)
matplot(tseq,log(AB[,2:3]),type="l")
legend("bottomright",legend=LETTERS[1:2],col=1:2,lty=1:2,bty="n")

# Abrams models -----------------------------------------------------------

require(deSolve)
tseq <- seq(0,10,length.out=1000)

dAB <- function(t,y,parms){
  A <- y[1]
  B <- y[2]
  with(parms, {
    dA <- A * (b - (A+B)/k - m1) + m2*B
    dB <- B * (b*0.01 - (A+B)/k - x - m2) + m1*A
    list(c(dA=dA,dB=dB))
  })
}

mypars <- list(
  b = 1000,
  k = 10,
  m1 = 10^8,
  m2 = 10^7,
  x = 100
)

par(mfrow=c(1,2))
AB <- ode(y=c(0.01,0.01),times=tseq,func=dAB,parms=mypars)
matplot(tseq,log(AB[,2:3]),type="l")
legend("bottomright",legend=LETTERS[1:2],col=1:2,lty=1:2,bty="n")
plot(log(AB[,2]/AB[,3])~tseq,type="l")

tail(AB)

# Refuge functions --------------------------------------------------------

par(mfrow=c(1,1))
trial <- function(R,k1,k2){
  plogis(k1 - k2 + R*0)
}
k1 <- 1
k2 <- 2
curve(trial(x,k1,k2),xlim=c(0,k1+k2))
abline(h=log(k1/k2),col="red",lty=2)

# Resource selection functions --------------------------------------------

curve(plogis(1 + 1*x)/(plogis(1 + 1*x)+plogis(1 + 1*x)),xlim=c(0,5))


# Predation isodars -------------------------------------------------------

mu <- 1 # mortality from consumers
Nseq <- 10^(seq(-2,2,length.out=100))
r <- 10
k <- 50
logistic <- function(N,r,k,mu){
  r * (1-N/k) - mu
}
par(mfrow=c(1,1))
plot(logistic(Nseq,r,k,mu=0)~Nseq,type="l")
lines(logistic(Nseq,r,k,mu=mu)~Nseq,lty=2)
abline(h=0,lty=3)
abline(v=k,lty=3)

N1fun <- function(N2,r1,r2,k1,k2,x){
  k1*(x/r1 - (r2/r1)*(1-N2/k2) + 1)
}

curve(
  N1fun(x,r1=5,r2=10,k1=50,k2=100,x=1)/(x+N1fun(x,r1=5,r2=10,k1=50,k2=100,x=1)),
                                        xlim=c(0,100),ylim=c(0,1)
                                        )
curve(N1fun(x,r1=5,r2=10,k1=50,k2=100,x=10)/x,add=T,lty=2)
abline(h=1,lty=3)

# oddfun <- function(N2,r1,r2,k1,k2,x){
#   plogis( log( 1 / ( ((r2-r1)*k2)/(r2*N2) + r1/r2 * k2/k1 ) ) )
# }
oddfun <- function(N2,r1,r2,k1,k2,x){
  plogis( log( k1/r1 * ( (r1 - r2 + x)/N2 + r2/k2 ) ) )
}
curve(oddfun(N2=x,r1=4,r2=5,k1=80,k2=100,x=0),xlim=c(0,100),ylim=c(0,1))
curve(oddfun(N2=x,r1=4,r2=5,k1=80,k2=100,x=0.5),add=T,lty=2)
curve(oddfun(N2=x,r1=4,r2=5,k1=80,k2=100,x=1),add=T,lty=3)
curve(oddfun(N2=x,r1=4,r2=5,k1=80,k2=100,x=2),add=T,lty=4)
curve(oddfun(N2=x,r1=4,r2=5,k1=80,k2=100,x=5),add=T,lty=5)
curve(oddfun(N2=x,r1=5,r2=5,k1=100,k2=100,x=0),add=T,lty=5)

curve(oddfun(N2=25,r1=5,r2=5,k1=100,k2=100,x=x),xlim=c(0,100),ylim=c(0,1))
curve(oddfun(N2=50,r1=5,r2=5,k1=100,k2=100,x=x),add=T,lty=2)
curve(oddfun(N2=75,r1=5,r2=5,k1=100,k2=100,x=x),add=T,lty=3)
curve(oddfun(N2=100,r1=5,r2=5,k1=100,k2=100,x=x),add=T,lty=4)

oddfun2 <- function(N2,r,k,x){
  plogis( log( 1 - (x*k)/(r*N2) ) )
}
curve(oddfun2(N2=x,r=5,k=100,x=1),xlim=c(0,100),ylim=c(0,1))
curve(oddfun2(N2=x,r=5,k=100,x=-1),add=T,lty=2)

N2isofun <- function(N1,r,k,x){
  N1 + x*k/r
}
par(mar=c(4,4,4,4))
curve(N2isofun(N1=x,r=5,k=100,x=1),xlim=c(0,100),ylim=c(0,100))
curve(N2isofun(N1=x,r=5,k=100,x=-1),add=T,lty=2)


# ESS ---------------------------------------------------------------------

ess <- function(N,r,k,x){
  p <- 1/2 - k*x / (4*r*N)
  p[p<=0] <- 0
  p[p>=1] <- 1
  return(p)
  # 1/2 - k*x / (4*r*N)
}
curve(ess(N=x,r=5,k=100,x=0.5),xlim=c(0,200),ylim=c(0,1))
curve(ess(N=x,r=5,k=100,x=1),add=T,lty=2)
curve(ess(N=x,r=5,k=100,x=2),add=T,lty=3)
curve(ess(N=x,r=5,k=100,x=5),add=T,lty=5)

curve(x*ess(N=x,r=5,k=100,x=1)/(1+1*x*ess(N=x,r=5,k=100,x=1)),
      xlim=c(0,100),ylim=c(0,1))

ess2 <- function(N,r,k,x,c){
  p <- (2*r*N - k*x) / (2*r*N*(1+c)) 
  p[p<=0] <- 0
  p[p>=1] <- 1
  return(p)
  # 1/2 - k*x / (4*r*N)
}
curve(ess2(N=x,r=5,k=100,x=2,c=0.25),xlim=c(0,100),ylim=c(0,1))
curve(ess2(N=x,r=5,k=100,x=2,c=0.5),add=T,lty=2)
curve(ess2(N=x,r=5,k=100,x=2,c=1),add=T,lty=3)
curve(ess2(N=x,r=5,k=100,x=2,c=2),add=T,lty=4)

curve(ess2(N=x,r=5,k=100,x=1,c=0.25),add=T,lty=1,col="red")
curve(ess2(N=x,r=5,k=100,x=1,c=0.5),add=T,lty=2,col="red")
curve(ess2(N=x,r=5,k=100,x=1,c=1),add=T,lty=3,col="red")
curve(ess2(N=x,r=5,k=100,x=1,c=2),add=T,lty=4,col="red")

curve(ess2(N=x,r=5,k=100,x=0,c=2),add=T,lty=4)

# Mixed DD functions ------------------------------------------------------

dN <- function(N,r,k,x,c){
  p <- (2*r*N - k*x) / (2*r*N*(1+c)) 
  p[p<=0] <- 0
  p[p>=1] <- 1
  dN <- p*r*(1-p*c*N/k) + (1-p)*(r*(1-(1-p)*N/k)+x)
  return(dN)
  # 1/2 - k*x / (4*r*N)
}

curve(dN(N=x,r=5,k=100,x=2,c=0.25),xlim=c(0,1000))
curve(dN(N=x,r=5,k=100,x=2,c=0.5),add=T,lty=2)
curve(dN(N=x,r=5,k=100,x=2,c=1),add=T,lty=3)
curve(dN(N=x,r=5,k=100,x=2,c=2),add=T,lty=4)
abline(h=0,col="red")


