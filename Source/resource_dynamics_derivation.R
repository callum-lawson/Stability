### Derive biotic resource dynamics using timescale separation ###

# Initial graphs ----------------------------------------------------------

curve(1/x-0.1,xlim=c(0.5,10)) # abiotic
curve(1/(1+x),add=T,lty=2) # BH
curve(1-x/10,add=T,lty=3) # logistic

curve(x*(1-x/10),xlim=c(0,10),lty=3) # logistic
curve(x*(1/x-0.1),add=T) # abiotic
curve(x*(1/(1+x)),add=T,lty=2) # BH

# Timescale separation ----------------------------------------------------

source("Source/predprey_functions_general.R")

### Joint R and C equilibrium

source("Source/grindR/grind.R")

# Set up params

x0 <- arrrate(0,e0["x"],e1["x"])

tmax <- 60^2 * 24 * 7 * 52 * 10 # maximum length of time in seconds
tf <- 1000
tseq <- seq(0,tmax,length.out=tf)

zmu <- 0 # =20°C 
zsig <- 5 # wave amplitude
zf <- 10 # wave frequency over whole time series
zl <- tmax/zf 

zparms <- list(zmu=zmu,zsig=zsig,zl=zl)
eparms <- list(e0=e0,e1=e1)

R0 <- 10^1
C0 <- 10^-2
y0 <- c(R=R0,C=C0)


dRCt_cont2 <- function(t,y,parms){
  with(parms, dRC_cont(y,m,r,k,a,h,x) )
}

# Equilibrium

fixparms <- list(m=10,r=0,k=100,a=0.1,h=0,x=0.1)

run(tmax=100, tstep=0.1, 
    state=y0, 
    parms=fixparms, 
    odes=dRCt_cont2, 
    ymin=0, ymax=NULL, log="", x=1, y=2, xlab="Time", ylab="Density", tmin=0, draw=lines, times=NULL, show=NULL, arrest=NULL, after=NULL, tweak=NULL, timeplot=TRUE, traject=FALSE, table=FALSE, add=FALSE, legend=TRUE, solution=FALSE, delay=FALSE, lwd=2, col="black", pch=20)

plane(xmin=0.1, xmax=1.1, ymin=0.1, ymax=1.1, log="", npixels=500, 
      state=y0, 
      parms=fixparms, 
      odes=dRCt_cont2, 
      tstep=1,
      x=1, y=2, time=0, grid=5, eps=0, show=NULL, portrait=TRUE, vector=FALSE, 
      add=FALSE, legend=TRUE, zero=TRUE, lwd=2, col="black", pch=20) 

RCstar <- newton(state=y0, 
       parms=fixparms, 
       odes=dRCt_cont2, 
       time=0, 
       x=1, 
       y=2, 
       positive=TRUE, 
       jacobian=FALSE, 
       vector=FALSE, 
       plot=TRUE
)

# method="runsteady"

q <- steady(y=y0, 
       parms=fixparms, 
       fun=dRCt_cont2, 
       time=c(0,Inf),
       method="runsteady"
       )

run(tmax=100, tstep=0.1, 
    state=y0, 
    parms=fixparms, 
    odes=dRCt_cont2, 
    ymin=0, ymax=NULL, log="", x=1, y=2, xlab="Time", ylab="Density", tmin=0, draw=lines, times=NULL, show=NULL, arrest=NULL, after=NULL, tweak=NULL, timeplot=TRUE, traject=FALSE, table=FALSE, add=FALSE, legend=TRUE, solution=FALSE, delay=FALSE, lwd=2, col="black", pch=20)

equ <- q$y
jac <- jacobian.full(y=equ,func=dRCt_cont2,parms=fixparms)
eig <- eigen(jac)
realeig <- Re(eig$values)
realeig[1]/realeig[2]
dom <- max(sort(Re(eig$values)))

continue(state=s, 
         parms=p, 
         odes=model, 
         x=1, step=0.01, xmin=0, xmax=1, y=2, ymin=0, ymax=1.1, log="", time=0, positive=FALSE, add=FALSE, ...)

# Light

fixparms <- list(m=10,r=0,k=Inf,a=0.1,h=0,x=0.1)

run(tmax=100, tstep=0.1, 
    state=y0, 
    parms=fixparms, 
    odes=dRCt_cont2, 
    ymin=0, ymax=NULL, log="", x=1, y=2, xlab="Time", ylab="Density", tmin=0, draw=lines, times=NULL, show=NULL, arrest=NULL, after=NULL, tweak=NULL, timeplot=TRUE, traject=FALSE, table=FALSE, add=FALSE, legend=TRUE, solution=FALSE, delay=FALSE, lwd=2, col="black", pch=20)

q <- steady(y=y0, 
            parms=fixparms, 
            fun=dRCt_cont2, 
            time=c(0,Inf),
            method="runsteady"
)

equ <- q$y
jac <- jacobian.full(y=equ,func=dRCt_cont2,parms=fixparms)
eig <- eigen(jac)
realeig <- Re(eig$values)
realeig[1]/realeig[2]
dom <- max(sort(Re(eig$values)))

# Feedback

dRC_rec <- function(y,a,h=0,x,alpha=0.85){
  R <- y[1]
  C <- y[2]
  Ft <- f(R,C,a,h)
  dR <- x*C/alpha - Ft
  dC <- alpha * Ft - x*C
  list(c(dR=dR,dC=dC))
}

dRCt_cont3 <- function(t,y,parms){
  with(parms, dRC_rec(y,a,h,x) )
}

run(tmax=100, tstep=0.1, 
    state=c(R=fixparms$m/2,C=fixparms$m/2), 
    parms=fixparms, 
    odes=dRCt_cont3, 
    ymin=0, ymax=NULL, log="", x=1, y=2, xlab="Time", ylab="Density", tmin=0, draw=lines, times=NULL, show=NULL, arrest=NULL, after=NULL, tweak=NULL, timeplot=TRUE, traject=FALSE, table=FALSE, add=FALSE, legend=TRUE, solution=FALSE, delay=FALSE, lwd=2, col="black", pch=20)


# Interference competition ------------------------------------------------

hi <- function(C,a,h,m){
  a*exp(-C) / (1 + a*h*exp(-C)) - m
}

curve(hi(x,1,0,1),xlim=c(0,10))

hi2 <- function(C,a,i,e,m){
  i*a / ((e+a)*C) - m
}

curve(hi2(x,1,1,1,1),xlim=c(0,10))

curve(x/(1+1),xlim=c(0,100))
curve(x/(1+2),add=T,lty=2)

# Beverton-Holt Consumer --------------------------------------------------

dC_BH <- function(C,i,e,alpha,a,mu){
  alpha * a * i / (e + a*C) - mu
}

curve(dC_BH(10^(x),i=0.5,e=1,alpha=1,a=2,mu=0.1),col="red",xlim=c(0,1))
curve(dC_BH(10^(x),i=0.5,e=1,alpha=1,a=1,mu=0.1),col="black",add=T)
curve(dC_BH(10^(x),i=0.5,e=1,alpha=1,a=0.5,mu=0.1),col="blue",add=T)
abline(h=0,lty=2,col="gray")

Cseq <- 10^(seq(0,1,length.out=100))
rmean <- apply(cbind(
  dC_BH(Cseq,i=0.5,e=1,alpha=1,a=2,mu=0.1),
  dC_BH(Cseq,i=0.5,e=1,alpha=1,a=1,mu=0.1),
  dC_BH(Cseq,i=0.5,e=1,alpha=1,a=0.5,mu=0.1)
  ),1,mean)
lines(log10(Cseq),rmean,col="gray")