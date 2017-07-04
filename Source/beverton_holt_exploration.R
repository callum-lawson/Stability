####################################################################
# Explore how random variation in Beverton-Holt mortality function #
# affects the distribution of population sizes                     #
####################################################################

# Functions ---------------------------------------------------------------

BH <- function(n_in,m0,m1){
  n_out <- n_in * exp(-m0) / ( 1 + (m1/m0)*(1-exp(-m0))*n_in )
  return(n_out)
}
# assuming T = 1

lnbh <- function(y,m0,m1,T3=0.3){
  c1 <- exp(-m0*T3)
  c2 <- (1-exp(-m0)*T3)*m1/m0
  log((y*c1-1)/(c2*y))
}

lnbhdor <- function(Y,m0,m1,G,T3=0.3){
  S0 <- exp(-m0)
  c1 <- exp(-m0*T3)
  c2 <- (1-exp(-m0)*T3)*m1/m0
  log((c1*G*Y-G*S0+S0-1) / (c2*G*Y*((G-1)*S0+1)))
}

Kcalc <- function(m0,m1,T=1){
  exp(-m0*T)/( (1-exp(-m0*T)) * (m1/m0) )
} 
# assuming T=1

sGBH <- function(n_in,Y,m0,m1,G,T3=0.3){
  So <- exp(-m0)
  Sn <- exp(-m0*T3)
  n_in*((1-G)*So + G*Y*Sn / (1 + G*Y*(1-Sn)*(m1/m0)*n_in))
}

rG <- function(Y,m0,G,T3=0.3){
  So <- exp(-m0)
  Sn <- exp(-m0*T3)
  log( (1-G)*So + G*Y*Sn )
}

# Exploration -------------------------------------------------------------

ln0_min <- -5
ln0_max <- 0
nseq <- 10^3
neps <- 10^3
alpha_m <- 0.1
beta_m <- 1

# n0_seq <- exp(seq(ln0_min,ln0_max,length.out=nseq))
# eps_m <- seq(-1,1,length.out=neps)
n0_seq <- exp(rnorm(nseq,0,1.5))
eps_m <- rnorm(neps,0,2)
  
m0 <- m1 <- n0 <- n1 <- matrix(NA,nr=nseq,nc=neps)

m0[] <- rep(exp(outer(alpha_m,eps_m,"+")),each=nseq)
m1[] <- rep(exp(outer(beta_m,eps_m,"+")),each=nseq)
n0[] <- rep(n0_seq,times=neps)
n1 <- BH(n0,m0,m1)

par(mfrow=c(2,1))
hist(n0,breaks=100,prob=T)
hist(n1,breaks=100,prob=T)

hist(n1[1,],breaks=100)
hist(n1[2,],breaks=100)
hist(n1[3,],breaks=100)

matplot(n0[,1],n1,type="l")
abline(0,1,col="purple",lty=3)

(K <- BH(10^10,m0[1,],m1[1,]))

# Plant example -----------------------------------------------------------
  # From Yodzis and work on desert annuals

### Without seed dormancy

curve(log(1-1/exp(x)),col="red",xlim=c(0,0.5))
curve(log(0.9-1/exp(x)),col="blue",add=T)
# drop-off occurs earlier at lower K

m0 <- 1.1
m1 <- 0.1
par(mfrow=c(1,1),mar=c(4,4,1,1))
curve(lnbh(x,m0,m1),col="blue",xlim=c(1,3))
m1 <- 0.01
curve(lnbh(x,m0,m1)-2.3,col="red",add=T)
# changing c2 only changes equilibrium N by constant factor (shifts intercept)

m0 <- 1
m1 <- 0.01
c1 <- exp(-m0)
c2 <- (1-exp(-m0))*m1/m0
curve(log((exp(x)*c1-1)/(c2*exp(x))),col="blue",xlim=c(0,3))
m0 <- 0.5
c1 <- exp(-m0)
c2 <- (1-exp(-m0))*m1/m0
curve(log((exp(x)*c1-1)/(c2*exp(x))),col="red",add=T)
# changing c1 -> drop-off occurs *earlier* at lower K

### With seed dormancy

G <- 0.4
m0 <- 1
m1 <- 0.01
curve(lnbhdor(x,m0,m1,G),col="blue",xlim=c(1,3))
m0 <- 0.5
curve(lnbhdor(x,m0,m1,G),col="red",add=T)
  # Same result with G - drop-off earlier for lower K

# How does Y affect BH curve? ---------------------------------------------

Y <- 3 # log scale
m0 <- 1
curve(log(BH(exp(x+Y),m0,m1)),col="blue",xlim=c(-1,5),
  xlab=expression(N[g]),ylab=expression(N[s])) # ignores dormancy
curve(log(BH(exp(x+Y-0.5),m0,m1)),col="red",add=T)
m0 <- 2
curve(log(BH(exp(x+Y),m0,m1)),col="blue",add=T,lty=2)
curve(log(BH(exp(x+Y-0.5),m0,m1)),col="red",add=T,lty=2)
abline(0,1,lty=3)
  # -0.5 -> reduce y by a freaction of exp(0.5)=1.6
  # decrease in y reduces equilibrium N by *more* when max popsize (K) is lower
  # Explantion:
  # - increasing Y makes a bigger difference to Ns when few seeds are being 
  # produced (bigger difference in [blue - red] at lower Ng (or Y)) 
  # - increasing m0 -> fewer seeds survive -> lower reproduction next year
  # -> bigger effect of change in Y
  # (but assumes that general shape of function stays the same? 
  # Increasing m0 makes curve steeper too, but this isn't important here,
  # as doesn't affect difference between blue and red much?)

# Y is just scaling up and down, changing c2 does the same thing?

# How does G affect BH curve? ---------------------------------------------

m0 <- 0.1
m1 <- 0.01
curve(lnbhdor(x,m0,m1,G=0.75),col="blue",xlim=c(0,5),ylim=c(-0.5,5))
curve(lnbhdor(x,m0,m1,G=0.5),col="red",add=T)
curve(lnbhdor(x,m0,m1,G=0.25),col="green",add=T)
  # above a certain Y, get more pop with high G, but as Y decreases, better to have low G

m0 <- 0.1
m1 <- 0.01
curve(lnbhdor(exp(x),m0,m1,G=0.75),col="blue",xlim=c(0,1.5),ylim=c(-0.5,5))
curve(lnbhdor(exp(x),m0,m1,G=0.25),col="red",add=T)
m0 <- 0.2
m1 <- 0.01
curve(lnbhdor(exp(x),m0,m1,G=0.75),col="blue",lty=2,add=T)
curve(lnbhdor(exp(x),m0,m1,G=0.25),col="red",lty=2,add=T)
  # When G is lower, reducing K has bigger effect on drop-off threshold for Y
  # (but will still always reduce it, so doesn't explain positive correlations)

m0 <- 0.1
m1 <- 0.01
curve(lnbhdor(x,m0,m1,G=0.75)-lnbhdor(x,m0,m1,G=0.25),col="blue",xlim=c(0,10))
m0 <- 0.2
m1 <- 0.01
curve(lnbhdor(x,m0,m1,G=0.75)-lnbhdor(x,m0,m1,G=0.25),col="red",add=T)
m0 <- 0.3
m1 <- 0.01
curve(lnbhdor(x,m0,m1,G=0.75)-lnbhdor(x,m0,m1,G=0.25),col="green",add=T)
abline(h=0,lty=3)
  # m1 -> lower germination favoured for given Y 

m0 <- 0.1
m1 <- 0.01
curve(lnbhdor(x,m0,m1,G=0.75)-lnbhdor(x,m0,m1,G=0.25),col="blue",xlim=c(0,10))
m0 <- 0.1
m1 <- 0.015
curve(lnbhdor(x,m0,m1,G=0.75)-lnbhdor(x,m0,m1,G=0.25),col="red",add=T)
abline(h=0,lty=3)
  # m2 -> proportional effect of G on dK is the same

# How does DI mortality relate to K? --------------------------------------

curve(Kcalc(exp(x),1),xlim=c(-5,5))
curve(Kcalc(exp(x),2),col="red",add=T)
  # higher m0 -> lower K

# Optimal G - constant environment ----------------------------------------

  # maximise equilibrium population size (ignoring extinction)

nY <- 5
nm0 <- 5
nm1 <- 5
nG <- 50

Yseq <- exp(seq(0,5,length.out=nY))
m0seq <- exp(seq(-1.5,0,length.out=nm0))
m1seq <- exp(seq(-4,4,length.out=nm1))
Gseq <- plogis(seq(-5,5,length.out=nG))

pd <- expand.grid(Y=Yseq,m0=m0seq,m1=m1seq,G=Gseq) # parameter data

lK_GY_plot <- function(m){
  require(reshape2)
  require(fields)
  pd$lK <- with(pd, lnbhdor(Y,m0,m1,G,T3=0.3))
  plotvars <- c("Y","G")
  pd1 <- subset(pd,m0==m & m1==m1seq[2],select=c(plotvars,"lK"))
  plotmat <- acast(melt(pd1,id=plotvars),G~Y)
  matplot(qlogis(Gseq),plotmat,type="l",lty=1,col=tim.colors(nY))
}

par(mfrow=c(2,2),mar=c(3,3,2,2))
lK_GY_plot(m0seq[1])  
lK_GY_plot(m0seq[2])  
lK_GY_plot(m0seq[3])  
lK_GY_plot(m0seq[4])  
  # m1 has no impact on optimal phenotype (only affects value of K)
  # G=0 is never optimal
  # Lower Y favours *higher* G -> should have negative G~z relationship?
  # (These plots show optimal G in years with different Y?)
  # Higher m0 favours higher G

lK_YG_plot <- function(m){
  require(reshape2)
  require(fields)
  pd$lK <- with(pd, lnbhdor(Y,m0,m1,G,T3=0.3))
  plotvars <- c("Y","G")
  pd1 <- subset(pd,m0==m & m1==m1seq[2],select=c(plotvars,"lK"))
  plotmat <- acast(melt(pd1,id=plotvars),Y~G)
  matplot(log(Yseq),plotmat,type="l",lty=1,col=tim.colors(nG))
}

par(mfrow=c(2,2),mar=c(3,3,2,2))
lK_YG_plot(m0seq[1])  
lK_YG_plot(m0seq[2])  
lK_YG_plot(m0seq[3])  
lK_YG_plot(m0seq[4])
  # low G -> 
  # - earlier drop-off in N with Y
  # - larger proportional drop-offs in N with Y
  # so *more* affected by CC if lower G (even without plasticity)
  # this tallies with "lower Y -> higher optimal G" results

# Optimal G - variable environment ----------------------------------------

np <- nrow(pd) # number of parameter combinations
nb <- 100
nt <- 1000 + nb
ns <- 10 # n y sig
Ysigseq <- exp(seq(-2,2,length.out=ns))
N0 <- 1000

Narr <- Yarr <- array(dim=c(nt,ns,np))
m0arr <- m1arr <- Garr <- array(dim=c(ns,np))
m0arr[] <- rep(pd$m0,each=ns)
m1arr[] <- rep(pd$m1,each=ns)
Garr[] <- rep(pd$G,each=ns)

Yarr[] <- exp(
  rep(log(pd$Y),each=nt*ns)
  + rep(outer(rnorm(n=nt,mean=0,sd=1),Ysigseq),times=np)
) # only works because Y varies the fastest of all params in pd

Narr[1,,] <- N0

for(t in 2:nt){
  Narr[t,,] <- sGBH(Narr[t-1,,],Yarr[t,,],m0arr,m1arr,Garr)
}

lNarr <- log(Narr)
lNmuarr <- apply(lNarr[-(1:nb),,],2:3,mean)
lNsdarr <- apply(lNarr[-(1:nb),,],2:3,sd)
lNmedarr <- apply(lNarr[-(1:nb),,],2:3,median)

Nd <- expand.grid(Ysig=Ysigseq,Y=Yseq,m0=m0seq,m1=m1seq,G=Gseq)
Nd$lNmu <- as.vector(lNmuarr)
Nd$lNsd <- as.vector(lNsdarr)
Nd$lNmed <- as.vector(lNmedarr)

lmu_Ysig_G_plot <- function(pY,pm0,pm1,xvar="G",...){
  require(reshape2)
  require(fields)
  plotvars <- c("Ysig","G")
  Nd1 <- subset(Nd,Y==pY & m0==pm0 & m1==pm1,select=c(plotvars,"lNmu"))
  if(xvar=="G"){
    plotmat <- acast(melt(Nd1,id=plotvars),G~Ysig)
    matplot(qlogis(Gseq),plotmat,type="l",lty=1,col=tim.colors(ns),
            ylim=c(-2,max(ceiling(plotmat))),
            ...)
    Gopt <- sapply(split(Nd1,Nd1$Ysig),function(d){
      qlogis(d$G[which(d$lNmu==max(d$lNmu))])
    })
    abline(v=Gopt,col=tim.colors(nY),lty=2)
  }
  if(xvar=="Ysig"){
    plotmat <- acast(melt(Nd1,id=plotvars),Ysig~G)
    matplot(Ysigseq,plotmat,type="l",lty=1,col=tim.colors(nG),
            ylim=c(-2,max(ceiling(plotmat))),
            ...)
  }
}

pdf(paste0("Plots/lmu_Ysig_G_",format(Sys.Date(),"%d%b%Y"),".pdf"),
    width=15,height=15)
iseq <- 1:5
ni <- length(iseq)
par(mfrow=c(ni,ni),mar=c(4,4,2,2),las=1,bty="l")
for(i in iseq){
  for(j in iseq){
    for(k in iseq){
      lmu_Ysig_G_plot(pY=Yseq[i],pm0=m0seq[j],pm1=m1seq[k],xvar="G",
                      xlab="G",ylab="ln(N)",
                      main=paste0("Y=",signif(Yseq[i],2),
                                  " m0=",signif(m0seq[j],2),
                                  " m1=",signif(m1seq[k],2)
                                  )
      )
    }
  }
}
dev.off()

par(mfrow=c(2,2),mar=c(4,4,2,2),las=1,bty="n")
lmu_Ysig_G_plot(pY=Yseq[4],pm0=m0seq[2],pm1=m1seq[2],xvar="G")
lmu_Ysig_G_plot(pY=Yseq[4],pm0=m0seq[4],pm1=m1seq[2],xvar="G")
lmu_Ysig_G_plot(pY=Yseq[2],pm0=m0seq[2],pm1=m1seq[4],xvar="G")
lmu_Ysig_G_plot(pY=Yseq[2],pm0=m0seq[2],pm1=m1seq[2],xvar="G")
  # - high m0 reduces effect of variance on optimal G, 
  # possibly even leads to env var favouring high G?
  # - m1 again just shifts population size without changing optimum
  # - low Y accentuates effects of env var on G?  
  # - low Y seems to favour higher G (again)

  # fluctuations favour G<1

par(mfrow=c(2,2))
lmu_Ysig_G_plot(pY=Yseq[5],pm0=m0seq[3],pm1=m1seq[3],xvar="Ysig")
lmu_Ysig_G_plot(pY=Yseq[5],pm0=m0seq[5],pm1=m1seq[3],xvar="Ysig")
lmu_Ysig_G_plot(pY=Yseq[5],pm0=m0seq[3],pm1=m1seq[5],xvar="Ysig")
lmu_Ysig_G_plot(pY=Yseq[2],pm0=m0seq[3],pm1=m1seq[3],xvar="Ysig")

lmu_Y_G_plot <- function(pYsig,pm0,pm1,xvar="G",...){
  require(reshape2)
  require(fields)
  plotvars <- c("Y","G")
  Nd1 <- subset(Nd,Ysig==pYsig & m0==pm0 & m1==pm1,select=c(plotvars,"lNmu"))
  if(xvar=="G"){
    plotmat <- acast(melt(Nd1,id=plotvars),G~Y)
    matplot(qlogis(Gseq),plotmat,type="l",lty=1,col=tim.colors(nY),
            ylim=c(-2,max(ceiling(plotmat))),
            ...
            )
    Gopt <- sapply(split(Nd1,Nd1$Y),function(d){
      qlogis(d$G[which(d$lNmu==max(d$lNmu))])
    })
    abline(v=Gopt,col=tim.colors(nY),lty=2)
  }
  if(xvar=="Y"){
    plotmat <- acast(melt(Nd1,id=plotvars),Y~G)
    matplot(Yseq,plotmat,type="l",lty=1,col=tim.colors(nG),
            ylim=c(-2,max(ceiling(plotmat))),
            ...
            )
  }
}

pdf(paste0("Plots/lmu_Y_G_",format(Sys.Date(),"%d%b%Y"),".pdf"),
    width=15,height=15)
iseq <- 1:5
ni <- length(iseq)
par(mfrow=c(ni,ni),mar=c(4,4,2,2),las=1,bty="l")
for(i in iseq*2){ # Ysig twice as long as others
  for(j in iseq){
    for(k in iseq){
      lmu_Y_G_plot(pYsig=Ysigseq[i],pm0=m0seq[j],pm1=m1seq[k],xvar="G",
                      xlab="G",ylab="ln(N)",
                      main=paste0("Ysig=",signif(Ysigseq[i],2),
                                  " m0=",signif(m0seq[j],2),
                                  " m1=",signif(m1seq[k],2)
                      )
      )
    }
  }
}
dev.off()

par(mfrow=c(2,2))
lmu_Y_G_plot(pYsig=Ysigseq[4],pm0=m0seq[2],pm1=m1seq[2],xvar="G")
lmu_Y_G_plot(pYsig=Ysigseq[4],pm0=m0seq[4],pm1=m1seq[2],xvar="G")
lmu_Y_G_plot(pYsig=Ysigseq[5],pm0=m0seq[2],pm1=m1seq[4],xvar="G")
lmu_Y_G_plot(pYsig=Ysigseq[2],pm0=m0seq[2],pm1=m1seq[2],xvar="G")
  # NB: outcome with median is almost exactly the same

# Analysis of high-variance scenario --------------------------------------

pd$K <- with(pd, Kcalc(m0,m1,T=0.3))

Darr <- Yearr <- array(dim=dim(Narr))
for(i in 2:nt){
  Darr[i,,] <- log(Narr[i,,]) - log(Narr[i-1,,])
  Yearr[i,,] <- log(Narr[i,,] - exp(-m0arr)*(1-Garr)*Narr[i-1,,]) - log(Garr*Narr[i-1,,])
}

pdsel <- with(pd, which(
  Y %in% Yseq[c(1,5)] & m0==m0seq[5] & m1==m1seq[2] & G %in% Gseq[c(10,20)]
  ))
Ysigsel <- 10

scol <- rep(c("blue","red"),each=2) # blue = low G, red = high G
slty <- rep(1:2,times=2) # - = low Y, --- = high Y

par(mfrow=c(1,1),mar=c(2,2,2,2))
matplot(log(Narr[900:1100,Ysigsel,pdsel]),type="l",col=scol,lty=slty)

plot(density(log(Yarr[,Ysigsel,pdsel[1]])))
lines(density(log(Yarr[,Ysigsel,pdsel[2]])),lty=2)

histmax <- 10^4
mybreaks <- c(seq(0,histmax,length.out=10^3),Inf)
hist(Yarr[,Ysigsel,pdsel[1]],xlim=c(0,histmax),breaks=mybreaks,col="blue",border=NA)
hist(Yarr[,Ysigsel,pdsel[2]],breaks=mybreaks,col="red",add=T,border=NA)

table(Yarr[,Ysigsel,pdsel[1]]>1)
table(Yarr[,Ysigsel,pdsel[2]]>1)
table(Yarr[,Ysigsel,pdsel[1]]>pd$K[pdsel[1]])
table(Yarr[,Ysigsel,pdsel[2]]>pd$K[pdsel[2]])

matplot(Darr[900:1100,Ysigsel,pdsel],type="l",col=scol,lty=slty)
matplot(Yearr[900:1100,Ysigsel,pdsel],type="l",col=scol,lty=slty)

plot(density(log(Narr[100:nt,Ysigsel,pdsel[4]])),col=scol[4],lty=slty[4])
lines(density(log(Narr[100:nt,Ysigsel,pdsel[2]])),col=scol[2],lty=slty[2])
plot(density(log(Narr[100:nt,Ysigsel,pdsel[3]])),col=scol[3],lty=slty[3])
lines(density(log(Narr[100:nt,Ysigsel,pdsel[1]])),col=scol[1],lty=slty[1])

plot(density(Darr[100:nt,Ysigsel,pdsel[4]]),col=scol[4],lty=slty[4])
lines(density(Darr[100:nt,Ysigsel,pdsel[2]]),col=scol[2],lty=slty[2])
abline(v=0,lty=3)
plot(density(Darr[100:nt,Ysigsel,pdsel[3]]),col=scol[3],lty=slty[3])
lines(density(Darr[100:nt,Ysigsel,pdsel[1]]),col=scol[1],lty=slty[1])
abline(v=0,lty=3)

plot(density(Yearr[100:nt,Ysigsel,pdsel[4]]),col=scol[4],lty=slty[4])
lines(density(Yearr[100:nt,Ysigsel,pdsel[2]]),col=scol[2],lty=slty[2])
abline(v=m0seq[5],lty=3)
plot(density(Yearr[100:nt,Ysigsel,pdsel[3]]),col=scol[3],lty=slty[3])
lines(density(Yearr[100:nt,Ysigsel,pdsel[1]]),col=scol[1],lty=slty[1])
abline(v=m0seq[5],lty=3)

# Density-independent -----------------------------------------------------

# > dim(lNmuarr)
# [1]   10 6250

pds <- expand.grid(Y=Yseq,m0=m0seq,G=Gseq) # parameter data
nps <- nrow(pds)

nlnY <- 10
lnYseq <- seq(-3,3,length.out=nlnY)
lnYarr <- rarr <- parr <- array(dim=c(nlnY,ns,nps))
lnYarr[] <- rep(outer(lnYseq,Ysigseq),times=nps) + rep(log(pds$Y),each=nlnY*ns)
rarr <- array(dim=c(nlnY,ns,nps))
rarr[] <- with(pds, rG(Y=exp(lnYarr),
                      m0=rep(m0,each=nlnY*ns),
                      G=rep(G,each=nlnY*ns)
                      ))
parr[] <- rep(dnorm(lnYseq,mean=0,sd=1),times=ns*nps) 
  # shape of p(norm) is same regardless of mean
sump <- sum(parr[,1,1]) 
  # again, distribution is the same for all, so the sum is too
rbararr <- apply(parr*rarr,c(2,3),sum)/sump

rd <- expand.grid(Ysig=Ysigseq,Y=Yseq,m0=m0seq,G=Gseq)
rd$rbar <- as.vector(rbararr)

Gopt_f <- function(d){
  d$G[which(d$rbar==max(d$rbar))]
}

rbar_Y_G_plot <- function(pYsig,pm0,xvar="G",...){
  require(reshape2)
  require(fields)
  plotvars <- c("Y","G")
  rd1 <- subset(rd,Ysig==pYsig & m0==pm0,select=c(plotvars,"rbar"))
  plotmat <- acast(melt(rd1,id=plotvars),G~Y)
  matplot(qlogis(Gseq),plotmat,type="l",lty=1,col=tim.colors(nY),
          # ylim=c(-2,max(ceiling(plotmat))),
          ...
          )
  Gopt <- sapply(split(rd1,rd1$Y),function(d){
    qlogis(d$G[which(d$rbar==max(d$rbar))])
  })
  abline(v=Gopt,col=tim.colors(nY),lty=2)
}
  # m1 actually unecessary as no DD

pdf(paste0("Plots/rbar_Y_G_",format(Sys.Date(),"%d%b%Y"),".pdf"),
    width=15,height=15)
iseq <- 1:5
ni <- length(iseq)
par(mfrow=c(ni,ni),mar=c(4,4,2,2),las=1,bty="l")
for(i in iseq*2){ # Ysig twice as long as others
  for(j in iseq){
    rbar_Y_G_plot(pYsig=Ysigseq[i],pm0=m0seq[j],xvar="G",
                  xlab="G",ylab="r",
                  main=paste0("Ysig=",signif(Ysigseq[i],2)," m0=",signif(m0seq[j],2))
                  )
    
  }
}
dev.off()


# Simple Royama plots -----------------------------------------------------

RoyBH <- function(G,lNmin,lNmax){
  Y1 <- 5 # 4 # log scale
  Y2 <- 1 # 3.5 # log scale
  dY <- 3.5
  m0 <- 1
  nlnN <- 100
  lnNseq <- seq(lNmin,lNmax,length.out=nlnN)
  Nseq <- exp(lnNseq)
  
  rd <- data.frame(
    Y = rep(c(Y1,Y2),each=3),
    dY = rep(c(-dY,0,dY),times=2)
  )
  rm <- matrix(nr=nlnN,nc=6)
  for(i in 1:nrow(rd)){
    rm[,i] <- with(rd[i,], log(sGBH(Nseq,exp(Y+dY),m0,m1,G)/Nseq))
  }
  
  par(mfrow=c(1,1))
  matplot(lnNseq,rm,col=rep(c("red","blue"),each=3),lty=rep(c(2,1,2),times=2),type="l")
  abline(h=0,lty=3)
}

RoyBH(G=0.15,lNmin=4,lNmax=8)
RoyBH(G=0.2,lNmin=4,lNmax=8)
RoyBH(G=0.5,lNmin=4,lNmax=8)

# Royama plots ------------------------------------------------------------

par(mfrow=c(2,2))
nN <- 100
Nseq <- exp(seq(0,5,length.out=nN))

r_sGBH_seq <- function(Y,m0,m1,G){
  log(sGBH(Nseq,Y,m0,m1,G)/Nseq)
}

Royama_matplot <- function(ppd,...){
  require(fields)
  require(plyr)
  colourvar <- ppd[,which(sapply(ppd,function(x) length(unique(x)))>1)]
  ncolours <- length(colourvar)
  rd <- mdply(ppd, r_sGBH_seq)
  rarr <- t(as.matrix(subset(rd,select=names(rd)[!names(rd) %in% names(pd)])))
  matplot(log(Nseq),rarr,type="l",lty=1,col=tim.colors(ncolours),...)
  abline(h=0,col="black",lty=3)
}

pdr <- pd[,names(pd)!="K"]

par(mfrow=c(2,2))
Royama_matplot(subset(pdr[,names(pd)!="K"], Y==Yseq[3] & m0==m0seq[3] & m1==m1seq[3]))
Royama_matplot(subset(pdr[,names(pd)!="K"], Y==Yseq[5] & m0==m0seq[3] & m1==m1seq[3]))
Royama_matplot(subset(pdr[,names(pd)!="K"], Y==Yseq[3] & m0==m0seq[5] & m1==m1seq[3]))
Royama_matplot(subset(pdr[,names(pd)!="K"], Y==Yseq[3] & m0==m0seq[3] & m1==m1seq[5]))
  # - increasing G increases stability but generally reduces K (unless at low G)
  # - when seed mortatlity is high, G doesn't affect K much

Royama_matplot(subset(pdr, m0==m0seq[3] & m1==m1seq[3] & G==Gseq[5]),ylim=c(-1,2))
Royama_matplot(subset(pdr, m0==m0seq[3] & m1==m1seq[3] & G==Gseq[20]),ylim=c(-1,2))
Royama_matplot(subset(pdr, m0==m0seq[3] & m1==m1seq[3] & G==Gseq[30]),ylim=c(-1,2))
Royama_matplot(subset(pdr, m0==m0seq[3] & m1==m1seq[3] & G==Gseq[50]),ylim=c(-1,2))
  # - env var has more positive effects when G is low and when population
  # size is low (DD-clim interaction)
  # - declines in reproduction can have *less* severe effects when G is high?

  # (for squared rainfall, can just compare what happens if return to low Y again)

Royama_matplot(subset(pdr, m0==m0seq[2] & m1==m1seq[3] & G==Gseq[25]))
Royama_matplot(subset(pdr, m0==m0seq[3] & m1==m1seq[3] & G==Gseq[25]))
Royama_matplot(subset(pdr, m0==m0seq[4] & m1==m1seq[3] & G==Gseq[25]))
Royama_matplot(subset(pdr, m0==m0seq[5] & m1==m1seq[3] & G==Gseq[25]))
  # higher mortality -> lower K but overall patterns pretty similar?

Royama_matplot(subset(pdr, m0==m0seq[3] & m1==m1seq[2] & G==Gseq[25]))
Royama_matplot(subset(pdr, m0==m0seq[3] & m1==m1seq[3] & G==Gseq[25]))
Royama_matplot(subset(pdr, m0==m0seq[3] & m1==m1seq[4] & G==Gseq[25]))
Royama_matplot(subset(pdr, m0==m0seq[3] & m1==m1seq[5] & G==Gseq[25]))
  # higher m2 -> higher stability, lower K 
  # (just shifts intercept, effect happens because of non-linear DD?)

# Optimal plastic G - variable environment --------------------------------

  # (no need for plasticity in constant environment)

alpha_G <- beta_Gz <- matrix(NA,nr=np,nc=np) 
alpha_G[] <- -tau_mu*sqrt(pi^2/(3*tau_sd^2))
beta_Gz[] <- sqrt(pi^2/(3*tau_sd^2))
# based on Godfray & Rees 2002

  # - plot fitness distribution over different environments, calculate width to
  # quantify degree of bet-hedging
  # - can also infer effects of plasticity by looking at plots showing response 
  # to environment for different germination fractions (effectively 3d plots for
  # cue / environment)
