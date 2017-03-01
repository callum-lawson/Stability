#################################################################################
# Population simulations under different climate DD effects, for grant proposal #
#################################################################################

source("Source/royama_functions.R")

###################
### PLOT GRAPHS ###
###################

# b0 = intercept
# b1 = linear climate 	
# b2 = squared climate 	
# b3 = density 
# b4 = climate*density

### PARAMS 

set.seed(5)

zmu <- c(1,2)
zsd <- c(0.25,0.25)
zdem <- zmu # demonstrated relationships
nz <- length(zmu)
nt <- 10^4
zsim <- matrix(NA,nr=nt,nc=nz)
zsim[] <- rnorm(nt*nz,rep(zmu,each=nt),sd=zsd)
# zsim[] <- filter(rnorm(nt*nz,rep(zmu,each=nt),sd=zsd), filter=rep(1,3), circular=TRUE)

myxmin <- 0.5 # 0
myxmax <- 5.5 # 7

cseq <- c("red","blue")
lseq <- 2:1

pars1 <- data.frame(b0=0,b1=1,b2=0,b3=-1/2,b4=0)
pars2 <- data.frame(b0=2,b1=-1/2,b2=0,b3=-5/4,b4=1/2)
pars3 <- data.frame(b0=-2,b1=2.5,b2=0,b3=1/4,b4=-1/2)
  # DD WAS -1; interaction WAS -0.5

K1 <- Kcalc(z=zmu,pars=pars1)
K2 <- Kcalc(z=zmu,pars=pars2)
K3 <- Kcalc(z=zmu,pars=pars3)

### SIMS

xmat1 <- xsim(zsim,pars1,nt=nt)
xmat2 <- xsim(zsim,pars2,nt=nt)
xmat3 <- xsim(zsim,pars3,nt=nt)

### PLOTS

# Run butterfly density code first!

# setwd("C:/Users/Callum/Dropbox/NIOO/Grants/Marie Curie")
# setwd("D:/Users/calluml/Dropbox/NIOO/Grants/Humboldt")

lineylim <- c(-1,1)
densylim <- c(0,0.75)

tiff(filename=paste0("mariecurie_theory_illustration_nodensity_",format(Sys.Date(),"%d%b%Y"),".tiff"),
  width=9.5, height=5.25, units="in", res=600)

yval <- 0.325
sarr <- c(0.065,0.065)
mycex <- 0.8
mymult <- 0.8

par(mfrow=c(2,4),mar=c(2,2,1,2)+0.1,oma=c(3,3,4,3),
  bty="l",xaxt="n",yaxt="n"
	)
rplot(zdem,pars1,ylim=lineylim,sarr=sarr,yval=yval)
lettlab(1)
addylab(expression(log~population~growth~rate~(yr^-1)))
addtitlab(expression(bold(independent)))
rplot(zdem,pars2,ylim=lineylim,sarr=sarr,yval=yval)
lettlab(2)
addtitlab(expression(bold(counteracting)))
rplot(zdem,pars3,ylim=lineylim,sarr=sarr,yval=yval)
lettlab(3)
addtitlab(expression(bold(compounding)))

par(mfg=c(2,1),xaxt="s",yaxt="s")

gatdat <- splotbyclim(spl[[13]],names(spl)[13],climname="temp4",incldat=T,thomas=F,ndens=ndens,
  axlabs=F)
lettlab(4)

addylab(expression(log~population~growth~rate~(yr^-1)))
addxlab(expression(log~population~density~(m^-2)))

tortdat <- splotbyclim(spl[[45]],names(spl)[45],climname="temp1",incldat=T,thomas=F,ndens=ndens,
  axlabs=F)
lettlab(5)

addxlab(expression(log~population~density~(m^-2)))

paintdat <- splotbyclim(spl[[29]],names(spl)[29],climname="temp2",incldat=T,thomas=F,ndens=ndens,
  axlabs=F)
lettlab(6)

addxlab(expression(log~population~density~(m^-2)))

par(mar=c(0,0,0,0),mfg=c(1,4))
plot(1,1,type="n",xlab="",ylab="n",xaxt="n",yaxt="n",bty="n",ann=F)
legend("left",pch=c(16,21),col=c("blue","red"),
  legend=c("cold","warm"),
  title="Temperature",bty="n",cex=1.25)

dev.off()

##################################################################

zfake <- seq(1,2,length.out=100)

Vcalc <- function(z,pars){
  with(pars, b3+b4*z)
  }

Ks1 <- Kcalc(z=zfake,pars=pars1)
Ks2 <- Kcalc(z=zfake,pars=pars2)
Ks3 <- Kcalc(z=zfake,pars=pars3)

Vs1 <- Vcalc(z=zfake,pars=pars1)
Vs2 <- Vcalc(z=zfake,pars=pars2)
Vs3 <- Vcalc(z=zfake,pars=pars3)

plot(abs(Vs1)~Ks1,type="b")
plot(abs(Vs2)~Ks2,type="b")
plot(abs(Vs3)~Ks3,type="b")

#####################
### SUMMARY STATS ###
#####################

apply(xmat,2,mean)
apply(xmat,2,sd)
	# order: clim, pars

apply(rmat,2,function(x) xseq[which(abs(x)==min(abs(x)))])

####################
### SPATIAL SIMS ###
####################

# simulate sites with different mean climates
# simulate different responses among those same sites
# increase the climates in all sites by a degree
# see if spatial patterns predicted temporal ones

nt_k <- 1100
nk <- 100

eps_k <- rnorm(nk,0,1)
zmain_k_1 <- rnorm(nt_k,mean=1,sd=1)
zmain_k_2 <- rnorm(nt_k,mean=2,sd=1)
zsim_k_1 <- outer(zmain_k_1,eps_k,"+")
zsim_k_2 <- outer(zmain_k_2,eps_k,"+")

# pars_k <- data.frame(b0=1,b1=1,b2=0,b3=-1,b4=rnorm(nk,-0.05,0.05))
pars_k <- data.frame(b0=0,b1=1+rnorm(nk,0,1),b2=0,b3=-1,b4=0)
pars_k <- data.frame(b0=0,b1=3+eps_k,b2=0,b3=-1,b4=0)
pars_k <- data.frame(b0=0,b1=1,b2=0,b3=-1,b4=eps_k/50-0.1)

xmat_k_1 <- xsim(zsim_k_1,pars=pars_k,nt=nt_k,outmat=F,matchedpars=T)
xmat_k_bar_1 <- apply(xmat_k_1,2,mean)
xmat_k_2 <- xsim(zsim_k_2,pars=pars_k,nt=nt_k,outmat=F,matchedpars=T)
xmat_k_bar_2 <- apply(xmat_k_2,2,mean)

plot(c(eps_k,eps_k+1),c(xmat_k_bar_1,xmat_k_bar_2),pch="+",col=rep(c("blue","red"),each=nk))
plot(eps_k,xmat_k_bar_2-xmat_k_bar_1,pch="+",col=rep(c("blue","red"),each=nk))

