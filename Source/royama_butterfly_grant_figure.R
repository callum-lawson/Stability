#################################################################################
# Population simulations under different climate DD effects, for grant proposal #
#################################################################################

library(reshape2)

#################
### FUNCTIONS ###
#################

mymatplot <- function(m,...) matplot(m,type="l",col=cseq,lty=lseq,...)

rcalc <- function(z,lN,b0,b1,b2,b3,b4){
	b0 + b1*z + b2*z^2 + b3*lN + b4*z*lN
	}

popsim <- function(zmat,lN0,nt,b0,b1,b2,b3,b4,warmup){
	nz <- ncol(zmat)
	lNt <- matrix(NA,nr=nt,nc=nz)
	lNt[1,] <- lN0 
	for(t in 2:nt){
		lNt[t,] <- lNt[t-1,] + rcalc(z=zmat[t,],lN=lNt[t-1,],b0=b0,b1=b1,b2=b2,b3=b3,b4=b4)
		}
	keep <- (warmup+1):nt
	return(lNt[keep,])
	}
	# Gompertz
	# zmat: rows=t, cols=zmu
	# no noise (eps term)
	# assumes zsd is fixed

Kcalc <- function(z,pars){
  with(pars, - (b0+b1*z+b2*z^2) / (b3+b4*z) )
  }
  # slope scales rate of change in K:  (b3+b4*z)

rplot <- function(zmu,pars,xmin=myxmin,xmax=myxmax,sarr,yval,nx=10^4,...){
  nz <- length(zmu)
  np <- nrow(pars)
  xseq <- seq(xmin,xmax,length.out=nx)
  rarr <- array(NA,dim=c(nx,nz,np))
  karr <- array(NA,dim=c(1,nz,np))
  for(i in 1:np){
    for(j in 1:nz){
      rarr[,j,i] <- with(pars[i,], rcalc(z=zmu[j],xseq,b0,b1,b2,b3,b4))
      karr[1,j,i] <- Kcalc(zmu[j],pars[i,])
      }
    }
  
  rmat <- acast(melt(rarr,varnames=c("x","p","z")), x ~ z + p)
  kmat <- acast(melt(karr,varnames=c("x","p","z")), x ~ z + p)
  
  mymatplot(xseq,rmat,xlab="",ylab="",...)
  
  points(kmat,rep(0,length(kmat)),pch=c(16,16),col=cseq,cex=1.3)
  text(x=kmat,
    y=rep(0,length(kmat)),
    labels=c(expression(K[2]),expression(K[1])),
    col=cseq,
    cex=1.3,
    pos=3
    )
  
  xbot <- xseq[apply(rmat,2,function(x) which(abs(x-yval)==min(abs(x-yval))))]
  xtop <- xseq[apply(rmat,2,function(x) which(abs(x+yval)==min(abs(x+yval))))]
  
  arrows(x0=xbot,
    x1=xtop,
    y0=rep(0,4),
    code=3,
    col=cseq,
    length=sarr
    )
  arrows(x0=c(xbot,xtop),
    x1=c(xbot,xtop),
    y0=rep(0,4),
    y1=rep(c(yval,-yval),each=2),
    length=0,
    col=rep(cseq,2),
    lty=3
    )
  abline(h=0,col="black",lty=3)
  
  }

### SIMULATE DYNAMICS

xsim <- function(zmat,pars,nt=nt,lN0=7,warmup=100,outmat=T,matchedpars=F){
	
	nz <- ncol(zmat)
	np <- nrow(pars)
	
	if(matchedpars==F){
		xarr <- array(NA,dim=c(nt-warmup,nz,np))
		for(i in 1:np){
			xarr[,,i] <- with(pars[i,], 
				popsim(zmat=zmat,lN0=lN0,nt=nt,b0,b1,b2,b3,b4,warmup=warmup)
				)
			}
			# zmu (but not pars) handled simultaneously
		if(outmat==T){
			xmat <- acast(melt(xarr,varnames=c("x","p","z")), x ~ z + p)
			return(xmat)
			}
		else{
			return(xarr)
			}
		}

	if(matchedpars==T){
		if(nz!=np) warning("clim and par dims differ")
		xarr <- array(NA,dim=c(nt-warmup,nz)) # nz=np
		for(i in 1:np){
			xarr[,i] <- with(pars[i,], 
				popsim(zmat=cbind(zmat[,i]),lN0=lN0,nt=nt,b0,b1,b2,b3,b4,warmup=warmup)
				)
			}
		return(xarr)
		}
	}

# TIME SERIES
# matplot(xmat,type="l")

# LN(N) DENSITY DISTRIBUTIONS

# dplot(xmat1,K1,harr=harr,sarr=sarr,bw=bw,ylim=densylim)

dplot <- function(xmat,K,xmin=myxmin,xmax=myxmax,harr,sarr,ndens=2^9,bw,...){
  xplot <- seq(xmin,xmax,length.out=ndens)
  xdens <- apply(xmat,2,function(x){
    density(x,from=xmin,to=xmax,n=ndens,bw=bw)$y
    })
  y1 <- xdens[xplot<=K,]
  y2 <- xdens[xplot>K,]
  yposf <- function(y) which(abs(y-harr)==min(abs(y-harr))[1])
  y1pos <- apply(y1,2,yposf)
  y2pos <- apply(y2,2,yposf)
  x1 <- c(xplot[xplot<K][y1pos[1]]+parr,xplot[xplot>K][y2pos[1]]-parr)
  x2 <- c(xplot[xplot<K][y1pos[2]]+parr,xplot[xplot>K][y2pos[2]]-parr)
  mymatplot(xplot,xdens,xlab="",ylab="",...)
  arrows(x0=x1[1],x1=x1[2],y0=harr,length=sarr[1],code=3,col=cseq[1])
  arrows(x0=x2[1],x1=x2[2],y0=harr,length=sarr[2],code=3,col=cseq[2])
  text(x=mean(x1),y=harr,labels=expression(sigma[1]),pos=3,col=cseq[1])
  text(x=mean(x2),y=harr,labels=expression(sigma[2]),pos=3,col=cseq[2])
  }

mycex <- 0.85

lettlab <- function(i){
	mtext(bquote(bold((.(letters[i])))),side=3,line=0.5,xpd=T,adj=0.05,cex=1*mymult)
	}

addxlab <- function(xname){
	mtext(xname,side=1,line=3,las=1,xpd=T,cex=mycex)
	}

addtitlab <- function(titname){
	mtext(titname,side=3,line=3,las=1,xpd=T,cex=0.9)
	}

addylab <- function(yname){
	mtext(yname,side=2,line=3,las=3,xpd=T,cex=mycex)
	}

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

