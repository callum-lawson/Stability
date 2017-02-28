### PACKAGES

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

rplot <- function(zmu,pars,xmin=0,xmax=10,nx=10^4,...){
	nz <- length(zmu)
	np <- nrow(pars)
	xseq <- seq(xmin,xmax,length.out=nx)
	rarr <- array(NA,dim=c(nx,nz,np))
	karr <- array(NA,dim=c(1,nz,np))
	for(i in 1:np){
		for(j in 1:nz){
			rarr[,j,i] <- with(pars[i,], rcalc(z=zmu[j],xseq,b0,b1,b2,b3,b4))
			karr[1,j,i] <- with(pars[i,], -(b0+b1*zmu[j])/(b3+b4*zmu[j]) )
			}
		}

	rmat <- acast(melt(rarr,varnames=c("x","p","z")), x ~ z + p)
	kmat <- acast(melt(karr,varnames=c("x","p","z")), x ~ z + p)

	mymatplot(xseq,rmat,xlab="",ylab="",...)
	points(kmat,rep(0,length(kmat)),pch=rep(c(16,21),2),col=cseq,cex=1.3)
	abline(h=0,col="black",lty=3)
	}

rplot(zmu,pars2)

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

dplot <- function(xmat,xmin=0,xmax=10,ndens=512,bw,...){
	xplot <- seq(xmin,xmax,length.out=ndens)
	xdens <- apply(xmat,2,function(x) density(x,from=xmin,to=xmax,n=ndens,bw=bw)$y)
	mymatplot(xplot,xdens,xlab="",ylab="",...)
	}

lettlab <- function(i){
	mtext(bquote(bold((.(letters[i])))),side=3,line=0.5,xpd=T,adj=0.05)
	}

addxlab <- function(xname){
	mtext(xname,side=1,line=1,las=1,xpd=T)
	}

addtitlab <- function(titname){
	mtext(titname,side=3,line=3,las=1,xpd=T)
	}

addylab <- function(yname){
	mtext(yname,side=2,line=1.5,las=3,xpd=T)
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

# zmu <- c(1,2)
zmu <- c(1,2)
zsd <- c(1,1)
nz <- length(zmu)
nt <- 10^5
zsim <- matrix(NA,nr=nt,nc=nz)
zsim[] <- rnorm(nt*nz,rep(zmu,each=nt),sd=zsd)

cols <- c("red","blue")
ltys <- 1:2
cseq <- rep(cols,each=length(ltys))
lseq <- rep(ltys,times=length(cols))

pars1 <- data.frame(b0=c(4.5,2.25),b1=c(1,1),b2=c(0,0),b3=c(-1,-1),b4=c(0,0))
pars2 <- data.frame(b0=c(2.75,4),b1=c(1,1),b2=c(0,0),b3=c(-0.7,-1.7),b4=c(0,0))
pars3 <- data.frame(b0=c(4.5,0.5),b1=c(1,1),b2=c(0,0),b3=c(-0.5,-0.5),b4=c(-0.3,0))

### SIMS

xmat1 <- xsim(zsim,pars1,nt=nt)
xmat2 <- xsim(zsim,pars2,nt=nt)
xmat3 <- xsim(zsim,pars3,nt=nt)

### PLOTS

setwd("C:/Users/Callum/Dropbox/NIOO/Grants/Ambizione")
tiff(filename="ambizione_theory_illustration.tiff", width=7.5, height=4.5, units="in", res=200)

par(mfrow=c(2,3),mar=c(1,1.5,3,1.5),oma=c(3,3,3,3),
	bty="l",xaxt="n",yaxt="n"
	)
rplot(zmu,pars1)
lettlab(1)
addylab("Population growth rate")
addtitlab(expression(bold(Different~intercepts)))
rplot(zmu,pars2)
lettlab(2)
addtitlab(expression(bold(Different~slopes)))
rplot(zmu,pars3)
lettlab(3)
addtitlab(expression(bold(Different~interactions)))

bw <- 0.5
dplot(xmat1,bw=bw)
lettlab(4)
addylab("Probability density")
addxlab("Population density")
dplot(xmat2,bw=bw)
lettlab(5)
addxlab("Population density")
dplot(xmat3,bw=bw)
lettlab(6)
addxlab("Population density")

dev.off()

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

