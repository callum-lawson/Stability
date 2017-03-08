######################################################################################## Functions to simulate and plot population dynamics with Gompertz density-dependence #
#######################################################################################

# General graphical functions ---------------------------------------------

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

# Analytical calculations -------------------------------------------------

mymatplot <- function(m,...) matplot(m,type="l",...)

rcalc <- function(z,lN,b0,b1,b2,b3,b4){
  b0 + b1*z + b2*z^2 + b3*lN + b4*z*lN
}

Kcalc <- function(z,pars){
  with(pars, - (b0+b1*z+b2*z^2) / (b3+b4*z) )
}
# slope scales rate of change in K:  (b3+b4*z)


# Simulation functions ----------------------------------------------------

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

xsim <- function(zmat,pars,nt=nt,lN0=7,warmup=100,outmat=F,matchedpars=F){
  
  require(reshape2)
  
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
      xmat <- acast(melt(xarr,varnames=c("t","p","z")), t ~ z + p)
      return(xmat)
    }
    else{
      return(xarr)
    }
  }
  
}

# Results-plotting functions ----------------------------------------------

rplot_3eg <- function(zmu,zsd,pars,xmin=myxmin,xmax=myxmax,nx=10^4,...){
  zmu <- zmu[1]
  zsd <- zsd[1]
  # using first climate mean and sd as example
  zseq <- c(zmu-zsd,zmu,zmu+zsd)
  nz <- length(zseq)
  np <- nrow(pars)
  cseq <- rep(c("blue","red","black")[1:np],each=3) # up to three line types
  lseq <- rep(c(2,1,2),times=np)
  xseq <- seq(xmin,xmax,length.out=nx)
  rarr <- array(NA,dim=c(nx,nz,np))
  karr <- array(NA,dim=c(1,nz,np))
  for(i in 1:np){
    for(j in 1:nz){
      rarr[,j,i] <- with(pars[i,], rcalc(z=zseq[j],xseq,b0,b1,b2,b3,b4))
      karr[1,j,i] <- Kcalc(zmu[j],pars[i,])
    }
  }
  
  rmat <- acast(melt(rarr,varnames=c("t","p","z")), t ~ z + p)
  kmat <- acast(melt(karr,varnames=c("t","p","z")), t ~ z + p)
  
  mymatplot(xseq,rmat,xlab="",ylab="",col=cseq,lty=lseq,...)
  
}

dplot <- function(xmat,nz,np,xmin=myxmin,xmax=myxmax,harr,sarr,bw=NULL,ndens=2^9,...){
  # nz <- length(grep("1_",colnames(xmat)))
  # np <- ncol(xmat)/nz
  xplot <- seq(xmin,xmax,length.out=ndens)
  xdens <- apply(xmat,2,function(x){
    if(is.null(bw)) density(x,from=xmin,to=xmax,n=ndens,na.rm=T)$y
    if(!is.null(bw)) density(x,from=xmin,to=xmax,n=ndens,bw=bw,na.rm=T)$y
  })
  matplot(xplot,xdens,type="l",
          ...
          )
}

