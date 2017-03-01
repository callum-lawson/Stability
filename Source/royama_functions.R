######################################################################################## Functions to simulate and plot population dynamics with Gompertz density-dependence #
#######################################################################################

# Graphical parameters ----------------------------------------------------

lineylim <- c(-1,1)
densylim <- c(0,0.75)

yval <- 0.325
sarr <- c(0.065,0.065)
mycex <- 0.8
mymult <- 0.8

# Simulation functions ----------------------------------------------------

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
