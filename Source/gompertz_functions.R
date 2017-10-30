# Functions to simulate and plot population dynamics with #
# Gompertz density-dependence                             #

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

rcalc <- function(z,lN,pars){
  with(pars,{
    if(b4==0){
      return(b0 + b1*z + b2*z^2 + b3*lN + b4*lN^2 + b5*z*lN)
    }
    if(b4!=0){
      lNmax <- -b3/(2*b4) # value of lN which maximises r
      return(ifelse(lN < lNmax,
                  b0 + b1*z + b2*z^2 + b3*lNmax + b4*lNmax^2 + b5*z*lNmax,
                  b0 + b1*z + b2*z^2 + b3*lN + b4*lN^2 + b5*z*lN
                  ))
    }
  })
}

Kcalc <- function(z,pars){
  if(pars$b4==0){
    with(pars, - (b0+b1*z+b2*z^2) / (b3+b5*z) )
  }
  if(pars$b4!=0){
    with(pars, -( sqrt( (b3+b5*z)^2 -4*b4*(b0+z*(b1+b2*z)) ) + b3 + b5*z ) / (2*b4) )
  }
}
# slope scales rate of change in K:  (b3+b5*z)

prcalc <- function(z,zmu,zsd,lN,pars){
  dnorm(z,zmu,zsd) * rcalc(z,lN,pars)
}

pKcalc <- function(z,zmu,zsd,pars){
  dnorm(z,zmu,zsd) * Kcalc(z,pars)
}

intrcalc <- function(lN,zmu,zsd,pars){
  integrate(prcalc,
            lower=-Inf,upper=Inf,
            zmu=zmu,zsd=zsd,
            lN=lN,
            pars=pars
  )$value
}

intKcalc <- function(zmu,zsd,pars){
  integrate(pKcalc,
            lower=-Inf,upper=Inf,
            zmu=zmu,zsd=zsd,
            pars=pars
  )$value
}

# Simulation functions ----------------------------------------------------

popsim <- function(zmat,lN0,nt,pars,warmup){
  nz <- ncol(zmat)
  lNt <- matrix(NA,nr=nt,nc=nz)
  lNt[1,] <- lN0 
  for(t in 2:nt){
    lNt[t,] <- lNt[t-1,] + rcalc(z=zmat[t,],lN=lNt[t-1,],pars=pars)
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
                        popsim(zmat=zmat,lN0=lN0,nt=nt,pars[i,],warmup=warmup)
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

rplot_3eg <- function(zmu,zsd,pars,xmin,xmax,averages=FALSE,nx=10^2,...){
  require(reshape2)
  
  zmu <- zmu[1]
  zsd <- zsd[1]
  # using first climate mean and sd as example
  zseq <- c(zmu-zsd,zmu,zmu+zsd)
  nz <- length(zseq)
  np <- nrow(pars)
  cvals <- c("blue","red","black")[1:np]
  cseq <- rep(cvals,each=3) # up to three line types
  lseq <- rep(c(2,1,2),times=np)
  xseq <- seq(xmin,xmax,length.out=nx)
  rarr <- array(NA,dim=c(nx,nz,np))
  kmat <- array(NA,dim=c(nz,np))
  armat <- array(NA,dim=c(nx,np))
  akseq <- rep(NA,np)
  for(i in 1:np){
    for(j in 1:nz){
      rarr[,j,i] <- rcalc(z=zseq[j],xseq,pars[i,])
      kmat[j,i] <- Kcalc(zmu[j],pars[i,])
    }
    if(averages==TRUE){
      for(j in 1:nx){
        armat[j,i] <- intrcalc(xseq[j],zmu,zsd,pars[i,])
      }
      akseq[i] <- intKcalc(zmu,zsd,pars[i,])
    }
  }
  
  rmat <- acast(melt(rarr,varnames=c("t","p","z")), t ~ z + p)

  mymatplot(xseq,rmat,xlab="",ylab="",col=cseq,lty=lseq,...)
  if(averages==TRUE){
    matplot(xseq,armat,type="l",lty=3,add=T,col=cvals)
    points(akseq,rep(0,np),col=cvals,pch=16)
  } 
  abline(h=0,lty=3)
  
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

