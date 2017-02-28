##############################################################################
### Examine ideas behind "nonlinear perturbations" in Royama (1992) p38-40 ###
##############################################################################

psurv <- function(fpi) 1-exp(-fpi) # food per indiv; taken from hole example
Tsize <- 100

fpop <- function(nstart,fstart,T=Tsize,vecout=F){

	nvec <- fvec <- vector()
	nvec[1] <- nstart
	fvec[1] <- fstart

	for(t in 2:T){
		nvec[t] <- nvec[t-1]*psurv(fvec[t-1]/nvec[t-1])
			# n indivs to get food
		fvec[t] <- fvec[t-1] - nvec[t] 
			# n food items eaten = nindivs to get food
		}
	
	if(vecout==T) list(nvec=nvec[T],fvec=fvec[T])
	if(vecout==F) nvec[T]

	}
	
nn <- 1000
nf <- 5
nsvec <- exp(seq(-1,8,length.out=nn))
fsvec <- exp(seq(5,10,length.out=nf))

nfinmat <- matrix(nr=nn,nc=nf)

for(i in 1:nn){
	for(j in 1:nf){
		nfinmat[i,j] <- fpop(nsvec[i],fsvec[j])
		}
	}

pdf(paste0("Plots/royama_nonlin_peturb_",format(Sys.Date(),"%d%b%Y"),".pdf"),width=5.25,height=5.5)
mycols <- 1:nf
matplot(log(nsvec),nfinmat/nsvec,type="l",col=mycols,lty=1,
	xlab="start N (log scale)",ylab="prop surv",bty="l")
legend("topright",legend=signif(fsvec,2),col=mycols,title="start food",lty=rep(1,nf),bty="n")
dev.off()




