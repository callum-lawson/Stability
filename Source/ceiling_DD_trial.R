###############################################################
### Exploration of DD relationship based on "ceiling" model ###
###############################################################

# Life-cycle:
# In summer, have births (b) according to climate
# In autumn, winter, spring, births are reduced to some max available spaces K
# c.f. lottery model?

N0 <- 1
K <- 2000
b <- 3
d <- 1
zmu <- 1
zsd <- 1

nt <- 10^5
t <- 1:nt
N <- vector(mode="numeric",length=nt)
N[1] <- N0
z <- rnorm(nt,zmu,zsd)
b <- exp(b*z)

for(i in 2:nt){
	Ndash <-  N[i-1]*b[i-1]
	N[i] <- ifelse(Ndash<K, Ndash, K)
	}

plot(tail(log(N),500)~tail(t,500),type="l")
plot(log(N[-1]/N[-nt])~log(N[-nt]))
abline(h=0,col="red",lty=2)
lines(supsmu(log(N[-nt]),log(N[-1]/N[-nt]),bass=10),col="blue")

### ANALYTICAL B-H

bh <- function(N,rmax=5,K=10^3){
	rmax*(1-N/K) # from wikipedia (/N)
	}
Nseq <- exp(seq(log(10),8,length.out=10^4))
rseq <- bh(Nseq)
plot(rseq~log(Nseq),type="l",xlab="ln(N)",ylab="r",las=1,bty="l")
abline(h=0,col="red",lty=2)




