model <- function(t, state, parms) {
  with(as.list(c(state,parms)), { 
    dN <- r*N*(1 - N) - N  
    return(list(c(dN)))  
  }) 
}  

p <- c(r=3.75)
s <- c(N=0.01)
run(100,method="euler")

# Make the famous bifurcation diagram
plot(range(2,4),range(0,1),type='n',xlab="r",ylab="N")
npoints <- 100; burnin <- 500
for (r in seq(2,4,0.01)) {
  p["r"] <- r
  data <- run(burnin+npoints,method="euler",table=TRUE,timeplot=FALSE)
  points(rep(r,npoints),data$N[(burnin+1):(burnin+npoints)],pch=".")
}

# Make a Takens reconstruction
p <- c(r=3.75)
s <- c(N=0.01)
data <- run(1000,method="euler",table=TRUE)
plot(data$N[1:999],data$N[2:1000],pch=".")


# PS here is a slightly faster version not using a for-loop:
plot(range(2,4),range(0,1),type='n',xlab="r",ylab="N")
fun <- function(r,npoints=50,burnin=500) return(run(burnin+npoints,parms=c(r=r),method="euler",table=TRUE,timeplot=FALSE)$N[(burnin+1):(burnin+npoints)])
r <- seq(2,4,0.01)
npoints <- 100; burnin <- 500
points(rep(r,each=npoints),sapply(r,fun,npoints,burnin),pch=".")
