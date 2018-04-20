model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    
    dR <- r*R*(1 - R/K) - a*R*N
    dN <- c*a*R*N - delta*N
    
    return(list(c(dR, dN)))  
  }) 
}  

p <- c(r=1,K=1,a=1,c=1,delta=0.5)
s <- c(R=1,N=0.01) 

run(after="parms[\"r\"]<-rnorm(1,mean=1,sd=0.1)")

run(after="if(t==20)state[\"N\"]<-0")
# Use arrest to handle events at time points within time steps:
run(50,arrest=33.14,after="if(t==33.14)state[\"N\"]<-0",table=T)

f <- newton(c(R=0.5,N=0.5))
run(state=f,after="state<-state+rnorm(2,mean=0,sd=0.01)",ymax=1)

# Here is an example of similar event handling in deSolve:
fun<-function(t, y, parms){y["N"]<-0;return(y)}
run(events=list(func=fun,time=20))

