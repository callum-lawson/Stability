model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    
    dR <- r*R*(1 - R/K) - a*R*N/(h+R)
    dN <- c*a*R*N/(h+R) - delta*N
    
    return(list(c(dR, dN)))  
  }) 
}  

p <- c(r=1,K=1,h=0.1,a=0.5,c=1,delta=0.4)
s <- c(R=1,N=0.01)

plane(eps=-0.001)
run(200,traject=TRUE)
