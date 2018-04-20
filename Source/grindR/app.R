source("grind.R")

model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    dR <- r*R*(1 - R/K) - a*R*N
    dN <- c*a*R*N - delta*N
    return(list(c(dR, dN)))  
  }) 
}  
p <- c(r=1,K=1,a=1,c=1,delta=0.5) 
s <- c(R=1,N=0.1)

library(shiny)

ui <- fluidPage(
  titlePanel("Lotka Volterra model"),
  sidebarLayout(
    sidebarPanel(
      sliderInput(inputId=names(s[1]),label=names(s[1]),min=0,max=10,step=0.01,value=as.numeric(s[1])),
      sliderInput(inputId=names(s[2]),label=names(s[2]),min=0,max=10,step=0.01,value=as.numeric(s[2])),
      sliderInput(inputId=names(p[1]),label=names(p[1]),min=0,max=p[1]*10,step=p[1]*0.01,value=p[1]),
      sliderInput(inputId=names(p[2]),label=names(p[2]),min=0,max=p[2]*10,step=p[2]*0.01,value=p[2]),
      sliderInput(inputId=names(p[3]),label=names(p[3]),min=0,max=p[3]*10,step=p[3]*0.01,value=p[3]),
      sliderInput(inputId=names(p[4]),label=names(p[4]),min=0,max=p[4]*10,step=p[4]*0.01,value=p[4]),
      sliderInput(inputId=names(p[5]),label=names(p[5]),min=0,max=p[5]*10,step=p[5]*0.01,value=p[5])
    ),

    mainPanel(
      h4("R '= r*R*(1 - R/K) - a*R*N"),
      h4("N '= c*a*R*N - delta*N"),
      plotOutput(outputId="grind"),
      radioButtons("radio", "Grind output:",c("Time plot"=0,"Nullclines"=1,"Trajectory"=2,"Portrait"=3,"Steady state"=4),selected=1,inline=TRUE),
      fluidRow(
      column(3,numericInput(inputId="tmax",label="Tmax",value=100)),
      column(3,numericInput(inputId="tstep",label="Tstep",value=0.1,step=0.1))),
      textOutput("log")
    )
  )
)

server <- function(input, output) {
  output$grind <- renderPlot({
    radiob <- input$radio
    for (i in names(s)) s[i] <- input[[i]]
    for (i in names(p)) p[i] <- input[[i]]
    tmax <- input$tmax
    tstep <- input$tstep
    xmax <- 1.1*max(p["K"],p["delta"]/(p["c"]*p["a"]),s["R"])
    ymax <- 1.1*max(p["r"]/p["a"],s["N"])
    output$log <- renderText("")
    if (radiob == 0) {
      f <- run(tmax=tmax,tstep=0.1,odes=model,state=s,parms=p)
      output$log <- renderText({paste("Ended in R =",round(f[1],5),"N =",round(f[2],5),sep=" ")})
    }
    else if (radiob <= 2) {
      plane(xmax=xmax,ymax=ymax,odes=model,state=s,parms=p,eps=-.001)
      if (radiob == 2) {
        f <- run(tmax=tmax,tstep=tstep,odes=model,state=s,parms=p,traject=TRUE)
        output$log <- renderText({paste("Ended in R =",round(f[1],5),"N =",round(f[2],5),sep=" ")})
      }
    }
    else if (radiob == 3) plane(xmax=xmax,ymax=ymax,tmax=tmax,tstep=tstep,odes=model,state=s,parms=p,eps=-.001,portrait=TRUE)
    else {
      plane(xmax=xmax,ymax=ymax,odes=model,state=s,parms=p,eps=-.001)
      f <- newton(state=s,parms=p,odes=model,plot=TRUE,positive=TRUE)
      if (is.null(f)) output$log <- renderText({"No convergence: start closer to a steady state"})
      else output$log <- renderText({paste("Converged into R =",round(f[1],5),"N =",round(f[2],5),sep=" ")})
    }
  })
}

shinyApp(ui = ui, server = server)

