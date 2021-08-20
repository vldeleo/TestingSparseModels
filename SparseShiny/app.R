#
# Run the application by clicking the 'Run App' button above.
# 

library(shiny)
library(monomvn)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Competing Sparse Models"),

    # User Specifications for Sim Data
        # Number of observations
        # Number of Predictors
        # Sparsity
    # User Specifications for Models
    
    # Output Plots
    
    # Sidebar with a slider input for number of observations 
#    sidebarLayout(
#        sidebarPanel(
            sliderInput(inputId = "obs",
                        "Number of Observations:",
                        min = 50,
                        max = 500,
                        value = 100),
            sliderInput(inputId = "preds",
                        "Number of Predictors:",
                        min = 5,
                        max = 500,
                        value = 50),
            radioButtons(inputId = "sparsity", # this might be better as a drop down box 
                         "Sparse Architecture",
                         choices = c("Dense", "0.33", "0.1", "0.02")), #),
        # Choose models   
            checkboxGroupInput(inputId = "models",
                          "Models to Compare:",
                          c("Ordinary Least Squares" = "ols",
                            "Partial Least Squares" ="pls",
                            "Lasso" = "lasso", 
                            "Ridge" = "ridge"
                            #, "blasso", "bridge", "bhs", "susie"
                            ),
                          selected = c("ols")
              ),
            
            submitButton("Update"), # button for refreshing calculations
            
#        ),

        # Show a output
#        mainPanel(
            plotOutput("histPlot"), #plot of the generated distribution
#            tableOutput("X") # this is a bit much to print to screen///
           plotOutput("coefPlot"), ### true vs model coef
            textOutput("testvar")
#           plotOutput("predPlot") ### true vs model predictions
#           plotOutput("timePlot") ### model diagnostics
#        )
#    )
)



# Define server logic to simulate and plot data
server <- function(input, output) {
    
    # Simulate Data
  n <- reactive({input$obs})
  p <- reactive({input$preds})
  
  Xfull <- reactive({
      n <- as.numeric(n())
      p <- as.numeric(p())
      return(matrix(rnorm((n)*p, 0, 1),nrow=(n),ncol=p))
    })
    
    output$X <- renderTable({
        Xfull()
#        Xfull_temp <- as.data.frame(Xfull())
#        X <- Xfull_temp[1:n,]
#        X # I'm having trouble evaluating and then setting aside a subset for prediction
      })
    
    beta <- reactive({
      n <- as.numeric(n())
      p <- as.numeric(p())
      beta_temp <- (rgamma(p, 0.02, 0.1) * sample(c(-1, 1), p, replace = T))
      return(beta_temp)
    })
    
    yfull <- reactive({
      n <- as.numeric(n())
      p <- as.numeric(p())
      beta <- as.numeric(beta())
      if (input$sparsity != "Dense"){
        sparsity <- as.numeric(input$sparsity)
        beta[sample(c(1:p), round(p*(1-sparsity),0), replace = F)] <- 0 # all but the fraction of predictors that we've decided should be true effects are set to 0
      }
      yfull_temp <- Xfull() %*% beta + rnorm(n) 
      y <- yfull_temp[1:n,] # this could be my smaller, "sampled" population, but in this case it's the same size as Xfull
      return(y)
    })
    
    
    output$histPlot <- renderPlot({
        hist(yfull(), col = 'darkgray', border = 'white', main = "Response Values", xlab = "")})
  
    
    # Fit models 
    mods <- reactive({input$models})
    
    # Plot models
    output$coefPlot <- renderPlot({
      modlist <- mods()
      if ("ols" %in% modlist){
        pt1 <- proc.time() 
        olsres <- regress(as.matrix(Xfull()), as.matrix(yfull()), method = c("lsr"), p = 0)
        olsTime <- as.numeric(proc.time() - pt1)[3]
#      olsResid <- c(olsres$b[1] + as.matrix(X) %*% olsres$b[-1]) - y
#      olsRMSE <- sqrt(mean((olsResid)^2))
      olsplot <- plot(as.numeric(beta()), as.numeric(olsres$b[-1]), xlab = "True Beta", ylab = "Model Beta", main = "OLS")
      }
      if ("pls" %in% modlist){
        pt1 <- proc.time() 
        plsres <- regress(as.matrix(Xfull()), as.matrix(yfull()), method = c("plsr"), p = 0)
        plsTime <- as.numeric(proc.time() - pt1)[3]
        #      plsResid <- c(plsres$b[1] + as.matrix(X) %*% plsres$b[-1]) - y
        #      plsRMSE <- sqrt(mean((plsResid)^2))
        plsplot <- plot(as.numeric(beta()), as.numeric(plsres$b[-1]), xlab = "True Beta", ylab = "Model Beta", main = "PLS")
        
      }
#      if ("lasso" %in% mods()){
        
#      }
#      if ("ridge" %in% mods()){
        
#      }
      if(exists("olsplot")){olsplot}
      if(exists("plsplot")){plsplot}
    })
        
#    output$effectSize <- renderPlot({
        
#    })
        
#    output$predAccuracy <- renderPlot({
        
#    })
        
    
#    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
#        x    <- faithful[, 2]
#        bins <- seq(min(x), max(x), length.out = input$bins + 1)
        # draw the histogram with the specified number of bins
#        hist(x, breaks = bins, col = 'darkgray', border = 'white')
#    })
}



# Run the application 
shinyApp(ui = ui, server = server)
