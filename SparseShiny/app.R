#
# Run the application by clicking the 'Run App' button above.
# 

library(shiny)
library(monomvn)
library(scales)


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
                        max = 500, # app will crash if there are too many predictors
                        value = 50),
            selectInput(inputId = "sparsity", # can also present this as radio buttons 
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
#           tableOutput("X") # this is a bit much to print to screen///
            plotOutput("coefPlot"), ### true vs model coef
            plotOutput("predAccuracy") ### true vs model predictions
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
      nbig <- signif(3*n/2, 0)
      return(matrix(rnorm((nbig)*p, 0, 1),nrow=(nbig),ncol=p))
    })
    
    output$X <- renderTable({
        Xbig <- Xfull()
        Xsmall <- as.data.frame(Xbig)[1:n,]
        return(Xsmall)
#        Xfull_temp <- as.data.frame(Xfull())
#        X <- Xfull_temp[1:n,]
#        X # I'm having trouble evaluating and then setting aside a subset for prediction
      })
    
    beta <- reactive({
      n <- as.numeric(n())
      p <- as.numeric(p())
      beta_temp <- (rgamma(p, 0.02, 0.1) * sample(c(-1, 1), p, replace = T))
      if (input$sparsity != "Dense"){
        sparsity <- as.numeric(input$sparsity)
        beta_temp[sample(c(1:p), round(p*(1-sparsity),0), replace = F)] <- 0 # all but the fraction of predictors that we've decided should be true effects are set to 0
      }
      return(beta_temp)
    })
    
    yfull <- reactive({
      n <- as.numeric(n())
      p <- as.numeric(p())
      beta <- as.numeric(beta())
      yfull_temp <- Xfull() %*% beta + rnorm(signif(3*n/2, 0)) 
      y <- yfull_temp # this could be my smaller, "sampled" population, but in this case it's the same size as Xfull
      return(y)
    })
    
    
    output$histPlot <- renderPlot({
      n <- as.numeric(n())
      ysamp <- yfull()[1:n]
        hist(ysamp, col = 'darkgray', border = 'white', main = "Response Values", xlab = "")})  # eventually, I think I want to make visualizing actual data optional
  
    
    # Fit models 
    mods <- reactive({input$models})
    
    olsreact <- reactive({
      n <- as.numeric(n())
      Xsamp <- as.data.frame(Xfull())[1:n,]
      ysamp <- yfull()[1:n]
    pt1 <- proc.time() 
    olsres <- regress(as.matrix(Xsamp), as.matrix(ysamp), method = c("lsr"), p = 0)
    olsres$olsTime <- as.numeric(proc.time() - pt1)[3]
    olsres$olsResid <- c(olsres$b[1] + as.matrix(Xsamp) %*% olsres$b[-1]) - ysamp
    olsres$olsRMSE <- signif(sqrt(mean((olsres$olsResid)^2)), 4)
    olsres
    })
    
    plsreact <- reactive({
      modlist <- mods()
      n <- as.numeric(n())
      Xsamp <- as.data.frame(Xfull())[1:n,]
      ysamp <- yfull()[1:n]
      if ("pls" %in% modlist){
      pt1 <- proc.time() 
      plsres <- regress(as.matrix(Xsamp), as.matrix(ysamp), method = c("plsr"), p = 0)
      plsres$plsTime <- as.numeric(proc.time() - pt1)[3]
      plsres$plsResid <- c(plsres$b[1] + as.matrix(Xsamp) %*% plsres$b[-1]) - ysamp
      plsres$plsRMSE <- signif(sqrt(mean((plsres$plsResid)^2)), 4)
      }
      plsres
      })
    
    lasreact <- reactive({
      modlist <- mods()
      n <- as.numeric(n())
      Xsamp <- as.data.frame(Xfull())[1:n,]
      ysamp <- yfull()[1:n]
      if ("lasso" %in% modlist){
        pt1 <- proc.time() 
        lasres <- regress(as.matrix(Xsamp), as.matrix(ysamp), method = c("lasso"), p = 0)
        lasres$lasTime <- as.numeric(proc.time() - pt1)[3]
        lasres$lasResid <- c(lasres$b[1] + as.matrix(Xsamp) %*% lasres$b[-1]) - ysamp
        lasres$lasRMSE <- signif(sqrt(mean((lasres$lasResid)^2)), 4)
      }
      lasres
    })
    
    ridreact <- reactive({
      modlist <- mods()
      n <- as.numeric(n())
      Xsamp <- as.data.frame(Xfull())[1:n,]
      ysamp <- yfull()[1:n]
      if ("ridge" %in% modlist){
        pt1 <- proc.time() 
        ridres <- regress(as.matrix(Xsamp), as.matrix(ysamp), method = c("ridge"), p = 0)
        ridres$ridTime <- as.numeric(proc.time() - pt1)[3]
        ridres$ridResid <- c(ridres$b[1] + as.matrix(Xsamp) %*% ridres$b[-1]) - ysamp
        ridres$ridRMSE <- signif(sqrt(mean((ridres$ridResid)^2)), 4)
      }
      ridres
    })
    
    # Plot models
    output$coefPlot <- renderPlot({
      modlist <- mods()
#      if ("ols" %in% modlist){ #I should probably have this run whether or not the user wants it...
      olsres <- olsreact()
      olsplot <- c(plot(as.numeric(beta()), as.numeric(olsres$b[-1]), xlab = "True Beta", ylab = "Model Beta", main = "OLS"), legend("bottom", legend = paste("Model RMSE = ", signif(olsres$olsRMSE, 3), sep = ""), bty = "n"))
#      }
      if ("pls" %in% modlist){
        plsres <- plsreact()
        plsplot <- c(plot(as.numeric(beta()), as.numeric(plsres$b[-1]), xlab = "True Beta", ylab = "Model Beta", main = "PLS"), legend("bottom", legend = paste("Model RMSE = ", signif(plsres$plsRMSE, 3), sep = ""), bty = "n"))
      }
      if ("lasso" %in% modlist){
        lasres <- lasreact()
        lassoplot <- c(plot(as.numeric(beta()), as.numeric(lasres$b[-1]), xlab = "True Beta", ylab = "Model Beta", main = "Lasso"),  legend("bottom", legend = paste("Model RMSE = ", signif(lasres$lasRMSE, 3), sep = ""), bty = "n"))
      }
      if ("ridge" %in% modlist){
        ridres <- ridreact()
        ridgeplot <- c(plot(as.numeric(beta()), as.numeric(ridres$b[-1]), xlab = "True Beta", ylab = "Model Beta", main = "Ridge"),  legend("bottom", legend = paste("Model RMSE = ", signif(ridres$ridRMSE, 3), sep = ""), bty = "n"))
      }
      
#      togPlot <- c(
        plot(as.numeric(beta()), as.numeric(olsres$b[-1]), xlab = "True Beta", ylab = "Model Beta", col = "white")
        abline(0, 1, col = "gray", lwd = 3)
        clip(min(beta()) - .1, max(beta()) + .1, par("usr")[3], par("usr")[4]) # I would like to fix this so that regression lines get clipped but points do not
          # at the moment, it cuts off my RMSE legend as well
        if(exists("olsres")){
          points(as.numeric(beta()), as.numeric(olsres$b[-1]), col = alpha("black", 0.75), lwd = 3)
          abline(lm(olsres$b[-1] ~ beta()), col = "black")
          }
        if(exists("plsres")){
          points(as.numeric(beta()), as.numeric(plsres$b[-1]), col = alpha("red", 0.75), lwd = 3)
          abline(lm(plsres$b[-1] ~ beta()), col = "red")
          }
        if(exists("lasres")){
          points(as.numeric(beta()), as.numeric(lasres$b[-1]), col = alpha("blue", 0.75), lwd = 3)
          abline(lm(lasres$b[-1] ~ beta()), col = "blue")
          }
        if(exists("ridres")){
          points(as.numeric(beta()), as.numeric(ridres$b[-1]), col = alpha("darkgoldenrod1", 0.75), lwd = 3)
          abline(lm(ridres$b[-1] ~ beta()), col = "darkgoldenrod1")
        }
        legend("bottomright", bty = "n", legend = c(if(exists("olsres")){paste("OLS RMSE =", olsres$olsRMSE, sep = "")}, if(exists("plsres")){paste("PLS RMSE =", plsres$plsRMSE, sep = "")}, if(exists("lasres")){paste("Lasso RMSE =", lasres$lasRMSE, sep = "")}, if(exists("ridres")){paste("Ridge RMSE =", ridres$ridRMSE, sep = "")}), text.col = c(if(exists("olsres")){"black"}, if(exists("plsres")){"red"}, if(exists("lasres")){"blue"}, if(exists("ridres")){"darkgoldenrod1"}))
    })
     

      
    
#    output$effectSize <- renderPlot({
        
#    })
        
    output$predAccuracy <- renderPlot({ 
     modlist <- mods()
     n <- as.numeric(n())
     yfull <- yfull()
     Xfull <- Xfull()
     XOOS <- Xfull[(n+1):nrow(Xfull),]
     OOSvals <- as.numeric(yfull[(n+1):length(yfull)])
     # now to get the model predictions
     olsres <- olsreact()
     olspred <- olsres$b[1] + as.matrix(XOOS) %*% olsres$b[-1]
     plot(OOSvals, olspred, xlab = "True Unobserved Values", ylab = "Model Predicted Values", col = "white")
     abline(0, 1, col = "gray", lwd = 3)
     if("ols" %in% modlist){
       points(OOSvals, olspred, col = alpha("black", 0.75), lwd = 3)
       abline(lm(olspred ~ OOSvals), col = "black")
     }
     if("pls" %in% modlist){
       plsres <- plsreact()
       plspred <- plsres$b[1] + as.matrix(XOOS) %*% plsres$b[-1]
       points(OOSvals, plspred, col = alpha("red", 0.75), lwd = 3)
       abline(lm(plspred ~ OOSvals), col = "red")
     }
     if("lasso" %in% modlist){
       lasres <- lasreact()
       laspred <- lasres$b[1] + as.matrix(XOOS) %*% lasres$b[-1]
       points(OOSvals, laspred, col = alpha("blue", 0.75), lwd = 3)
       abline(lm(laspred ~ OOSvals), col = "blue")
     }
     if("ridge" %in% modlist){
       ridres <- ridreact()
       ridpred <- ridres$b[1] + as.matrix(XOOS) %*% ridres$b[-1]
       points(OOSvals, ridpred, col = alpha("darkgoldenrod1", 0.75), lwd = 3)
       abline(lm(ridpred ~ OOSvals), col = "darkgoldenrod1")
     }


    
        
    })
        
    

}



# Run the application 
shinyApp(ui = ui, server = server)
