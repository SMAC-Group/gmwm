shinyUI(fluidPage(
  
  titlePanel("GMWM Model Specification"),
  
  sidebarLayout(
    sidebarPanel(
      h3("AR1"),
      sliderInput("phi", "Phi:",
                  min = 0.0001, max = 0.99999, value = 0.5, step = 0.001),
      sliderInput("sigma2", "Sigma^2:",
                  min = 0.0001, max = 20, value = 0.5, step = 0.001),
      
      sliderInput("n", "Sample Size:",
                  min = 10, max = 1000000, value = 1000, step = 1),
      
      submitButton("Get New Data"),
      br(),
      h5("Created by:"),
      a("SMAC-Group", 
             href="http://www.smac-group.com"),
      br(),
      h5("Source available at:"),
      a("GitHub Repository:", 
             href="https://github.com/smac-group/gmwm")
    ),
    
    mainPanel(
      plotOutput("modelPlot"),
      verbatimTextOutput("modelText")
    )
  )
))