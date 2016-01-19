shinyServer(
  function(input, output) {
  
  # If new data is requested
  gendata = reactive({
    gen.gts(AR1(phi=input$phi, sigma2=input$sigma2))
  })
  
  # Cache the model object
  gmwm_results = reactive({
    gmwm(AR1(), data=gendata(), model.type="imu")
  })
  
  # Handle model text
  output$modelText = renderPrint({
    print(gmwm_results())
  })
  
  # Render the plot
  output$modelPlot = renderPlot({
    autoplot(gmwm_results())
  })
})