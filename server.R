shinyServer(function(input, output, session) {
  
  geneNames <- reactive({
    type <- input$accessionType
    #isolate(input$genes %>% strsplit(' ') %>% unlist() %>% tolower())
    unique(data[,type])
  })
  
  observe({
    updateSelectInput(session, 'gene',
                      choices = geneNames()
    )})
  
  output$rnaExpression <- renderPlot({
    if(length(input$gene > 0)){
      plotrna(input$gene,data) 
    } else{
      NULL
    }
    
  })
  
})
