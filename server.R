shinyServer(function(input, output, session) {
  geneNames <- reactive({
    type <- input$accessionType
    unique(data[, type])
  })
  
  transcriptID <- reactive({
    data %>% filter(gene.symbol == input$gene) %>%
      dplyr::select(transcript_id) %>% unlist() %>% as.vector() %>% unique()
  })
  
  observe({
    updateSelectInput(session, 'gene',
                      choices = geneNames())
  })
  
  observe({
    if(is.null(input$gene)){
      updateSelectInput(session, 'transcriptID',
                        choices = '')
    } else{
      updateSelectInput(session, 'transcriptID',
                        choices = transcriptID())
    }
  })
  
  rnaPlot <- reactive({
    plotrna(input$gene, data)
  })
  
  output$rnaExpression <- renderPlot({
    if (length(input$gene > 0)) {
      rnaPlot()
    } else{
      NULL
    }
    
  })
  
  
  plotCount <- reactive({
    if (length(input$gene > 0)) {
      dim = plotrna(input$gene, data)$data$transcript_id
      ceiling(length(unique(dim)) / 4)
    } else{
      1
    }
  })
  
  plotHeight <- reactive(350 * plotCount())
  
  output$rnaExpression.ui <- renderUI({
    plotOutput('rnaExpression', height = plotHeight())
  })
  
  
  #outputOptions(output, 'dims', suspendWhenHidden = FALSE)
  
  gvizCoords <- reactive({
    if(length(input$gene) > 0){
      getGvizCoords(input$gene)
    }
  })
  
  observe({
    
    updateSliderInput(session, 'xrange',
                      min = gvizCoords()[[2]] - 1E5 - SLOP,
                      max = gvizCoords()[[3]] + 1E5 + SLOP,
                      value = c(gvizCoords()[[2]] - SLOP,gvizCoords()[[3]] + SLOP))
    
  })
  
  
  gvizPlot <- reactive({
    
    plotGenomeView(
      gene.symbol = input$gene,
      coverage.list = coverage.list,
      ylims = c(0, input$ymax),
      coords = gvizCoords(),
      chr = gvizCoords()[[1]],
      beg = input$xrange[1],
      END = as.numeric(input$xrange[2]),
      transcriptID = input$transcriptID,
      cor_cut = input$cor_cut,
      cluster_id = input$clusterID
    )
  })
  
  output$gviz <- renderPlot({
    if(sum(input$xrange) > 0) {
      gvizPlot()
    } else{
      NULL
    }
    
  })
  
  
  ### Download Report
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "report.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(gene = input$gene,
                     ymax = input$ymax,
                     gvizPlot = gvizPlot(),
                     rnaPlot = rnaPlot()
      )
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(
        tempReport,
        output_file = file,
        params = params,
        envir = new.env(parent = globalenv())
      )
    }
  )
  
  
})
