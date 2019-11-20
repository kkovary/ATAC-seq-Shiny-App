shinyServer(function(input, output, session) {
  geneNames <- reactive({
    #type <- input$accessionType
    type <- 'gene.symbol'
    unique(data[, type])
  })
  
  transcriptID <- reactive({
    data %>% filter(gene.symbol == input$gene) %>%
      dplyr::select(transcript_id) %>% unlist() %>% as.vector() %>% unique()
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
  
  rnaPlot <- eventReactive(input$plot_button, {
    plotrna(input$gene, data)
  })
  
  output$rnaExpression <- renderPlot({
    if(input$gene != '') {
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
    if(input$gene !=''){
      updateSliderInput(session, 'xrange',
                        min = gvizCoords()[[2]] - 2E5 - SLOP,
                        max = gvizCoords()[[3]] + 2E5 + SLOP,
                        value = c(gvizCoords()[[2]] - SLOP,gvizCoords()[[3]] + SLOP))
    }
  })
  
  
  gvizPlot <- eventReactive(input$plot_button, {
    
    plotGenomeView(
      gene.symbol = input$gene,
      coverage.list = coverage.list,
      ylims = c(0, input$ymax),
      coords = gvizCoords(),
      chr = gvizCoords()[[1]],
      beg = input$xrange[1],
      END = as.numeric(input$xrange[2]),
      transcriptID = input$transcriptID,
      corCut = input$cor_cut,
      pval_cut = input$pval_cut,
      cluster_id = input$clusterID
    )
  })
  
  # output$gviz <- renderPlot({
  #   if(!is.null(input$gene)) {
  #     gvizPlot()
  #   } else{
  #     NULL
  #   }
  #   
  # })
  output$gviz <- renderPlot({
    gvizPlot()
  })
  
  # Peaks table
  cor.gr_table <- eventReactive(input$plot_button, {
    if(is.na(input$cor_cut)){
      cor.gr.subset <- cor.gr[(elementMetadata(cor.gr)$transcript_id == input$transcriptID) &
                                (elementMetadata(cor.gr)$cluster.name %in% input$clusterID) &
                                (elementMetadata(cor.gr)$vs.null.p.value <= input$pval_cut)]
    } else{
      cor.gr.subset <- cor.gr[(elementMetadata(cor.gr)$transcript_id == input$transcriptID) &
                                (elementMetadata(cor.gr)$estimate >= input$cor_cut) & 
                                (elementMetadata(cor.gr)$cluster.name %in% input$clusterID) &
                                (elementMetadata(cor.gr)$vs.null.p.value <= input$pval_cut)]
    }
    dplyr::as_tibble(cor.gr.subset)
  })
  
  output$peaks_table <- renderDataTable(
    cor.gr_table(),
    extensions = 'Buttons',
    options = list(
      paging = TRUE,
      searching = TRUE,
      fixedColumns = TRUE,
      autoWidth = TRUE,
      ordering = TRUE,
      dom = 'Bfrtip',
      buttons = c('csv', 'excel'),
      pageLength = nrow(cor.gr_table())
    ),
    class = "display"
  )
  
  
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