shinyServer(function(input, output, session) {
  
  output$cluster_select <- renderUI(my_checkboxGroupInput('clusterID', h5('Cluster ID'),
                                                          choices = my_names,
                                                          selected=my_selected, 
                                                          colors=my_colors))
  
  search_types <- reactive(return(input$search_type == 'SNP ID'))
  
  output$h <- reactive(1000)
  
  # Logic to control if SNP ID search should be displayed or not
  observe({
    shinyjs::toggle(id = 'snp_search', condition = search_types())
  })
  

  # Update gene names
  geneNames <- reactive({
    if(search_types()){
      if(sum(input$snp_search %in% snp_table$SNP.ID) > 0){
        snp_gene(snp_id = input$snp_search)
      } else{
        'Genes near selected SNP'
      }
    } else{
      gene_names
    }
  })
  
  observe({
    if(search_types()){
      updatePickerInput(
        session,
        inputId = 'gene',
        choices = geneNames(),
        selected = geneNames()[1]
      )
    } else{
      updatePickerInput(
        session,
        inputId = 'gene',
        choices = geneNames(),
        selected = 'NEUROG2'
      )
    }
    
  })
  
  
  transcriptID <- eventReactive(input$plot_button, {
    data %>% filter(gene.symbol == input$gene) %>%
      dplyr::select(transcript_id) %>% unlist() %>% as.vector() %>% unique()
  }, ignoreNULL = FALSE)
  
  observeEvent(input$plot_button, {
    updateSelectInput(session, 'transcriptID', choices = c('Any',transcriptID()))
  }, ignoreNULL = FALSE, priority = 1)
  
  
  rnaPlot <- reactive({
    plotrna(gene(), data)
  })
  
  output$rnaExpression <- renderPlot({
    rnaPlot()
  }, height = function(){
    if(length(rnaPlot()$data) > 0){
      ceiling(length(unique(rnaPlot()$data$transcript_id))/4) * 350
    } else{
      350
    }
  }
  
  
  )
  
  output$rnaExpression.ui <- renderUI({
    
    plotOutput('rnaExpression') %>% withSpinner(type = 6)
    
    
  })
  
  
  #outputOptions(output, 'dims', suspendWhenHidden = FALSE)
  
  gvizCoords <- reactive({
    getGtfCoords(input$gene, ENSEMBL_hg38_local_fromGTF)
  })
  
  
  observe({
    updateSliderInput(session, 'xrange',
                      min = gvizCoords()[[2]] - 3E5 - SLOP,
                      max = gvizCoords()[[3]] + 3E5 + SLOP,
                      value = c(gvizCoords()[[2]] - SLOP,gvizCoords()[[3]] + SLOP))
  }, priority = 9)
  
  tfMotifFilter <- reactive({
    filtered = motifs[seqnames(motifs) == chr() & 
                        start(motifs) >= xrange()[1] & 
                        end(motifs) <= xrange()[2]]
    
    tfPresent <- lapply(motifs@colData@rownames, function(x) {
      filtered[assay(filtered)[,x]] %>% length() > 0
    }) %>% unlist()
    
    names <- motifs@colData@rownames[tfPresent]
    names[order(names)]
  })
  
  observe({
    updatePickerInput(
      session,
      inputId = 'tf_motifs',
      selected = NULL,
      choices = tfMotifFilter()
    )
  })
  
  
  # Define and update input variables  
  gene <- reactiveVal('NEUROG2')
  ymax <- reactiveVal(200)
  xrange <- reactiveVal(c(112463516, 112566180))
  chr <- reactiveVal(unlist(getGtfCoords('NEUROG2', ENSEMBL_hg38_local_fromGTF)[1]))
  
  observeEvent(input$plot_button, {
    gene(input$gene)
    ymax(input$ymax)
    xrange(input$xrange)
    chr(unlist(getGtfCoords(input$gene, ENSEMBL_hg38_local_fromGTF)[[1]]))
  })
  
  
  # Throttle response to dynamic inputs by using
  # the debounce function
  values_d <- reactiveValues()
  
  observe({
    values_d$transcriptID <- debounce(function(){
      if(input$transcriptID == 'Any'){
        return(transcriptID())
      } else{
        return(input$transcriptID)
      }
    },0)
    values_d$clusterID <- debounce(function(){input$clusterID},2500)
    values_d$cor_cut <- debounce(function(){input$cor_cut},2500)
    values_d$pval_cut <- debounce(function(){input$pval_cut},2500)
    values_d$tf_motifs <- debounce(function(){input$tf_motifs},2500)
    values_d$selectedRows <- debounce(function(){input$peaks_table_rows_selected}, 2500)
    values_d$greater_less <- debounce(function(){input$greater_less}, 1000)
  })
  
  gvizPlot <- reactive({
    
    plotGenomeView(
      gene.symbol = gene(),
      coverage.list = coverage.list,
      ylims = c(0, ymax()),
      coords = NULL, #gvizCoords(),
      chr = chr(),
      beg = xrange()[1],
      END = xrange()[2],
      transcriptID = values_d$transcriptID(),
      corCut = values_d$cor_cut(),
      pval_cut = values_d$pval_cut(),
      cluster_id = values_d$clusterID(),
      motifs_list = values_d$tf_motifs(),
      selected_rows = values_d$selectedRows(),
      greater_less = values_d$greater_less()
    )
  })
  
  
  output$gviz <- renderPlot({
    gvizPlot()
  })
  
  outputOptions(output, suspendWhenHidden = FALSE)
  outputOptions(output, 'gviz', suspendWhenHidden = FALSE)
  outputOptions(output, 'rnaExpression.ui', suspendWhenHidden = FALSE)
  
  # Plot clusters and motifs tracks
  # gvizPlotClust <- reactive({
  # 
  #   plot_clust_motif(
  #     gene.symbol = gene(),
  #     coverage.list = coverage.list,
  #     ylims = c(0, ymax()),
  #     coords = NULL, #gvizCoords(),
  #     chr = chr(),
  #     beg = xrange()[1],
  #     END = xrange()[2],
  #     transcriptID = input$transcriptID,
  #     corCut = values_d$cor_cut(),
  #     pval_cut = values_d$pval_cut(),
  #     cluster_id = values_d$clusterID(),
  #     motifs_list = values_d$tf_motifs()
  #   )
  # })
  # 
  # output$gvizClust <- renderPlot({
  #   gvizPlotClust()
  # })
  
  # Peaks table
  cor.gr_table <- reactive({
    
    corFilter(cor_table = cor.gr, greater_less = values_d$greater_less(), cor_cut = values_d$cor_cut(),
              transcript_ID = values_d$transcriptID(), cluster_id = values_d$clusterID(), pval_cut = values_d$pval_cut(),
              beg = xrange()[1], END = xrange()[2]) %>% dplyr::as_tibble() %>% dplyr::select(colors, everything())
  })
  
  
  output$peaks_table <- DT::renderDataTable({ 
    dat <- datatable(cor.gr_table(), 
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
                     )
    ) %>% formatStyle('colors',  
                      color = styleEqual(c('#ABDDA4'), c('#ABDDA4')), 
                      backgroundColor = styleEqual(c('#ABDDA4'), c('#ABDDA4'))
    )
    return(dat)
  }, class = "display"
  )
  
  # TF Motifs
  output$selected_tf_motifs <- renderText({
    if(!is.null(input$tf_motifs)){
      paste0('Selected TFs: ', paste0(input$tf_motifs, collapse = ', '))
    } else{
      NULL
    }
  })
  
  tf_legend <- reactive({
    if(!is.null(input$tf_motifs)){
      
      motifs_legend_df <- tibble(`TF Motifs` = as.factor(input$tf_motifs),
                                 x = 1:length(input$tf_motifs)) %>%
        mutate(`TF Motifs` = fct_reorder(`TF Motifs`, 1:length(`TF Motifs`)))
      
      p <- ggplot(motifs_legend_df, aes(x = x, fill = `TF Motifs`)) + 
        geom_bar() + scale_fill_jco() + theme(legend.position="bottom")
      
      leg <- get_legend(p)
      as_ggplot(leg)
      
    } else{
      NULL
    }
  })
  
  output$tf_legend <- renderPlot({
    tf_legend()
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
      params <- list(gene = gene(),
                     ymax = ymax(),
                     xrange = xrange(),
                     gvizPlot = gvizPlot(),
                     tf_legend = tf_legend(),
                     rnaPlot = rnaPlot(),
                     transcriptID = input$transcriptID,
                     chr = chr()
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
