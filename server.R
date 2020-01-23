shinyServer(function(input, output, session) {
  
  output$cluster_select <- renderUI(my_checkboxGroupInput('clusterID', h5('Cluster ID'),
                                                          choices = my_names,
                                                          selected=my_selected, 
                                                          colors=my_colors))
  
  
  output$h <- reactive(1000)
  
  # Logic for search inputs
  #snp_search <- function(){return(F)}#reactive(return(input$search_type == 'SNP ID'))
  observe({
    shinyjs::toggle(id = 'snp_search', condition = input$search_type == 'SNP ID')
    shinyjs::toggle(id = 'coord_chr', condition = input$search_type == 'Coordinates' & input$gene == 'Any')
    shinyjs::toggle(id = 'coord_min', condition = input$search_type == 'Coordinates' & input$gene == 'Any')
    shinyjs::toggle(id = 'coord_max', condition = input$search_type == 'Coordinates' & input$gene == 'Any')
    shinyjs::toggle(id = 'xrange', condition = input$search_type != 'Coordinates' | input$gene != 'Any')
  })
  
  
  
  # Update gene names
  geneNames <- reactive({
    
    switch(input$search_type,
           'Gene Name' = gene_names,
           'SNP ID' = {
             if(sum(input$snp_search %in% snp_table$SNP.ID) > 0){
               snp_gene(snp_id = input$snp_search)
             } else{
               'Input SNP ID'
             }
           },
           'Coordinates' = {
             if(is.numeric(input$coord_max) & is.numeric(input$coord_min)){
               if(input$coord_max - input$coord_min >= 0){
                 genes_xrange(xmin = input$coord_min, xmax = input$coord_max, chr = input$coord_chr)
               } else{
                 c()
               }
             } else{
               c()
             }
           }
    )
  })
  
  geneNamesUpdate <- eventReactive(input$plot_button, {
    geneNames()
  })
  
  snp_xrange <- reactive({
    snp_xrange_temp <- filter(snp_table, SNP.ID == input$snp_search)
    if(nrow(snp_xrange_temp) > 0){
      snp_xrange_temp
    } else{
      NULL
    }
  })
  
  observe({
    if(input$search_type != 'Gene Name'){
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
    getGtfCoords(input$gene)
  })
  
  
  observe({
    if(input$search_type == 'Coordinates' & input$gene == 'Any'){
      updateSliderInput(session, 'xrange',
                        min = input$coord_min - 3E5 - SLOP,
                        max = input$coord_max + 3E5 + SLOP,
                        value = c(input$coord_min, input$coord_max))
    } else if(input$search_type == 'SNP ID' & input$gene == 'Any'){
      updateSliderInput(session, 'xrange',
                        min = snp_xrange()$start - 3E5 - SLOP,
                        max = snp_xrange()$end + 3E5 + SLOP,
                        value = c(snp_xrange()$start - SLOP, snp_xrange()$end + SLOP)
                        )
    } else{
      updateSliderInput(session, 'xrange',
                        min = gvizCoords()[[2]] - 3E5 - SLOP,
                        max = gvizCoords()[[3]] + 3E5 + SLOP,
                        value = c(gvizCoords()[[2]] - SLOP,gvizCoords()[[3]] + SLOP))
    }
    
  }, priority = 100)
  
  tfMotifFilter <- reactive({
    filtered = motifs[seqnames(motifs) == chrom() & 
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
  chrom <- reactiveVal(unlist(getGtfCoords('NEUROG2')[1]))
  
  
  observeEvent(input$plot_button, {
    gene(input$gene)
    ymax(input$ymax)
    
    xrange({
      if(input$gene == 'Any' & input$search_type == 'Coordinates'){
        c(input$coord_min, input$coord_max)
      } else{
        input$xrange
      }
    })
    
    chrom({
      if(input$gene == 'Any' & input$search_type == 'Coordinates'){
        input$coord_chr
      } else if(input$gene == 'Any' & !is.null(snp_xrange) & input$search_type == 'SNP ID'){
         as.character(snp_xrange()$chrom)
      } else{
        unlist(getGtfCoords(input$gene)[[1]])
      }
    })
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
  
  # Removing loading message
  delay(0,hide(id = "loading-content", anim = TRUE, animType = "fade", time = 0.5))
  
  gvizPlot <- reactive({
    if(gene() == 'Any'){
      plotGenomeView(
        gene.symbol = geneNamesUpdate(),
        coverage.list = coverage.list,
        ylims = c(0, ymax()),
        chr = chrom(),
        beg = xrange()[1],
        END = xrange()[2],
        transcriptID = values_d$transcriptID(),
        corCut = values_d$cor_cut(),
        pval_cut = values_d$pval_cut(),
        cluster_id = values_d$clusterID(),
        motifs_list = values_d$tf_motifs(),
        selected_rows = values_d$selectedRows(),
        greater_less = values_d$greater_less(),
        snp_xrange = snp_xrange()
      )
    } else{
      plotGenomeView(
        gene.symbol = gene(),
        coverage.list = coverage.list,
        ylims = c(0, ymax()),
        chr = chrom(),
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
    }
    
  })
  
  
  output$gviz <- renderPlot({
    gvizPlot()
  })
  
  outputOptions(output, suspendWhenHidden = FALSE)
  outputOptions(output, 'gviz', suspendWhenHidden = FALSE)
  outputOptions(output, 'rnaExpression.ui', suspendWhenHidden = FALSE)
  
  
  
  # Peaks table
  cor.gr_table <- reactive({
    if(gene() == 'Any'){
      corFilter(cor_table = cor.gr, greater_less = values_d$greater_less(), cor_cut = values_d$cor_cut(),
                transcript_ID = values_d$transcriptID(), cluster_id = values_d$clusterID(), pval_cut = values_d$pval_cut(),
                beg = xrange()[1], END = xrange()[2], chr = chrom(), gene.name = geneNamesUpdate()) %>% 
        dplyr::as_tibble() %>% dplyr::select(colors, everything())
    } else{
      corFilter(cor_table = cor.gr, greater_less = values_d$greater_less(), cor_cut = values_d$cor_cut(),
                transcript_ID = values_d$transcriptID(), cluster_id = values_d$clusterID(), pval_cut = values_d$pval_cut(),
                beg = xrange()[1], END = xrange()[2], chr = chrom(), gene.name = gene()) %>% 
        dplyr::as_tibble() %>% dplyr::select(colors, everything())
    }
    
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
                     chr = chrom()
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
