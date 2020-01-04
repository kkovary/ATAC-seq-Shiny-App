shinyServer(function(input, output, session) {
  
  output$cluster_select <- renderUI(my_checkboxGroupInput('clusterID', h5('Cluster ID'),
                                                   choices = my_names,
                                                   selected=my_selected, 
                                                   colors=my_colors))
  
  # geneNames <- reactive({
  #   #type <- input$accessionType
  #   type <- 'gene.symbol'
  #   names <- data[, type] %>% unique() %>% unlist() %>% as.vector()
  #   names[order(names)]
  # })
  
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
  
  rnaPlot <- reactive({
    plotrna(input$gene, data)
  })
  
  output$rnaExpression <- renderPlot({
    
      rnaPlot()

  })
  
  
  plotCount <- reactive({
    if (length(input$gene > 0)) {
      dim = plotrna(input$gene, data)$data$transcript_id
      ceiling(length(unique(dim)) / 4)
    }
  })
  
  plotHeight <- reactive({
    
      350 * plotCount()
    
    })
  
  output$rnaExpression.ui <- renderUI({
    
      plotOutput('rnaExpression', height = plotHeight()) %>% withSpinner(type = 6)
    
    
  })
  
  
  #outputOptions(output, 'dims', suspendWhenHidden = FALSE)
  
  gvizCoords <- reactive({
    if(length(input$gene) > 0){
      # Switch to getBM
      getGvizCoords(input$gene)
      
    }
  })
  

  observe({
    if(input$gene !=''){
      
      updateSliderInput(session, 'xrange',
                        min = gvizCoords()[[2]] - 3E5 - SLOP,
                        max = gvizCoords()[[3]] + 3E5 + SLOP,
                        value = c(gvizCoords()[[2]] - SLOP,gvizCoords()[[3]] + SLOP))
    }
  })
  
  tfMotifFilter <- reactive({
      filtered = motifs[seqnames(motifs) == gvizCoords()[[1]] & 
                          start(motifs) >= input$xrange[1] & 
                          end(motifs) <= input$xrange[2]]
      
      tfPresent <- lapply(motifs@colData@rownames, function(x) {
        filtered[assay(filtered)[,x]] %>% length() > 0
      }) %>% unlist()
      
      names <- motifs@colData@rownames[tfPresent]
      names[order(names)]
  })

  observe({
    if(input$gene !=''){
      updatePickerInput(
        session,
        inputId = 'tf_motifs',
        selected = NULL,
        choices = tfMotifFilter()
      )
    }

  })

  values <- reactiveValues()
  observe({
    values$tf_motifs <- function(x){input$tf_motifs}
    values$tf_motifs_d <- debounce(values$tf_motifs, 2000)
    
    values$ymax <- function(x){input$ymax}
    values$ymax_d <- debounce(values$ymax, 2000)
    
    values$xrange <- function(x){input$xrange}
    values$xrange_d <- debounce(values$xrange, 2000)
    
    values$cluster_id <- function(x){input$clusterID}
    values$cluster_id_d <- debounce(values$cluster_id, 2000)
  })
  
  gvizPlot <- reactive({
    
    plotGenomeView(
      gene.symbol = input$gene,
      coverage.list = coverage.list,
      ylims = c(0, values$ymax_d()),
      coords = gvizCoords(),
      chr = gvizCoords()[[1]],
      beg = values$xrange_d()[1],
      END = as.numeric(values$xrange_d()[2]),
      transcriptID = input$transcriptID,
      corCut = input$cor_cut,
      pval_cut = input$pval_cut,
      cluster_id = values$cluster_id_d(),
      motifs_list = values$tf_motifs_d()
    )
  })
  

  output$gviz <- renderPlot({
    gvizPlot()
  })
  
  outputOptions(output, 'gviz', suspendWhenHidden = FALSE)
  outputOptions(output, 'rnaExpression.ui', suspendWhenHidden = FALSE)
  
  # # Plot clusters and motifs tracks
  # gvizPlotClust <- reactive({
  # 
  #   plot_clust_motif(
  #     gene.symbol = input$gene,
  #     coverage.list = coverage.list,
  #     ylims = c(0, input$ymax),
  #     coords = gvizCoords(),
  #     chr = gvizCoords()[[1]],
  #     beg = input$xrange[1],
  #     END = as.numeric(input$xrange[2]),
  #     transcriptID = input$transcriptID,
  #     corCut = input$cor_cut,
  #     pval_cut = input$pval_cut,
  #     cluster_id = input$clusterID,
  #     motifs_list = input$tf_motifs
  #   )
  # })
  # 
  # output$gvizClust <- renderPlot({
  #   gvizPlotClust()
  # })
  
  # Peaks table
  cor.gr_table <- reactive({
    if(is.na(input$cor_cut)){
      cor.gr.subset <- cor.gr[(elementMetadata(cor.gr)$transcript_id == input$transcriptID) &
                                (elementMetadata(cor.gr)$cluster.name %in% input$clusterID) &
                                (elementMetadata(cor.gr)$vs.null.p.value <= input$pval_cut) &
                                start(cor.gr) > input$xrange[1] & 
                                end(cor.gr) < as.numeric(input$xrange[2])]
    } else{
      cor.gr.subset <- cor.gr[(elementMetadata(cor.gr)$transcript_id == input$transcriptID) &
                                (elementMetadata(cor.gr)$estimate >= input$cor_cut) & 
                                (elementMetadata(cor.gr)$cluster.name %in% input$clusterID) &
                                (elementMetadata(cor.gr)$vs.null.p.value <= input$pval_cut) &
                                start(cor.gr) > input$xrange[1] & 
                                end(cor.gr) < as.numeric(input$xrange[2])]
    }
    dplyr::as_tibble(cor.gr.subset) %>% 
      dplyr::select(cluster.color, gene.symbol, transcript_id, 
                    cluster.name, estimate, vs.null.p.value, dplyr::everything(),
                    -width, -strand, -mean.gene.corr, -sd.gene.corr, 
                    -ncorrs, -KM3.ord)
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
    ) %>% formatStyle('cluster.color',  
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
  
  output$tf_legend <- renderPlot({
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
