ui <- fluidPage(
  
  sidebarPanel(
    titlePanel(strong('ATAC-Seq Explorer')),
    width = 2,
    #setSliderColor(c('#428bca','#428bca'),c(1,3)),
    chooseSliderSkin("Modern"),
    #shiny::tags$head(shiny::tags$style(type='text/css', ".slider-animate-button { font-size: 1pt !important; }")),
    # selectInput(
    #   'accessionType',
    #   h4('Accession Type'),
    #   c('Gene Symbol' = 'gene.symbol', 'Transcript ID' = 'transcript_id')
    # ),
    pickerInput(
      inputId = "gene",
      label = h5('Gene Search'),
      choices = as.vector(unique(data$gene.symbol)),
      options = list(
        `live-search` = TRUE,
        title = 'Gene Name'
        )
    ),
    # selectizeInput(
    #   'gene',
    #   label = h5('Gene Search'),
    #   choices = '',
    #   options = list(maxOptions = 5, maxItems = 1),
    #   multiple = TRUE
    # ),
    sliderInput(
      'xrange',
      label = h5('Genome Coordinates'),
      min = 0,
      max = 0,
      value = c(0,0),
      round = T
    ),
    sliderInput(
      'ymax',
      label = h5('Normalized ATAC-seq Coverage'),
      min = 0,
      max = 1000,
      value = 100
    ),
    shiny::actionButton("plot_button", "Plot", icon = icon("refresh"), 
                        style="color: #fff; background-color: #28a745; border-color: #28a745"),
    shiny::tags$p(h6(em('please wait for genome coordinates to update'))),
    hr(style="border-color: grey"),
    h4('Peak-gene links'),
    selectInput(
      'transcriptID',
      label = h5('Transcript ID'),
      choices = '',
      multiple = FALSE
    ),
    numericInput(
      'cor_cut',
      h5('Correlation Cutoff'),
      value = NA,
      min = -1,
      max = 1,
      step = 0.05
    ),
    numericInput(
      'pval_cut',
      h5('p.value Cutoff'),
      value = 1,
      min = 0,
      max = 1,
      step = 0.01
    ),
    checkboxGroupInput('clusterID', label = h5('Cluster ID'),
                       choiceNames = c('MG1','MG2','GP','PL','EA1','EA2','PN1','PN2','LN1','LN2','CS'),
                       choiceValues = c('MG1','MG2','GP','PL','EA1','EA2','PN1','PN2','LN1','LN2','CS'),
                       selected = c('MG1','MG2','GP','PL','EA1','EA2','PN1','PN2','LN1','LN2','CS'),
                       inline = F),
    
    pickerInput(
      inputId = "myPicker",
      label = "Select/deselect TF Motifs",
      choices = motifs@colData@rownames,
      options = list(
        `actions-box` = TRUE,
        size = 10,
        `selected-text-format` = "count > 3",
        `live-search` = TRUE
      ),
      multiple = TRUE
    ),

    downloadButton("report", h6("Download Summary"))
  ),
  

  # withSpinner(
  #   uiOutput("rnaExpression.ui"),
  #   #plotOutput('rnaExpression', height = 400*plotLen),
  #   type = getOption('spinner.type', default = 4)
  # ),
  mainPanel(
    width = 10,
    tabsetPanel(
      tabPanel('Plots',
               plotOutput('gviz', height = 1100),
               hr(),
               busyIndicator(text = ''),
               uiOutput("rnaExpression.ui")
               ),
      tabPanel('Table',
               dataTableOutput('peaks_table')
               )
    ),
    #h5('ATAC-seq Tracks'),
    
    
  )
  
  
)
