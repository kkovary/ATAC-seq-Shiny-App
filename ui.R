

ui <- fluidPage(
  
  sidebarPanel(
    titlePanel(strong('ATAC-Seq Explorer')),
    width = 2,
    # selectInput(
    #   'accessionType',
    #   h4('Accession Type'),
    #   c('Gene Symbol' = 'gene.symbol', 'Transcript ID' = 'transcript_id')
    # ),
    selectizeInput(
      'gene',
      label = h4('Gene Search'),
      choices = '',
      options = list(maxOptions = 5, maxItems = 1),
      multiple = TRUE
    ),
    sliderInput(
      'xrange',
      label = h4('Track x-axis limits'),
      min = 0,
      max = 0,
      value = c(0,0)
    ),
    sliderInput(
      'ymax',
      label = h4('Track y-axis limit'),
      min = 0,
      max = 1000,
      value = 100
    ),
    hr(style="border-color: grey"),
    h5('Peak-gene links'),
    selectInput(
      'transcriptID',
      label = h4('Transcript ID'),
      choices = '',
      multiple = FALSE
    ),
    numericInput(
      'cor_cut',
      h4('Correlation Cutoff'),
      value = 0,
      min = -1,
      max = 1,
      step = 0.05
    ),
    selectInput(
      'clusterID',
      label = h4('Cluster ID'),
      choices = c('All',as.vector(unique(cor.gr$cluster.name))),
      multiple = FALSE
    ),
    downloadButton("report", "Download Summary")
  ),
  

  # withSpinner(
  #   uiOutput("rnaExpression.ui"),
  #   #plotOutput('rnaExpression', height = 400*plotLen),
  #   type = getOption('spinner.type', default = 4)
  # ),
  mainPanel(
    #h5('ATAC-seq Tracks'),
    width = 10,
    plotOutput('gviz', height = 1100),
    hr(),
    h5('Expression'),
    busyIndicator(text = ''),
    uiOutput("rnaExpression.ui")
  )
  
  
)