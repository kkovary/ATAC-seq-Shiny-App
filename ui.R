

ui <- fluidPage(
  
  sidebarPanel(
    titlePanel(strong('ATAC-Seq Explorer')),
    selectInput(
      'accessionType',
      h4('Accession Type'),
      c('Gene Symbol' = 'gene.symbol', 'Transcript ID' = 'transcript_id')
    ),
    selectizeInput(
      'gene',
      label = h4('Search'),
      choices = '',
      options = list(maxOptions = 5, maxItems = 1),
      multiple = TRUE
    ),
    selectInput(
      'transcriptID',
      label = h4('Transcript ID'),
      choices = '',
      multiple = FALSE
    ),
    numericInput(
      'cor_cut',
      'Correlation Cutoff',
      value = 0.2,
      min = -1,
      max = 1,
      step = 0.05
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
    downloadButton("report", "Download Summary")
  ),
  

  # withSpinner(
  #   uiOutput("rnaExpression.ui"),
  #   #plotOutput('rnaExpression', height = 400*plotLen),
  #   type = getOption('spinner.type', default = 4)
  # ),
  mainPanel(
    #h5('ATAC-seq Tracks'),
    plotOutput('gviz', height = 1100),
    hr(),
    h5('Expression'),
    busyIndicator(text = ''),
    uiOutput("rnaExpression.ui")
  )
  
  
)
