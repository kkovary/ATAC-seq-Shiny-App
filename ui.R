

ui <- fluidPage(
  fluidRow(
    column(4, h4('ATAC-Seq Explorer'),
           downloadButton("report", "Download Summary")
           ),
    
    column(4,
           h5(''),
           selectInput(
             'accessionType',
             h4('Accession Type'),
             c('Gene Symbol' = 'gene.symbol', 'Transcript ID' = 'transcript_id')
           ),
           sliderInput(
             'xrange',
             label = h4('Track x-axis limits'),
             min = 0,
             max = 0,
             value = c(0,0)
           )
           
    ),
    column(
      4,
      h5(''),
      selectizeInput(
        'gene',
        label = h4('Search'),
        choices = '',
        options = list(maxOptions = 5, maxItems = 1),
        multiple = TRUE
      ),
      sliderInput(
        'ymax',
        label = h4('Track y-axis limit'),
        min = 0,
        max = 1000,
        value = 100
      )
    )
  ),
  
  hr(),
  textOutput('xlims'),
  
  # withSpinner(
  #   uiOutput("rnaExpression.ui"),
  #   #plotOutput('rnaExpression', height = 400*plotLen),
  #   type = getOption('spinner.type', default = 4)
  # ),
  h5('ATAC-seq Tracks'),
  plotOutput('gviz', height = 2500),
  hr(),
  h5('Expression'),
  busyIndicator(text = ''),
  uiOutput("rnaExpression.ui")
)
