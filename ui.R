

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
      label = h5('Gene Search'),
      choices = '',
      options = list(maxOptions = 5, maxItems = 1),
      multiple = TRUE
    ),
    sliderInput(
      'xrange',
      label = h5('Genome Coordinates'),
      min = 0,
      max = 0,
      value = c(0,0)
    ),
    sliderInput(
      'ymax',
      label = h5('Normalized ATAC-seq Coverage'),
      min = 0,
      max = 1000,
      value = 100
    ),
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
      value = 0,
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

    downloadButton("report", h6("Download Summary"))
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
