shinyUI(
  navbarPage('ATACseqExplorer',
             tabPanel('Graphs',
                      sidebarPanel(
                        titlePanel(strong("Title")),
                        
                        h6(em("Authors")),
                        
                        selectInput('accessionType', 'Accession Type',
                                    c('Gene Symbol' = 'gene.symbol','Transcript ID' = 'transcript_id')),
                        
                        selectizeInput('gene', label = 'Search', choices = '',
                                       options = list(maxOptions = 5,maxItems = 1),
                                       multiple = TRUE)
                        
                      ),
                      
                      mainPanel(
                        plotOutput("rnaExpression")
                      )
             )
  )
)
