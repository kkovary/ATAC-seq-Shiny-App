ui <- navbarPage(title = 'ATAC-Seq',
                 tabPanel('Introduction',
                          shiny::includeHTML('HTML/ATACBrowserText.html')
                 ),
                 tabPanel('Explore Data',
                          useShinyjs(),
                          sidebarPanel(
                            titlePanel(h2('ATAC-Seq Explorer', align = 'center')),
                            width = 2,
                            chooseSliderSkin("Modern", color = '#428bca'),
                            shiny::tags$style(type='text/css', '#selected_tf_motifs {background-color: #DCDCDC;}',
                                              ".shiny-output-error { visibility: hidden; }",
                                              ".shiny-output-error:before { visibility: hidden; }"), 
                            #shiny::tags$head(shiny::tags$style(type='text/css', ".slider-animate-button { font-size: 1pt !important; }")),
                            # selectInput(
                            #   'accessionType',
                            #   h4('Accession Type'),
                            #   c('Gene Symbol' = 'gene.symbol', 'Transcript ID' = 'transcript_id')
                            # ),
                            
                            h4(''),
                            hr(style="border-color: grey"),
                            
                          
                            
                            radioButtons('search_type', 
                                         label = 'Search Options', 
                                         choices = c('Gene Name','SNP ID'),
                                         selected = 'Gene Name',
                                         inline = TRUE
                                         ),
                            textInput(
                              "snp_search",
                              label = 'SNP ID Search',
                              value = ''
                            ),
                            verbatimTextOutput('snp_text'),
                            
                            pickerInput(
                              inputId = "gene",
                              label = h5('Gene Search'),
                              selected = 'NEUROG2',
                              choices = gene_names,
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
                              label = h5('Y-axis scale - ATAC-seq signal'),
                              min = 0,
                              max = 1000,
                              value = 200
                            ),
                            conditionalPanel(condition = "!$('html').hasClass('shiny-busy')",
                                             shiny::actionButton("plot_button", h6("Refresh\nPlot"), icon = icon("refresh"),
                                                                 style="color: #fff; background-color: #64dd17; border-color: #64dd17; white-space:normal; width:45%"),
                                             downloadButton("report", h6("Download\nSummary"),
                                                            style = "white-space:normal; width:45%")
                            ),
                            hr(style="border-color: grey"),
                            h4('Peak-gene links'),
                            h6('We defined putative enhancer-gene linkages on the basis of the correlation of accessibility and gene expression across cell lineages. These parameters filter the data displayed by features of those correlations.'),
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
                            uiOutput("cluster_select"),
                            # checkboxGroupInput('clusterID', label = h5('Cluster ID'),
                            #                    choiceNames = c('MG1','MG2','GP','PL','EA1','EA2','PN1','PN2','LN1','LN2','CS'),
                            #                    choiceValues = c('MG1','MG2','GP','PL','EA1','EA2','PN1','PN2','LN1','LN2','CS'),
                            #                    selected = c('MG1','MG2','GP','PL','EA1','EA2','PN1','PN2','LN1','LN2','CS'),
                            #                    inline = F),
                            h4(''),
                            hr(style="border-color: grey"),
                            pickerInput(
                              inputId = "tf_motifs",
                              label = "Select TF motifs to display",
                              choices = motifs@colData@rownames,
                              options = list(
                                `actions-box` = TRUE,
                                size = 10,
                                `selected-text-format` = "count > 3",
                                `live-search` = TRUE,
                                `max-options` = 10,
                                `max-options-text` = "Max of 10 selected"
                              ),
                              multiple = TRUE
                            ),
                            textOutput("selected_tf_motifs")
                          ),
                          # withSpinner(
                          #   uiOutput("rnaExpression.ui"),
                          #   #plotOutput('rnaExpression', height = 400*plotLen),
                          #   type = getOption('spinner.type', default = 4)
                          # ),
                          mainPanel(
                            width = 10,
                            tabsetPanel(
                              tabPanel('RNA-seq',
                                       uiOutput("rnaExpression.ui")
                                       #plotOutput('rnaExpression') %>% withSpinner(type = 6)
                              ),
                              tabPanel('Table',
                                       dataTableOutput('peaks_table') %>% withSpinner(type = 6)
                              ),
                              tabPanel('ATAC-seq',
                                       plotOutput('gviz', height = 1100) %>% withSpinner(type = 6),
                                       #plotOutput('gvizClust') %>% withSpinner(type = 6),
                                       plotOutput('tf_legend', height = 100)
                                       #busyIndicator(text = '')
                              )
                            )
                            #h5('ATAC-seq Tracks'),
                            
                            
                          )
                 )
                 
                 
                 
)

