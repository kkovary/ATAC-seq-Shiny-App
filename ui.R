appCSS <- "
#loading-content {
  position: absolute;
  background: #000000;
  opacity: 0.9;
  z-index: 100;
  left: 0;
  right: 0;
  height: 100%;
  text-align: center;
  color: #FFFFFF;
}
"

ui <- navbarPage(title = 'ATAC-Seq',
                 tabPanel('Introduction',
                          shiny::tags$head(includeHTML(('HTML/googleAnalytics.html'))),
                          shiny::includeHTML('HTML/ATACBrowserText.html')
                 ),
                 tabPanel('Explore Data',
                          shiny::tags$head(includeHTML(('HTML/googleAnalytics.html'))),
                          useShinyjs(),
                          
                          # Create loading message that hides app until
                          # it is loaded
                          inlineCSS(appCSS),
                          div(
                            id = "loading-content",
                            h2("Loading...")
                          ),
                          sidebarPanel(
                            titlePanel(h2('ATAC-Seq Explorer', align = 'center')),
                            width = 2,
                            chooseSliderSkin("Modern", color = '#428bca'),
                            shiny::tags$style(type='text/css', '#selected_tf_motifs {background-color: #DCDCDC;}'), 
                            #shiny::tags$head(shiny::tags$style(type='text/css', ".slider-animate-button { font-size: 1pt !important; }")),
                            # selectInput(
                            #   'accessionType',
                            #   h4('Accession Type'),
                            #   c('Gene Symbol' = 'gene.symbol', 'Transcript ID' = 'transcript_id')
                            # ),
                            
                            h4(''),
                            hr(style="border-color: grey"),
                            
                            selectInput('search_type', 
                                         label = 'Search Type', 
                                         choices = c('Gene Name','Coordinates','SNP ID'),
                                         selected = 'Gene Name'
                            ),
                            textInput(
                              "snp_search",
                              label = 'SNP ID Search',
                              value = ''
                            ),
                            selectInput('coord_chr', 
                                        label = 'Chromosome Number', 
                                        choices = chromosomes
                                        ),
                            numericInput('coord_min', label = 'Min', value = 0),
                            numericInput('coord_max', label = 'Max', value = 0),
                            pickerInput(
                              inputId = "gene",
                              label = 'Gene Search',
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
                              label = 'Genome Coordinates',
                              min = 0,
                              max = 0,
                              value = c(0,0),
                              round = T
                            ),
                            sliderInput(
                              'ymax',
                              label = 'Y-axis scale - ATAC-seq signal',
                              min = 0,
                              max = 1000,
                              value = 200
                            ),
                            shiny::actionButton("plot_button", h6("Refresh\nPlot"), icon = icon("refresh"),
                                                style="color: #fff; background-color: #64dd17; white-space:normal; width:45%; height:75px"),
                            downloadButton("report", h6("Download\nSummary"),
                                           style = "white-space:normal; width:45%; height:75px"),
                            hr(style="border-color: grey"),
                            h4('Peak-gene links'),
                            h6('We defined putative enhancer-gene linkages on the basis of the correlation of accessibility and gene expression across cell lineages. These parameters filter the data displayed by features of those correlations.'),
                            selectInput(
                              'transcriptID',
                              label = h5('Transcript ID'),
                              choices = '',
                              multiple = FALSE
                            ),
                            h5('Correlation Cutoff'),
                            radioButtons('greater_less',
                                         label = NULL,
                                         choices = c('Greater than','Less than'),
                                         inline = TRUE),
                            numericInput(
                              'cor_cut',
                              label = NULL,
                              value = NA,
                              min = -1,
                              max = 1,
                              step = 0.05
                            ),
                            numericInput(
                              'pval_cut',
                              'p.value Cutoff',
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
                            textOutput("selected_tf_motifs"),
                            h6('Peaks containing the selected motifs will be displayed on the ATAC-seq plot.'),
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
                                       shiny::tags$style(type="text/css",
                                                         ".shiny-output-error { visibility: hidden; }",
                                                         ".shiny-output-error:before { visibility: hidden; }"
                                       ),
                                       dataTableOutput('peaks_table') %>% withSpinner(type = 6)
                              ),
                              tabPanel('ATAC-seq',
                                       shiny::tags$style(type="text/css",
                                                         ".shiny-output-error { visibility: hidden; }",
                                                         ".shiny-output-error:before { visibility: hidden; }"
                                       ),
                                       plotOutput('gviz', height = 1200) %>% withSpinner(type = 6),
                                       #plotOutput('gvizClust') %>% withSpinner(type = 6),
                                       plotOutput('tf_legend', height = 100)
                                       #busyIndicator(text = '')
                              )
                            )
                            #h5('ATAC-seq Tracks'),
                            
                            
                          )
                 )
                 
                 
                 
)

