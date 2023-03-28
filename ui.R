library(shiny)
library(markdown)

shinyUI(navbarPage("Baseline Correction for Copy Number Data from Cancer Samples", id="baseCN",
  tabPanel("Description",
          fluidPage(titlePanel("Description"),
                    fluidRow(column(5,includeMarkdown("description.md")),
                             column(4,uiOutput('LoadDataButtons')))  
          )),

  tabPanel("Upload region",
           fluidPage(titlePanel("Upload Regions data"),
                     sidebarLayout(
                     sidebarPanel(
                         fileInput('file1', 'Choose Regions CSV File', accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
                         tags$hr(),
                         tags$div(class="header", checked=NA,
                                  tags$p("CSV manipulation")),

                         checkboxInput('header', 'Header', TRUE),radioButtons('sep', 'Separator',c(Comma=',',Semicolon=';',Tab='\t'),','),
                         radioButtons('quote', 'Quote',c(None='','Double Quote'='"','Single Quote'="'"),'"'),
                         tags$div(tags$p("Select which column is:")),
                         selectInput("RegionSample", "Sample:",NULL),
                         selectInput("RegionChromosome", "Chromosome:",NULL),
                         selectInput("Regionbpstart", "bp.Start:",NULL),
                         selectInput("Regionbpend", "bp.End:",NULL),
                         selectInput("RegionNumMark", "Num.of.Markers:",NULL),
                         selectInput("RegionMean", "Mean:",NULL)
                     ),

                     mainPanel(tableOutput('csvtableRegions'),
                               uiOutput('regionsbuttonsGo2Sample'),
                               uiOutput('regionsbuttonsGo2PlotRaw'))
                     )
                )),
  tabPanel("Upload sample list",
           fluidPage(titlePanel("Optional Upload sample list"),
                     sidebarLayout(
                         sidebarPanel(
                             fileInput('file2', '*Optional* Choose samplenamne  CSV File',accept=c('text/csv', 'text/comma-separated-values,text/plain','.csv')),
                             tags$hr(),
                             tags$div(class="headersamp", checked=NA,
                                      tags$p("CSV manipulation")),

                             checkboxInput('headersamp', 'Header', TRUE),radioButtons('sepsamp', 'Separator',c(Comma=',',Semicolon=';',Tab='\t'),','),
                             radioButtons('quotesamp', 'Quote',c(None='','Double Quote'='"','Single Quote'="'"),'"'),
                             tags$div(tags$p("Select which column is:")),
                             selectInput("SampleNumber", "Number:",NULL),
                             selectInput("SampleSample", "Sample:",NULL),
                             selectInput("SampleComment", "Comment:",NULL)
                         ),

                         mainPanel(tableOutput('csvtableSample'),
                                   uiOutput('sampleButtonG2Raw'))
                     )
             )),

  tabPanel("TCGA",
           fluidPage(titlePanel("TCGA"),
	                 sidebarLayout(
                         sidebarPanel(uiOutput('tcga'), uiOutput('tcgaSamplenumber')),
		                 mainPanel(tableOutput('tableTCGA'),
		                           uiOutput('tcgaButtonG2Raw'))
                     )
           )),

  tabPanel("Plot raw",
           fluidPage(titlePanel("Raw plot"),
                     sidebarLayout(
                         sidebarPanel(
                             sliderInput("NumberSampleSlider",
                                         "Number of Samples:",
                                         min=0,
                                         max=0,
                                         value=1,
                                         ticks=FALSE ),
                             sliderInput("NumberMarkerSlider",
                                         "Number of Markers:",
                                         min = 0,
                                         max = 100,
                                         value = 20,
                                         step=20),
                             sliderInput("NumberCutoffSlider",
                                         "Cutoff:",
                                         min = 0,
                                         max = 1,
                                         value = 0.1),
                             checkboxInput('ShowComments', 'Show Comments?', TRUE),
                             uiOutput('plotbuttonsGo2Correct')),
                             mainPanel(uiOutput("plotraw"))
                     )
           )),

  tabPanel("AutoCorrection",
           fluidPage(titlePanel("AutoCorrection"),         
                     wellPanel(fluidRow(uiOutput("downloadButtons"))),
                     fluidRow(column(12,uiOutput("autocorrection")))          
           ))
))
