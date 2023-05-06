library(shiny)
library(shinyWidgets)
library(shinyjs)
library(shinybrowser)
library(colourpicker)
library (data.table)
library (magrittr)
if (!require(R.utils)){ warning("\nWARNING: R.utils is not installed. R.utils is required by data.table to read zipped files (ex. \".gz\" extensions).\n")}
#library (ComplexHeatmap)
#library (ggplot2)

source("./helpers.R")
source("./MSstats_V4_Functions.R")
source("./MSstats_Helper_Functions.R")

options(shiny.maxRequestSize = 1000 * 1024^2)


ui <- fluidPage(
  shinybrowser::detect(),
  useShinyjs(),
  tags$head(
    tags$style(HTML("hr {border-top: 1px solid #000000;}"))
  ),
  
  titlePanel("Proteomics Data Preparation with MSstats"),
  sidebarLayout(
    sidebarPanel( id = "sidebar",
                  strong("Upload DIA or DDA peptide data to summarize as sites or proteins and calculate l2fc."),# Preliminary QC included to entirely exclude runs with poor data quality."),
                  h3(""),
                  fileInput("file", "Upload MaxQuant evidence.txt or Spectronaut report file here."),
                  fileInput("keys", "Upload keys file here."), 
                  
                  hr(),
                  h3("Data Filtering Settings"),
                  div(style="display:inline-block", numericInput("magn", 
                               label = "Absolute Magnitude Threshold (log scale):",
                               value = 5, step = 1)),
                  div(style="display:inline-block", actionButton("magnFilt", "Filter")),
                  div(style="display:inline-block", helpText(" ") ), div(style="display:inline-block", helpText(" ") ), div(style="display:inline-block", helpText(" ") ), div(style="display:inline-block", helpText(" ") ),
                  div(style="display:inline-block",numericInput("ptms", 
                               label = "PTM Confidence Threshold:",
                               min = 0.001, max = 1, value = 0.75, step = 0.05, width = "282px")),
                  div(style="display:inline-block", actionButton("ptmsFilt", "Filter")),
                  div(style= "display:block", helpText(" ")),
                  textAreaInput("exclude", "Runs to Exclude (by label or rawfile)", resize = "vertical", row = 3 ),
                  
                
                  hr(),
                  h3("MSstats Settings"),
                  selectInput("species", "Select Data Type", c("Protein", "PTM"), selected = "Protein"),
                  checkboxInput('bioRep', "Treat bioreplicates as related"),
                  textAreaInput("regContrast", "RegEx Contrast Labels", resize = "vertical", row = 3 ),
                  fileInput("contrasts", "Upload contrast file here (Optional).")
                  
                  
    ),
    
    mainPanel(plotOutput("intensityHist"),
     id="main" ) 
    )
)

server <- function(input, output, session) {
  
  ### UI Functions ###
  observe({
    if (!is.null(input$contrasts)){
      updateTextAreaInput(session, "regContrast", value = sourcetools::read(input$contrasts$datapath) )
    } 
  })
  
  ### Data Filtering ###
  inputData <- reactive({
    if(!is.null(input$file)){
        applyLabels(fread(input$file$datapath), fread(input$keys$datapath))
    } else { NULL }
  })
  
  inputDatalf <- reactive({
    if(!is.null(inputData())){
      
    } else { NULL }
  })
  
  #observeEvent({ input$magnFilt,
    
  # })
  
  ### Running MSStats ###
  
  
  ### Rendering in Main Window ###
  output$intensityHist <- renderPlot({
    intensityHistPlot(inputData())
  })
  
}

shinyApp(ui = ui, server = server)
