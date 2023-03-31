library(shiny)
library(shinyWidgets)
library(shinyjs)
library(colourpicker)
library (data.table)
if (require(R.utils)){ warning("\nWARNING: R.utils is not installed. R.utils is required by data.table to read zipped files (ex. \".gz\" extensions).\n")}
library (ComplexHeatmap)
library (ggplot2)
source("./helpers.R")
options(shiny.maxRequestSize = 1000 * 1024^2)


ui <- fluidPage(
  
  useShinyjs(),
  tags$head(
    tags$style(HTML("hr {border-top: 1px solid #000000;}"))
  ),
  
  titlePanel("Proteomics Summary Analyses"),
  sidebarLayout(
    sidebarPanel( id = "sidebar",
                  strong("Render and download a Volcano Plot, PCA Plot, Effect Summary Stacked Bar Chart, and Complex Heatmap from your MSStats output files."),
                  
                  hr(),
                  h3("General Settings"),
                  fileInput("results", "Upload MSStats group comparison results file here."),
                  fileInput("intensities", "Upload MSStats protein level data file here."),
                  
                  selectInput("species", "Select Species", c("HUMAN", "MOUSE", "RAT", "OTHER"), selected = "HUMAN"),
                  
                  numericInput("fdr", 
                               label = "FDR adjusted p value threshold:",
                               min = 0, max = 1, value = 0.05, step = 0.01),
                  
                  numericInput("foldchange", 
                               label = "Absolute log2 fold change threshold:",
                               min = 0.001, max = 5, value = 1.0, step = 0.5),
                
                  div(style="display:inline-block", colourInput("pcol", "Select positive color", "#DB4F4F")),
                  div(style="display:inline-block", helpText(" ") ), div(style="display:inline-block", helpText(" ") ),
                  div(style="display:inline-block",colourInput("ncol", "Select negative color", "#4E53D9")),
                  div(style= "display:block", helpText(" ")),
                  
                  textInput("dlprepend", "Label for plot downloads."),
                  
                  hr(),
                  h3("Plot Settings"),
                  tabsetPanel(type = "tabs", id = "tabside",
                              tabPanel("Volcano Plot",
                                       checkboxInput("separate", "Separate volcano plots by contrast?", value = TRUE),
                                       checkboxInput("volcLabel", "Label Significant Outliers?"),
                                       downloadButton("dlvolcano", "Download")),
                              tabPanel("Effect Summary Bar Chart", 
                                       #textInput("barchartCondLabel", "Positive Condition Label (ex. 'Infected')"),
                                       div(style="display:inline-block", colourInput("ipcol", "Infinite positive color", "#9E0000")),
                                       div(style="display:inline-block", helpText(" ") ), div(style="display:inline-block", helpText(" ") ),
                                       div(style="display:inline-block",colourInput("incol", "Infinite negative color", "#000294")),
                                       div(style="display:block", helpText(" ") ),
                                       downloadButton("dlbarchart", "Download")),
                              tabPanel("PCA Plot", 
                                       checkboxInput("pcaCircleConditions", "Circle replicates in same condition?"),
                                       #checkboxInput("pcaLines", "Draw lines between time points for each condition?"),
                                       #helpText("(requires time series formatted as 'conditionxyz-01h')"),
                                       checkboxInput("pcaNoLabel", "Remove Labels?"),
                                       checkboxInput("pcaNoLegend", "Remove Legend?"),
                                       downloadButton("dlpca", "Download")),
                              tabPanel("Heatmap Plot",
                                       selectInput("selectMat", "Select normalization for heatmap:", 
                                                   c("Median Swept", "Normalized to Control", "No normalization"), selected = "Median Swept"),
                                       checkboxInput("hmcluster", "Order heatmap columns by hierarchical clustering?"),
<<<<<<< HEAD
                                       checkboxInput("hmsignificant", "Exclude insignificant proteins/sites from heatmap?"),
=======
                                       checkboxInput("hmsignificant", "Exclude insignificant proteins/sites heatmap?"),
>>>>>>> 379365f47b8fa8d462506f3a1ddee78b27261821
                                       checkboxInput("hmsample", "Render entire heatmap instead of randomly sampling 1000 rows? \n(this may take a few minutes depending on the size of your data)"),
                                       downloadButton("dlheatmap", "Download")
                              ))
                  
    ),
    
    mainPanel(
      materialSwitch("tglSide", label = "Toggle Sidebar", value = TRUE, status = "danger"), 
      #textOutput("selected_fdr"),
      #textOutput("selected_l2fc"),
      #textOutput("intMatText"),
      hr(),
      tabsetPanel(type = "tabs", id = "tabset",
                  tabPanel("Volcano Plot", plotOutput("volcanoPlot")),
                  tabPanel("Effect Summary Bar Chart", plotOutput("stackedBarPlot")),
                  tabPanel("PCA Plot", plotOutput("pcaPlot")),
                  tabPanel("Heatmap Plot", plotOutput("heatmapPlot"))
      ), id="main"
    )
    
  )
)

server <- function(input, output) {
  
  # OVERHEAD
  observeEvent(input[["tglSide"]], {
    if(input[["tglSide"]] == FALSE){
      hideElement(selector = "#sidebar")
      removeCssClass("main", "col-sm-8")
      addCssClass("main", "col-sm-12")
    }else{
      showElement(selector = "#sidebar")
      removeCssClass("main", "col-sm-12")
      addCssClass("main", "col-sm-8")
    }
  })
  
  # DATA PROCESSING
  inFile <- reactive({
    input$results
  })
  
  inFile2 <- reactive({
    input$intensities
  })
  
  dataType <- reactive({
    if(is.null(input$results)){
      "Abundance"
    } else{
      checkDataType(fread(inFile()$datapath))
    }
  })
  
  resultsForm <- reactive({       # Load and format results with names, only runs when new results file is uploaded
    if (is.null(inFile())){
      data.table(log2FC = c(0), adj.pvalue = c(0.99), posGROUP = c("Waiting For Input Data"), effectSign = c("notSig"))
    } else {
      formatResults( fread(inFile()$datapath),  dataType(), input$species)
    }
  })
  
  resultsData <- reactive({       # Apply simple thresholds to results data for volc plot
    if (is.null(inFile())){
      data.table(log2FC = c(0), adj.pvalue = c(0.99), posGROUP = c("Waiting For Input Data"), effectSign = c("notSig"), dataGene = c(" "))
    } else {
      threshResults( resultsForm(),  input$fdr, input$foldchange)
    }
  })
  
  intensitiesData <- reactive({   # Load and format intensities data and output data structure with matrices
    if (is.null(inFile2())){
      list("data" = data.table(LogIntensities = c(0), adj.pvalue = c(0.99), GROUP = c("Waiting For Input Data")),
           "int.mat" =  matrix(rep(c(1,-1), 5000), nrow = 1000, ncol = 10, dimnames = list(rep(".Waiting For Input Data",1000),rep(".Waiting For Input Data",10))),
           "swept.mat" =  matrix(rep(c(1,-1), 5000), nrow = 1000, ncol = 10, dimnames = list(rep("Waiting For Input Data",1000),rep("Waiting For Input Data",10))))
    } else {
      formatIntensities( fread(inFile2()$datapath), dataType(), input$species)   #"$data" = intensities, "$int.mat" = intensity.mat, "$norm.mat" = normIntensity.mat, 
      #"$swept.mat" = sweptIntensity.mat, "$ctlGroups" = c("ctl group names")
    }
  })
  
  resultsDataPlus <- reactive({
    if (is.null(inFile()) | is.null(inFile2())){
      data.table(log2FC = c(0), adj.pvalue = c(0.99), Label = c("Waiting For Results Input Data"), effectType = c("notSig"), dataGene = c(" "))
    } else {
      threshResultsPlus( resultsData(), intensitiesData()$data, input$fdr, input$foldchange)
    }
  })
  
  # OUTPUT RENDERING
  output$selected_fdr <- renderText({ 
    paste("FSR Adjusted p Values Less than", input$fdr) 
  })
  output$selected_l2fc <- renderText({ 
    paste("Absolute Log2 Fold Changes Greater than", as.character(input$foldchange) )
  })
  output$intMatText <- renderText({ 
    if (is.null(inFile2())){
      "Waiting For Intensities Input Data"
    } else {
      switch(input$selectMat,
             "Median Swept" = "Heatmap Intensity Data normalized by Row Median", 
             "Normalized to Control" = paste("Control Conditions Used to Normalize Intensity Data for Heatmap:", paste(intensitiesData()$ctlGroups, collapse = ", ")), 
             "No normalization" = "No Normalization Used for Heatmap Intensity Data.")
    }
  })
  
  # Volcano Tab
  output$volcanoPlot <- renderPlot({
    plotProteomicsVolcano(resultsData(), input$pcol, input$ncol, input$separate, input$volcLabel)
  }, height = 700)
  
  output$dlvolcano <- downloadHandler(
    filename = function(){
      paste(input$dlprepend, "_volcanoPlot_", Sys.Date() , ".pdf", sep="")
    },
    content = function(file) {
      makePdf(plotProteomicsVolcano(resultsData(), input$pcol, input$ncol, input$separate, input$volcLabel), file)
    }
  )
  
  # Bar Chart Tab
  output$stackedBarPlot <- renderPlot({
    plotStackedBarChart(resultsDataPlus(), dataType(), input$pcol, input$ncol, input$ipcol, input$incol)
  }, height = 700)
  
  output$dlbarchart <- downloadHandler(
    filename = paste(input$dlprepend, "_effectBarChart_", Sys.Date(), ".pdf", sep=""),
    content = function(file) {
      makePdf(plotStackedBarChart(resultsDataPlus(), dataType(), input$pcol, input$ncol, input$ipcol, input$incol), file)
    }
  )
  
  # PCA Tab
  output$pcaPlot <- renderPlot({
    plotPCA(intensitiesData()$int.mat, dataType(), input$pcaCircleConditions, input$pcaNoLabel, input$pcaNoLegend)
  }, height = 700)
  
  output$dlpca <- downloadHandler(
    filename = paste(input$dlprepend, "_PCA_", Sys.Date(), ".pdf", sep=""),
    content = function(file) {
      makePdf(plotPCA(intensitiesData()$int.mat, dataType(), input$pcaCircleConditions, input$pcaNoLabel, input$pcaNoLegend), file)
    }
  )
  
  # Heatmap Tab
  premat <- reactive({
    switch(input$selectMat,
           "Median Swept" = intensitiesData()$swept.mat, 
           "Normalized to Control" = intensitiesData()$normIntensity.mat, 
           "No normalization" = intensitiesData()$intensity.mat) })
  
  mat <- reactive({
    if(input$hmsignificant){
      filterSigIntensities(premat(), resultsDataPlus())
    } else {premat()} })
  
  sample <- reactive({if (input$hmsample == T){
    NULL
  } else {1000}})
  
  output$heatmapPlot <- renderPlot({
    plotHeatmap(mat(), sample(), input$hmcluster )
  }, height = 700)
  
  output$dlheatmap <- downloadHandler(
    filename = paste(input$dlprepend, "_Heatmap_", Sys.Date(), ".pdf", sep=""),
    content = function(file) {
      makePdf(plotHeatmap(mat(), sample(), input$hmcluster), file)
    }
  )
  
  
  
  
}

shinyApp(ui = ui, server = server)
