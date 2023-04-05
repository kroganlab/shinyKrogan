library(shiny)
library(shinyWidgets)
library(shinyjs)
library(shinybrowser)
library(colourpicker)
library (data.table)
if (!require(R.utils)){ warning("\nWARNING: R.utils is not installed. R.utils is required by data.table to read zipped files (ex. \".gz\" extensions).\n")}
library (ComplexHeatmap)
library (ggplot2)
source("./helpers.R")
options(shiny.maxRequestSize = 1000 * 1024^2)


ui <- fluidPage(
  shinybrowser::detect(),
  useShinyjs(),
  tags$head(
    tags$style(HTML("hr {border-top: 1px solid #000000;}"))
  ),
  
  titlePanel("Proteomics Summary Analyses"),
  sidebarLayout(
    sidebarPanel( id = "sidebar",
                  strong("Render and download a Volcano Plot, PCA Plot, Effect Summary Stacked Bar Chart, and Complex Heatmap from your MSStats output files. Both files must be uploaded."),
                  h3(""),
                  fileInput("results", "Upload MSStats group comparison results file here."),
                  fileInput("intensities", "Upload MSStats protein level data file here."),
                  actionButton("loadSampleData", "Load Sample Data"), 
                  hr(),
                  h3("General Settings"),
                  
                  textInput("experimentName", "Experiment/data name for figure title and download."),
                  
                  selectInput("species", "Select Species", c("HUMAN", "MOUSE", "RAT", "OTHER"), selected = "HUMAN"),
                  
                  div(style="display:inline-block", numericInput("fdr", 
                               label = "FDR adjusted p value threshold:",
                               min = 0, max = 1, value = 0.05, step = 0.01)),
                  div(style="display:inline-block", helpText(" ") ), div(style="display:inline-block", helpText(" ") ), div(style="display:inline-block", helpText(" ") ), div(style="display:inline-block", helpText(" ") ),
                  div(style="display:inline-block",numericInput("foldchange", 
                               label = "Absolute log2 fold change threshold:",
                               min = 0.001, max = 5, value = 1.0, step = 0.5)),
                  div(style= "display:block", helpText(" ")),
                
                  div(style="display:inline-block", colourInput("pcol", "Select positive color", "#DB4F4F")),
                  div(style="display:inline-block", helpText(" ") ), div(style="display:inline-block", helpText(" ") ), div(style="display:inline-block", helpText(" ") ), div(style="display:inline-block", helpText(" ") ),
                  div(style="display:inline-block",colourInput("ncol", "Select negative color", "#4E53D9")),
                  div(style= "display:block", helpText(" ")),
                  
                  hr(),
                  h3("Plot Settings"),
                  tabsetPanel(type = "tabs", id = "tabside",
                              tabPanel("Volcano Plot",
                                       checkboxInput("separate", "Separate volcano plots by contrast?", value = TRUE),
                                       checkboxInput("volcLabel", "Label Significant Outliers?"),
                                       checkboxInput("volcLines", "Show log2FC and pValue cutoffs?"),
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
                                       checkboxInput("hmsignificant", "Exclude insignificant proteins/sites from heatmap?"),
                                       checkboxInput("hmsample", "Render entire heatmap instead of randomly sampling 1000 rows? \n(this may take a few minutes depending on the size of your data)"),
                                       downloadButton("dlheatmap", "Download")
                              ))
                  
    ),
    
    mainPanel(
      div(style="display:inline-block", materialSwitch("tglSide", label = "Toggle Sidebar", value = TRUE, status = "danger")), 
      div(style="display:inline-block", helpText(" ") ), div(style="display:inline-block", helpText(" ") ), div(style="display:inline-block", helpText(" ") ), div(style="display:inline-block", helpText(" ") ),
      div(style="display:inline-block; color:red", textOutput("errorMsg")),
      div(style="display:block", helpText("")),
      #textOutput("selected_fdr"),
      #textOutput("selected_l2fc"),
      #textOutput("intMatText"),
      hr(),
      tabsetPanel(type = "tabs", id = "tabset",
                  tabPanel("Volcano Plot",       div(style="display:block", helpText("")),
                                                 plotOutput("volcanoPlot")),
                  tabPanel("Effect Summary Bar Chart", div(style="display:block", helpText("")),
                                                       plotOutput("stackedBarPlot")),
                  tabPanel("PCA Plot", div(style="display:block", helpText("")),
                                                       plotOutput("pcaPlot")),
                  tabPanel("Heatmap Plot", div(style="display:block", helpText("")),
                                                       plotOutput("heatmapPlot"))
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
  
  
  dimensionh <- reactive({
    if (input$tglSide == T| is.null(input$tglSide) == T){
      0.6 * shinybrowser::get_height()
    } else{
      0.9 * shinybrowser::get_height()
    }
      })
  
  dimensionw <-  reactive({
    if (input$tglSide == T| is.null(input$tglSide) == T){
      0.6 * shinybrowser::get_width()
    } else{
      0.9 * shinybrowser::get_width()
    }
  })
  
  
  
  # DATA PROCESSING
  inFile <- reactive({
    if (is.null(input$results)){
      if (input$loadSampleData == 0){
        NULL
      } else{
        list(datapath = "./sample_data/2022_12_27_Sample_Phospho_GroupComparisonResult.csv", name = "2022_12_17_Sample Phospho ")
      }
    } else {
      input$results
    }
  })
  
  inFile2 <- reactive({
    if (is.null(input$intensities)){
      if (input$loadSampleData == 0){
        NULL
      } else{
        list(datapath = "./sample_data/2022_12_27_Sample_Phospho_ProteinLevelData.csv.gz", name = "2022_12_17_Sample Phospho ")
      }
    } else {
      input$intensities
    }
  })
  
  #inFile2 <- reactive({
  #  input$intensities
  #})
  
  dataType <- reactive({
    if(is.null(inFile())){
      if(is.null(inFile2())){
        "Abundance"
        
      } else{
      checkDataType(fread(inFile2()$datapath))
      }
    } else{
      checkDataType(fread(inFile()$datapath))
    }
  })
  
  dataName <- reactive({
    if(is.null( inFile() )){
      if(is.null( inFile2() )){
        ""
      } else{
        temp <- gsub("_ProteinLevelData.csv.gz", " ", inFile2()$name)
        return(substr(temp, 12, nchar(temp)))
        }
    } else {
      temp <- gsub("_GroupComparisonResult.csv", " ", inFile()$name)
      return(substr(temp, 12, nchar(temp)))
      }
  })
  
  titleName <- reactive({
    if (input$experimentName == ""){
      dataName()
    } else{
      input$experimentName
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
  
  # Text 
  output$errorMsg <- renderText({
    if ( sum(unique(resultsData()$posGROUP) %in% unique(intensitiesData()$data$GROUP)) < length(unique(resultsData()$posGROUP)) ){
      "ERROR: Condition names in data files do not match! Make sure you're uploading BOTH correct files!"
    } else{ "" }
  })
  
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
  threshs <- reactive({
    if (input$volcLines == T){
      c(input$foldchange, input$fdr)
    }else{
      NULL
    }
  })
  
  output$volcanoPlot <- renderPlot({
    plotProteomicsVolcano(resultsData(), input$pcol, input$ncol, input$separate, input$volcLabel, titleName(), thresholds = threshs())
  }, width = dimensionw, height = dimensionh)
  
  output$dlvolcano <- downloadHandler(
    filename = function(){
      paste(titleName(), "_volcanoPlot_", Sys.Date() , ".pdf", sep="")
    },
    content = function(file) {
      makePdf(plotProteomicsVolcano(resultsData(), input$pcol, input$ncol, input$separate, input$volcLabel, titleName(), thresholds = threshs()), file)
    }
  )
  
  # Bar Chart Tab
  output$stackedBarPlot <- renderPlot({
    plotStackedBarChart(resultsDataPlus(), dataType(), input$pcol, input$ncol, input$ipcol, input$incol, titleName = titleName())
  }, width = dimensionw, height = dimensionh)
  
  output$dlbarchart <- downloadHandler(
    filename = paste(titleName(), "_effectBarChart_", Sys.Date(), ".pdf", sep=""),
    content = function(file) {
      makePdf(plotStackedBarChart(resultsDataPlus(), dataType(), input$pcol, input$ncol, input$ipcol, input$incol, titleName = titleName()), file)
    }
  )
  
  # PCA Tab
  output$pcaPlot <- renderPlot({
    plotPCA(intensitiesData()$int.mat, dataType(), input$pcaCircleConditions, input$pcaNoLabel, input$pcaNoLegend, titleName = titleName ())
  }, width = dimensionw, height = dimensionh)
  
  output$dlpca <- downloadHandler(
    filename = paste(titleName(), "_PCA_", Sys.Date(), ".pdf", sep=""),
    content = function(file) {
      makePdf(plotPCA(intensitiesData()$int.mat, dataType(), input$pcaCircleConditions, input$pcaNoLabel, input$pcaNoLegend, titleName = titleName ()), file)
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
    plotHeatmap(mat(), sample(), input$hmcluster, titleName() )
  }, width = dimensionw, height = dimensionh)
  
  output$dlheatmap <- downloadHandler(
    filename = paste(titleName(), "_Heatmap_", Sys.Date(), ".pdf", sep=""),
    content = function(file) {
      makePdf(plotHeatmap(mat(), sample(), input$hmcluster, titleName()), file)
    }
  )
  
}

shinyApp(ui = ui, server = server)
