library(shiny)
library(artMS)
library(seqinr)
# library(shinyFiles)

# Define UI for data upload app ----
ui <- fluidPage(
    
    # App title ----
    titlePanel("PPI Scoring"),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            
            # Input: Select a file ----
            fileInput("evidence", "Choose Evidence File",
                      multiple = FALSE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv", 
                                 ".txt")),
            
            fileInput("keys", "Choose Keys File",
                      multiple = FALSE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv", 
                                 ".txt")),
            
            fileInput("fasta", "Choose fasta File",
                      multiple = FALSE,
                      accept = ".fasta"),
            
            # Horizontal line ----
            tags$hr(),
            
            # Input: Checkbox if file need QC ----
            checkboxInput("qc_basic", "QC Basic", FALSE),
            checkboxInput("qc_extend", "QC Extend", FALSE),
            
            # Button
            downloadButton("downloadqc", "Download QC"), 
            
            # Horizontal line ----
            tags$hr(),
            
            # PPI scoring
            
            # SAINTexpress
            checkboxInput("saint", "SAINTexpress", FALSE),
            # MiST
            checkboxInput("mist", "MiST", FALSE),
            # CompPASS
            checkboxInput("comppass", "CompPASS", FALSE),
            
            # Button
            downloadButton("downloadppi", "Download PPI"), 
            
            dataTableOutput("ppi_datatable")
            
            # # Input: Select separator ----
            # radioButtons("sep", "Separator",
            #              choices = c(Comma = ",",
            #                          Semicolon = ";",
            #                          Tab = "\t"),
            #              selected = ","),
            # 
            # # Input: Select quotes ----
            # radioButtons("quote", "Quote",
            #              choices = c(None = "",
            #                          "Double Quote" = '"',
            #                          "Single Quote" = "'"),
            #              selected = '"'),
            # 
            # # Horizontal line ----
            # tags$hr(),
            # 
            # # Input: Select number of rows to display ----
            # radioButtons("disp", "Display",
            #              choices = c(Head = "head",
            #                          All = "all"),
            #              selected = "head")
            
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
            
            # Output: Data file ----
            tableOutput("contents")
            
        )
        
    )
)

# Define server logic to read selected file ----
server <- function(input, output, session) {
    options(shiny.maxRequestSize=2000*1024^2)
    
    # #read data
    # evidence <- reactive({
    #     req(input$evidence)
    #     vroom::vroom(input$evidence$datapath)
    # })
    # keys <- reactive({
    #     req(input$keys)
    #     vroom::vroom(input$keys$datapath)
    # })
    # ref_proteome <- reactive({read.fasta(
    #     file = input$fasta,
    #     seqtype = "AA",
    #     as.string = TRUE,
    #     set.attributes = TRUE,
    #     legacy.mode = TRUE,
    #     seqonly = FALSE,
    #     strip.desc = FALSE
    # )})
    # write fasta file
    
    qc_basic <- reactive({
        if(input$qc_basic){
           artmsQualityControlEvidenceBasic(evidence_file = input$evidence$datapath, keys_file = input$keys$datapath, prot_exp="APMS") 
        }
    })
    
    qc_extend <- reactive({
        if(input$qc_extend){
           artmsQualityControlEvidenceExtended(evidence_file = input$evidence$datapath, keys_file = input$keys$datapath, plotPCA = F) 
        }
    })
    # zip QC file
    zip_file <- reactive({
        if(input$qc_basic & input$qc_extend){
        qc_basic()
        qc_extend()
        zip("qc.zip", c("qc_basic", "qc_extended")) 
    }else if(input$qc_basic & !input$qc_extend){
        qc_basic()
        zip("qc.zip", "qc_basic") 
    }else if(!input$qc_basic & input$qc_extend){
        qc_extend()
        zip("qc.zip", "qc_extended") 
    }
    })
    
    output$downloadqc <- downloadHandler(
        filename <- function() {
            paste("QC", "zip", sep=".")
        },
        
        content <- function(file) {
            zip_file()
            file.copy("qc.zip", file)
        },
        contentType = "application/zip"
    )
    
    SAINTexpress_spc <- reactive({
        artmsEvidenceToSaintExpress(evidence_file=input$evidence$datapath,
                                    keys_file=input$keys$datapath,
                                    ref_proteome_file=input$fasta$datapath,
                                    quant_variable="msspc",output_file="spectral_counts.txt", verbose = TRUE)
    })
    
    SAINTexpress_int <- reactive({
        artmsEvidenceToSaintExpress(evidence_file=input$evidence$datapath,
                                    keys_file=input$keys$datapath,
                                    ref_proteome_file=input$fasta$datapath,
                                    quant_variable="msint",output_file= "MS_Intensity.txt", verbose = TRUE)
    })
    
    SAINTexpress <- reactive({
        if(input$saint){
            #SAINTexpress-spc
            dir.create("msspc")
            setwd("msspc")
            SAINTexpress_spc()
            a <- getwd()
            system(paste("cd", a, sep = " "))
            system("SAINTexpress-spc spectral_counts-saint-interactions.txt spectral_counts-saint-preys.txt spectral_counts-saint-baits.txt")
            saint_spc <- read.table("list.txt", sep = "\t",header=T)
            
            #SAINTexpress-int
            dir.create("../msint")
            setwd("../msint")
            SAINTexpress_int()
            a <- getwd()
            system(paste("cd", a, sep = " "))
            system("SAINTexpress-int MS_Intensity-saint-interactions.txt MS_Intensity-saint-preys.txt MS_Intensity-saint-baits.txt")
            saint_int <- read.table("list.txt", sep = "\t",header=T)
            saint_spc_int <- merge(saint_spc, saint_int, by = c("Bait", "Prey"), all = T)
            return(saint_spc_int)
        }
    })
    
    # merge PPI data
    ppi <- reactive({
        if(input$saint){
            SAINTexpress()
        }
        
        
        
        
        
    })
    
    output$ppi_datatable <- renderDataTable(ppi())
    
    
    
    
    
    
    
    
    
    
    
    
    
    # output$contents <- renderTable({
    #     
    #     # input$file1 will be NULL initially. After the user selects
    #     # and uploads a file, head of that data file by default,
    #     # or all rows if selected, will be shown.
    #     
    #     req(input$file1)
    #     
    #     # when reading semicolon separated files,
    #     # having a comma separator causes `read.csv` to error
    #     tryCatch(
    #         {
    #             df <- read.csv(input$file1$datapath,
    #                            header = input$header,
    #                            sep = input$sep,
    #                            quote = input$quote)
    #         },
    #         error = function(e) {
    #             # return a safeError if a parsing error occurs
    #             stop(safeError(e))
    #         }
    #     )
    #     
    #     if(input$disp == "head") {
    #         return(head(df))
    #     }
    #     else {
    #         return(df)
    #     }
    #     
    # })
    
}

# Create Shiny app ----
shinyApp(ui, server)