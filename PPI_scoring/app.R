library(shiny)
library(artMS)
library(seqinr)
library(data.table)
library(yaml)
library(cRomppass)
library(dplyr)
# library(shinyFiles)
#MIST
my.load_evidencekey=function(evidence,key)
{
  data=evidence
  colnames(data) <- gsub('\\s','.',colnames(data))
  data <- subset(data, trimws(data$Proteins) != "") # remove white Proteins ids
  colnames(data)[grep(colnames(data), pattern="raw.file", ignore.case = TRUE)] <- "RawFile"
  data$Intensity[data$Intensity<1] <- NA
  if(!'IsotopeLabelType' %in% colnames(data)) data$IsotopeLabelType <- 'L'
  
  keys <- key
  
  # check
  unique_data <- unique(data$RawFile)
  unique_keys <- unique(keys$RawFile)
  keys_not_found = setdiff(unique_keys, unique_data)
  data_not_found = setdiff(unique_data, unique_keys)
  if (length(keys_not_found) > 0) {
    message(sprintf("Keys: %d/%d RawFiles not found: %s", length(keys_not_found), length(unique_keys), paste(keys_not_found, collapse=",")))
  }
  if (length(data_not_found) > 0) {
    message(sprintf("Data: %d/%d RawFiles not found: %s", length(data_not_found), length(unique_data), paste(data_not_found, collapse=",")))
  }
  
  # combine
  data = merge(data, keys, by=c('RawFile','IsotopeLabelType'))
  return(list(data=data,keys=keys))
}

filterMaxqData <- function(data){
  data_selected <- data[grep("CON__|REV__",data$Proteins, invert=T),]
  blank.idx <- which(data_selected$Proteins =="")
  if(length(blank.idx)>0)  data_selected = data_selected[-blank.idx,]
  return(data_selected)
}

filterData <- function(data){
  data_f = data
  #if(config$filters$protein_groups == 'remove') 
  data_f <- data_f[grep(";",data_f$Proteins, invert=T),]
  #if(config$filters$contaminants) 
  data_f <- filterMaxqData(data_f)
  msg <- sprintf("Data Filter: %d/%d (%s%%) records remained", 
                 nrow(data_f), nrow(data), round(100*nrow(data_f)/nrow(data),1))
  message(msg)
  return(data_f)
}

# ******************************************** #
# FUNCTIONS for AP-MS scoring
# ******************************************** #

# Filtering and aggregate intensity/spectral_count for each protein
## rm_itself: remove interactions of bait and itself
## fix_itself: for interactions of bait and itself, fix prey name
my.preprocessAPMS <- function(evidence,key, rm_itself = TRUE, fix_itself = TRUE) {
  # Load data/key
  #config <- my.load_config(dataset_id)
  data <- my.load_evidencekey(evidence,key)$data
  
  # Filter
  data_f <- data
  # if (!is.null(config$filters$resolve_baitgroup)) {
  #   data_f <- resolve_baitgroup(data=data_f, bait_pattern=trimws(config$filters$resolve_baitgroup))
  #   data_f <- remove_bait_contaminant(data_f, bait_pattern=trimws(config$filters$resolve_baitgroup))
  # }
  data_f <- filterData(data_f)
  colnames(data_f)[grep(pattern="ms.ms.count", x = colnames(data_f), ignore.case = TRUE)] <- 'spectral_counts'
  
  if (rm_itself) data_f <- subset(data_f, Condition != Proteins)
  if (!rm_itself & fix_itself) {
    idx <- which(data_f$Condition == data_f$Proteins)
    data_f[idx,"Proteins"] <- paste0(data_f[idx,"Proteins"], "prey")
  }
  
  # Aggregate Intensity
  data_f_agg <- aggregate(Intensity ~ TestControl+BaitName+RawFile+BioReplicate+Run+Condition+Proteins+Sequence+Charge, data=data_f, FUN = max)
  data_f_agg <- aggregate(Intensity ~ TestControl+BaitName+RawFile+BioReplicate+Run+Condition+Proteins, data=data_f_agg, FUN = sum)
  data_f_agg <- subset(data_f_agg, !is.na(Intensity))
  
  # Aggregate SPC
  data_f_spc <- aggregate(spectral_counts ~ TestControl+BaitName+RawFile+BioReplicate+Run+Condition+Proteins+Sequence+Charge,data=data_f,FUN = max)
  data_f_spc <- aggregate(spectral_counts ~ TestControl+BaitName+RawFile+BioReplicate+Run+Condition+Proteins,data=data_f_spc,FUN = sum)
  data_f_spc <- subset(data_f_spc, !is.na(spectral_counts))
  
  return( list(data_f = data_f, agg_intensity = data_f_agg, agg_spc = data_f_spc) )
}

my.MaxQToMIST <- function(evidence,key,ref_proteome_file, outdir = "/bf2/smb-data/tnguyen/projects/fluomics/tempdata", quant_variable="spc") {
  # Load and aggregate data/key
  datalist <- my.preprocessAPMS(evidence,key)
  
  quant_variable <- trimws(tolower(quant_variable))
  if (! quant_variable %in% c("spc","intensity")) stop("Please input quant_variable as 'spc' or 'intensity'")
  data_sel <- datalist$agg_spc
  data_sel$ms_spectral_counts <- data_sel$spectral_counts
  quant_col <- 'ms_spectral_counts'
  if (quant_variable == "intensity") {
    data_sel <- datalist$agg_intensity
    data_sel$ms_intensity <- data_sel$Intensity
    quant_col <- 'ms_intensity'
  }
  keysout <- unique(data_sel[,c("RawFile","Condition")])
  
  # Select columns
  data_sel$ms_unique_pep = ""
  data_sel <- data_sel[,c('RawFile','Proteins','ms_unique_pep', quant_col)]
  colnames(data_sel) = c('id','ms_uniprot_ac','ms_unique_pep', quant_col)
  
  # Uniprot annotate
  #conn <- my.bf2_connect()
  #uniprot <- dbGetQuery(conn, "select distinct Entry as ms_uniprot_ac, Length from view_uniprot")
  #dbDisconnect(conn)
  
  # library(dplyr)
  # uniprot <- distinct(uniprot,Entry,Length)
  # colnames(uniprot)[which(colnames(uniprot)=="Entry")]="ms_uniprot_ac"
  # d <- setdiff(data_sel$ms_uniprot_ac, uniprot$ms_uniprot_ac)
  # if (length(d)>0) {
  #   msg <- sprintf("These proteins are not found in uniprot db, please check: %s", paste(d, collapse=","))
  #   #stop(msg)
  #   message(msg)
  # }
  
  ref_proteome <- read.fasta(
    file = ref_proteome_file,
    seqtype = "AA",
    as.string = TRUE,
    set.attributes = TRUE,
    legacy.mode = TRUE,
    seqonly = FALSE,
    strip.desc = FALSE
  )
  p_lengths <- c()
  p_names <- c()
  for (e in ref_proteome) {
    p_lengths <- c(p_lengths, nchar(e[1]))
    p_names <- c(p_names, attr(e, 'name'))
  }
  ref_table <- data.table(names = p_names, lengths = p_lengths)
  ref_table[, uniprot_ac := gsub('([a-z,0-9,A-Z]+\\|{1})([A-Z,0-9,\\_]+)(\\|[A-Z,a-z,0-9,_]+)',
                                 '\\2',
                                 names)]
  ref_table[, uniprot_id := gsub('([a-z,0-9,A-Z]+\\|{1})([a-z,0-9,A-Z]+\\|{1})([A-Z,a-z,0-9,_]+)',
                                 '\\3',
                                 names)]
  colnames(ref_table)[which(colnames(ref_table)=="uniprot_ac")]="ms_uniprot_ac"
  d <- setdiff(data_sel$ms_uniprot_ac, ref_table$ms_uniprot_ac)
  if (length(d)>0) {
    msg <- sprintf("These proteins are not found in fasta, please check: %s", paste(d, collapse=","))
    #stop(msg)
    message(msg)
  }
  
  # Get mass
  data_sel <- base::merge(data_sel, ref_table, by = "ms_uniprot_ac")
  data_sel$Mass <- 110*data_sel$lengths
  
  # Write out
  ## data
  outdir <- paste0(outdir, "/mist/", quant_variable)
  dir.create(outdir, show=FALSE, recursive = TRUE)
  data_outfile <- sprintf("%s/mist-data.%s.txt", outdir, quant_variable)
  write.table(data_sel, file=data_outfile, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
  
  ## keys
  key_outfile <- sprintf("%s/mist-key.%s.txt", outdir, quant_variable)
  write.table(keysout, file=key_outfile, eol='\n', sep='\t', quote=F, row.names=F, col.names=F)
  
  return(list(data_file=data_outfile, keys_file=key_outfile))
}

myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

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
      downloadButton("downloadqc", "Run and Download QC"), 
      
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
      
      p(strong("RawFile"), ": The name of the RAW-file for which the mass spectral data was derived from."), 
      p(strong("IsotopeLabelType"), ": 'L' for label free experiments ('H' will be used for SILAC experiments)"), 
      p(strong("Condition"), ": The conditions names must follow these rules:"), 
      p("\t Use only letters (A - Z, both uppercase and lowercase) and numbers (0 - 9). The only special character allowed is underscore (_)."), 
      p("\t Very important: A condition name cannot begin with a number (R limitation)."), 
      p(strong("BioReplicate"), ": biological replicate number. It is based on the condition name. Use as prefix the corresponding Condition name, and add as suffix dash (-) plus the biological replicate number. For example, if condition H1N1_06H has too biological replicates, name them H1N1_06H-1 and H1N1_06H-2"), 
      p(strong("Run"), ": a unique number for all the MS runs (from 1 to the total number of raw files). It will be especially useful when having technical replicates. For example, in the table below, there are 2 technical replicates of the same biological replicate (Cal_33-1, technical replicates 1 and 2). A special case is SILAC experiments (H and L label are run simultaneously)"), 
      img(src = "keys_sample.png", height = 130, width = 380)
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      verbatimTextOutput("errormessage"), 
      # Output: Data file ----
      dataTableOutput("ppi_datatable")
      
    )
    
  )
)

# Define server logic to read selected file ----
server <- function(input, output, session) {
  options(shiny.maxRequestSize=2000*1024^2)
  
  # creating tmp dir
  tempdirectory <- reactive({
    # Change the relative directory
    relativePath <- tempfile(tmpdir = "/home/yzhou5/ShinyApps/R_shiny_for_PPI/PPI_scoring")
    # relativePath <- tempfile(tmpdir = "/Users/yzhou/Downloads/R_shiny_for_PPI/PPI_scoring")
    if(!exists(relativePath)){
      dir.create(relativePath)
    }
    return(relativePath)
  })
  
  qc_basic <- reactive({
    if(input$qc_basic){
      setwd(tempdirectory())
      a <- myTryCatch(artmsQualityControlEvidenceBasic(evidence_file = input$evidence$datapath, keys_file = input$keys$datapath, prot_exp="APMS")) 
      return(a$error)
    }
  })
  
  qc_extend <- reactive({
    if(input$qc_extend){
      setwd(tempdirectory())
      a <- myTryCatch(artmsQualityControlEvidenceExtended(evidence_file = input$evidence$datapath, keys_file = input$keys$datapath, plotPCA = F)) 
      return(a$error)
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
  
  err <- reactive({
    if(!is.null(qc_basic()))
      return(qc_basic())
    if(!is.null(qc_extend()))
      return(qc_extend())
  })
  
  output$errormessage <- renderPrint(err())
  
  
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
      setwd(tempdirectory())
      if("msspc" %in% list.files("."))
      {
        setwd("msspc")
        a <- getwd()
        system(paste("cd", a, sep = " "))
        system("SAINTexpress-spc spectral_counts-saint-interactions.txt spectral_counts-saint-preys.txt spectral_counts-saint-baits.txt")
        saint_spc <- read.table("list.txt", sep = "\t",header=T)
      }else{
        dir.create("msspc")
        setwd("msspc")
        SAINTexpress_spc()
        a <- getwd()
        system(paste("cd", a, sep = " "))
        system("SAINTexpress-spc spectral_counts-saint-interactions.txt spectral_counts-saint-preys.txt spectral_counts-saint-baits.txt")
        saint_spc <- read.table("list.txt", sep = "\t",header=T)
      }
      
      #SAINTexpress-int
      setwd(tempdirectory())
      dir.create("msint")
      setwd("msint")
      SAINTexpress_int()
      a <- getwd()
      system(paste("cd", a, sep = " "))
      system("SAINTexpress-int MS_Intensity-saint-interactions.txt MS_Intensity-saint-preys.txt MS_Intensity-saint-baits.txt")
      saint_int <- read.table("list.txt", sep = "\t",header=T)
      saint_spc_int <- merge(saint_spc, saint_int, by = c("Bait", "Prey"), all = T)
      return(saint_spc_int)
    }
  })
  
  MiST <- reactive({
    if(input$mist)
    {
      setwd(tempdirectory())
      # evidence <- vroom::vroom(input$evidence$datapath)
      # keys <- vroom::vroom(input$keys$datapath)
      evidence <- read.table(file = input$evidence$datapath, header = T, sep = "\t", stringsAsFactors = F)
      keys <- read.table(file = input$keys$datapath, header = T, sep = "\t", stringsAsFactors = F)
      colnames(keys)[which(colnames(keys)=="Raw.file")] = "RawFile"
      colnames(keys)[which(colnames(keys)=="SAINT")]="TestControl"
      keys$BaitName <- keys$Condition
      #uniprot <- read.delim(file="/Users/yzhou/Downloads/view_uniprot.txt",header=T,sep = "\t",stringsAsFactors=F)
      # library(seqinr)
      # library(data.table)
      my.MaxQToMIST(evidence,keys,
                    input$fasta$datapath, 
                    outdir = tempdirectory(), quant_variable="spc")
      mist_yaml=read_yaml(file="/home/yzhou5/ShinyApps/R_shiny_for_PPI/PPI_scoring/mist.yaml")
      # mist_yaml=read_yaml(file="/Users/yzhou/Downloads/R_shiny_for_PPI/PPI_scoring/mist.yaml")
      
      mist_yaml$files$data=paste(tempdirectory(), "/mist/spc/mist-data.spc.txt", sep = "")
      mist_yaml$files$keys=paste(tempdirectory(), "/mist/spc/mist-key.spc.txt", sep = "")
      mist_yaml$files$remove=paste(tempdirectory(), "/mist/spc/a.txt", sep = "")
      mist_yaml$files$collapse=paste(tempdirectory(), "/mist/spc/a.txt", sep = "")
      mist_yaml$files$specificity_exclusions=paste(tempdirectory(), "/mist/spc/a.txt", sep = "")
      mist_yaml$files$output_dir=paste(tempdirectory(), "/mist/spc/result", sep = "")
      
      mist_yaml$preprocess$filter_contaminants=0
      mist_yaml$preprocess$contaminants_file=NULL
      mist_yaml$preprocess$id_colname="id"
      mist_yaml$preprocess$prey_colname="ms_uniprot_ac"
      mist_yaml$preprocess$pepcount_colname="ms_spectral_counts"
      mist_yaml$preprocess$mw_colname="Mass"
      
      mist_yaml$qc$matrix_file=NULL
      
      mist_yaml$mist_yaml$matrix_file=NULL
      mist_yaml$mist_yaml$weights="fixed"
      
      mist_yaml$comppass$enabled=0
      
      mist_yaml$annotate$enabled=0
      mist_yaml$annotate$species="HUMAN"
      mist_yaml$annotate$uniprot_dir="/home/yzhou5/ShinyApps/R_shiny_for_PPI/PPI_scoring/uniprot-human-filtered-organism__Homo.tab"
      
      mist_yaml$enrichment$enabled=0
      
      write_yaml(mist_yaml,file=paste(tempdirectory(), "/mist/mist.yaml", sep = ""))
      
      system(paste("Rscript /home/yzhou5/ShinyApps/R_shiny_for_PPI/PPI_scoring/private.mist/main.R -c ", paste(tempdirectory(), "/mist/mist.yaml", sep = ""), sep = ""))
      # system(paste("Rscript /Users/yzhou/Downloads/R_shiny_for_PPI/PPI_scoring/private.mist/main.R -c ", paste(tempdirectory(), "/mist/mist.yaml", sep = ""), sep = ""))
      mist <- read.table(paste(tempdirectory(), "/mist/spc/result/preprocessed_NoC_MAT_MIST.txt", sep = ""), sep = "\t",header=T)
      return(mist)
    }
  })
  
  compPASS <- reactive({
    if(input$comppass){
      setwd(tempdirectory())
      if("msspc" %in% list.files("."))
      {
        interactions=read.table(file=paste(tempdirectory(), "/msspc/spectral_counts-saint-interactions.txt", sep = ""), header=F,sep="\t",stringsAsFactors = F)
      }else{
        # dir.create("msspc")
        # setwd("msspc")
        SAINTexpress_spc()
        interactions=read.table(file=paste(tempdirectory(), "/spectral_counts-saint-interactions.txt", sep = ""), header=F,sep="\t",stringsAsFactors = F)
      }
      interactions <- reshape2::dcast(interactions, V1 + V2 ~ V3, value.var = "V4", fill = 0)
      interactions <- reshape2::melt(interactions, id.vars = c("V1", "V2"), measure.vars = colnames(interactions)[3:ncol(interactions)])
      interactions_agg <- aggregate(interactions$value, by=list(V2 = interactions$V2, variable = interactions$variable), FUN=sum)
      interactions_agg <- interactions_agg[which(interactions_agg$x != 0), ]
      interactions_agg$V2_variable <- paste(interactions_agg$V2, interactions_agg$variable, sep = "-")
      interactions$V2_variable <- paste(interactions$V2, interactions$variable, sep = "-")
      interactions <- interactions[which(interactions$V2_variable %in% interactions_agg$V2_variable), ]
      interactions <- interactions[, -which(colnames(interactions) == "V2_variable")]
      
      colnames(interactions)=c("Replicate","Bait","Prey","Spectral.Count")
      Experiment.ID=interactions$Bait
      interactions=cbind(Experiment.ID,interactions)
      interactions$Prey <- as.character(interactions$Prey)
      
      
      SPCwStats.results.noNaN <- comppass(interactions, stats = NULL, norm.factor = 0.98)
      
      #it's just a few lines of code to do the percentile by bait analysis on the compPASS output file.
      #Note, sometimes you get NaN values in the compPASS output file, these need to be replaced with Inf prior to running these few lines of code.
      
      #SPCwStats.results.noNaN=read.table(file="/Users/yzhou/Downloads/Roche/JGZ04/RSV-P/comppass_scores.tsv",header=T,sep="\t")
      SPCwStats.results.noNaN[is.na(SPCwStats.results.noNaN)]=Inf
      # Create columns that give the Percentile for each WD and Z score
      SPCwStats.results.noNaN.wide = SPCwStats.results.noNaN %>% mutate(wd_percentile = percent_rank(WD), z_percentile = percent_rank(Z) ) %>%
        
        # group by the Baits so we can get the Percentile PER Bait for WD and Z scores
        group_by(Bait) %>%
        mutate(wd_percentile_perBait = percent_rank(WD), z_percentile_perBait = percent_rank(Z) ) %>%
        data.frame(stringsAsFactors=F)
      
      # write.table(SPCwStats.results.noNaN.wide, paste(tempdirectory(), "/spc_comppass_scores_percentile.txt",sep = ""), quote=F, row.names=F, sep='\t')
      return(SPCwStats.results.noNaN.wide)
    }
  })
  
  # # merge PPI data
  # ppi <- reactive({
  #     if(input$saint & !input$mist & !input$comppass){
  #         result <- SAINTexpress()
  #     }else if(!input$saint & input$mist & !input$comppass){
  #         result <- MiST()
  #     }else if(!input$saint & !input$mist & input$comppass){
  #         result <- compPASS()
  #     }else if(input$saint & input$mist & !input$comppass){
  #         result <- merge(SAINTexpress(), MiST(), all = T, by = c("Bait", "Prey"))
  #     }else if(input$saint & !input$mist & input$comppass){
  #         result <- merge(SAINTexpress(), compPASS(), all = T, by = c("Bait", "Prey"))
  #     }else if(!input$saint & input$mist & input$comppass){
  #         result <- merge(MiST(), compPASS(), all = T, by = c("Bait", "Prey"))
  #     }else{
  #         # result <- merge(SAINTexpress(), MiST(), all = T, by = c("Bait", "Prey"))
  #         result <- merge(merge(SAINTexpress(), MiST(), all = T, by = c("Bait", "Prey")), compPASS(), all = T, by = c("Bait", "Prey"))
  #     }
  #     return(result)
  # })
  
  # merge PPI data
  ppi <- reactive({
    result <- data.frame()
    if(input$saint){
      if(nrow(result) == 0){
        result <- SAINTexpress()
      }else{
        result <- merge(result, SAINTexpress(), by = c("Bait", "Prey"), all = T)
      }
    }
    if(input$mist){
      if(nrow(result) == 0){
        result <- MiST()
      }else{
        result <- merge(result, MiST(), by = c("Bait", "Prey"), all = T)
      }
    }
    if(input$comppass){
      if(nrow(result) == 0){
        result <- compPASS()
      }else{
        result <- merge(result, compPASS(), by = c("Bait", "Prey"), all = T)
      }
    } 
    return(result)
  })
  
  output$ppi_datatable <- renderDataTable(ppi(), options = list(scrollX = T))
  
  # Downloadable txt of selected dataset ----
  output$downloadppi <- downloadHandler(
    filename = function() {
      paste("ppi_result", ".txt", sep = "")
    },
    content = function(file) {
      write.table(ppi(), file, row.names = FALSE, sep = "\t", quote = F)
    }
  )
  
  
}

# Create Shiny app ----
shinyApp(ui, server)
