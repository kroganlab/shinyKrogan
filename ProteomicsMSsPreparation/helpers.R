library (data.table)
library(ggplot2)

###
### Data Processing Functions ###
###

applyLabels <- function(data, keys){
  
  
  if(!is.null(input$keys)){
      
    keyConditionCol <- "Condition"
    keyRawCol <- "RawFile"
    keyBioRepCol <- "BioReplicate"
    
    setnames(keys, c(keyConditionCol, keyRawCol, keyBioRepCol), c("Sample", "Raw", "BioReplicate"))
    keys[, Raw :=as.character(Raw)]
    
    data [, Run := tstrsplit(data$Run, "\\.", keep = 1) ] [keys, c("BioReplicate", "Condition") := .(i.BioReplicate, i.Sample), on = .(Run = RawFile)]
  } else { 
      if(!"Condition" %in% colnames(data)){
      data[, Condition := Run]
      }
    }
  
  return(labelledData)
}


###
### Rendering Functions ###
###

intensityHistPlot <- function(inputFile){
  return(hist(log2(inputFile$Intensity)))
}

###
### DOWNLOAD FUNCTIONS
###
SaveAsPDF <- function(graphics, prefix = "", dimensions = NULL){
  #path <- PDFBackupFileName(prefix, subDir)
  path <- prefix
  
  if (is.null(dimensions))
    dimensions <- dev.size(units = "in")
  
  #print (sprintf("Writing image to  %s", path))
  cairo_pdf(path, width = dimensions[1], height = dimensions[2])
  
  # handle functions, my enrichment heatmaps that are part of a list
  if ("function" %in% class(graphics)){
    graphics()
    g <- "finished" # something to print to console instead of graphics to device 
  }else if (! ("ggplot" %in% class(graphics) | "grob" %in% class(graphics) | "Heatmap" %in% class(graphics) | "HeatmapList" %in% class(graphics))){
    g <- graphics$hmList    
  }  else{
    g <- graphics
  }
  print (g)
  
  dev.off()
}

makePdf <- function(x, file, dimensions = NULL){
  if (is.null(dimensions)){
    dimensions <- 2*dev.size(units = "in")
  } else{
    dimensions <- dimensions / 75
    }
  pdf(file = file, width = dimensions[1], height = dimensions[2])
  plot(x)
  dev.off()
}


### Scratch Space to format sample data 
#results <- fread("sample_data\\2022_12_27_Sample_Phospho_ProteinLevelData.csv")
#results[, GROUP := moreGsub(c("_A549", "MOI5", "USA", "Mock", "_2h", "_6h"), c("", "", "", "Ctl", "_02h", "_06h" ), results$GROUP) ]
#results[, SUBJECT := paste0(GROUP, ".", SUBJECT) ]
#results <- results[!grepl("THP1", GROUP)]
#fwrite( results,"sample_data\\2022_12_27_Sample_Phospho_ProteinLevelData.csv.gz" )
