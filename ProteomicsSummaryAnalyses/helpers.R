library (data.table)
library(ggplot2)
source ("./UniprotIDMapping.R")

###
### Data Processing Functions ###
###

formatResults <- function(results, dataType, spc){
  
  # name PhSite and Protein columns appropriately
  if (dataType != "Abundance"){
    results<- setnames(results, "Protein", "Site")
    results[, Protein := tstrsplit(Site, "_")[[1]]]
    if (spc == "OTHER"){
      results[, gene := Protein]
      results[, geneSite := Site]
    } else {
      results[, gene := multiUniprots2multiGenes(results$Protein, species = spc)]
      results[, geneSite := multiUniprotSites2multiGeneSites(results$Site, species = spc)] 
    }
    results[, dataID := Site]
    results[, dataGene := geneSite]
  } else{
    if (spc == "OTHER"){
      results[, gene := Protein]
    } else {
      results[, gene := multiUniprots2multiGenes(results$Protein, species = spc)]
      
    }
    results[, dataID := Protein]
    results[, dataGene :=  gene]
  }
  
  results <- results[, posGROUP := tstrsplit(Label, "-")[[1]]] [, negGROUP := tstrsplit(Label, "-")[[2]]]
  
  return(results)
}

threshResults <- function(results, pThresh, fcThresh){
  
  results[, effectSign := "missing"]
  results[!is.na(log2FC), effectSign := "notSig"]
  results[log2FC > fcThresh & adj.pvalue < pThresh, effectSign := "sigUp"]
  results[log2FC < -1*fcThresh & adj.pvalue < pThresh, effectSign := "sigDown"]
  
  return(results)
}

formatIntensities <- function(intensities, dataType, spc){
  # Pull out replicate numbers
  intensities[, REPLICATE := tstrsplit(SUBJECT, split = "\\.")[[2]]]
  # Convert UPIDs to gene names
  if (dataType != "Abundance"){
    intensities<- setnames(intensities, "Protein", "Site")
    intensities[, Protein := tstrsplit(Site, "_")[[1]]]
    if (spc == "OTHER"){
      intensities[, gene := Protein]
      intensities[, geneSite := Site]
    } else{
      intensities[, gene := multiUniprots2multiGenes(intensities$Protein, species = spc)]
      intensities[, geneSite := multiUniprotSites2multiGeneSites(intensities$Site, species = spc)] 
    }
    intensities[, dataID := Site]
    intensities[, dataGene := geneSite]
  } else{
    if (spc == "OTHER"){
      intensities[, gene := Protein]
      
    } else{
      intensities[, gene := multiUniprots2multiGenes(intensities$Protein, species = spc)]
    }
    intensities[, dataID := Protein]
    intensities[, dataGene :=  gene]
  }

  # Find Ctl Group for calculating normalized intensity
  # will look for anything with "control" or a shortening of "control" or "mock" first, then "veh", then "0h", then "WT"
  intensities[, controlGroup := F]
  intensities[grepl("([Cc][OoNn]{0,2}[Tt][RrOo]{0,2}[Ll])|([Mm][Oo][Cc][Kk])", GROUP), controlGroup := T] 
  if (sum(intensities$controlGroup) == 0){
    intensities[grepl("([Vv][Ee][Hh])", GROUP), controlGroup := T]
  }
  if (sum(intensities$controlGroup) == 0){
    intensities[grepl("(0[Hh])", GROUP), controlGroup := T]
  }
  if (sum(intensities$controlGroup) == 0){
    intensities[grepl("([Ww][Tt])", GROUP), controlGroup := T]
  }
  if (sum(intensities$controlGroup) == 0){
    normIntensity.mat <- NULL
  } else{
    intensities[, vsCTL := LogIntensities - mean(LogIntensities[controlGroup == T]), by = Protein]
    y <- dcast(intensities, dataGene ~ SUBJECT, fun.aggregate = mean, value.var = "vsCTL")
    normIntensity.mat <- as.matrix(y[,2:(ncol(y)-1)], rownames = y$dataGene)
  }
  
  x <- dcast(intensities, dataGene ~ SUBJECT, fun.aggregate = mean, value.var = "LogIntensities")
  intensity.mat <- as.matrix(x[,2:(ncol(x)-1)], rownames = x$dataGene)
  
  sweptIntensity.mat <- sweep (intensity.mat, 1, apply (intensity.mat, 1, median, na.rm = TRUE))
  
  intensitiesStructure <- list("data" = intensities, "int.mat" = intensity.mat, "norm.mat" = normIntensity.mat, "swept.mat" = sweptIntensity.mat, 
                               "ctlGroups" = unique(intensities[controlGroup == T, GROUP]))
  return(intensitiesStructure)
}

checkDataType <- function(results){
  if (sum(grepl("_", results[, Protein])) == 0){
    return("Abundance")
  } 
  if (sum(c("S", "T", "Y") %in% substr(tstrsplit(results[, Protein], "_", keep = 2)[[1]], 1,1)) > 0){
    return("Phospho")  
  }
  if (sum(c("K") %in% substr(tstrsplit(results[, Protein], "_", keep = 2)[[1]], 1,1)) > 0){
    return("Ubitquitin")
  } else{ stop("PTM cannot be determined, currently only Abundance, Phospho, and Ubiquitin datasets allowed.") }
  
}
  
threshResultsPlus <- function(results, intensities, pThresh, fcThresh){
  
  infProteinCheck <- intensities[!is.na(LogIntensities), .N, by = .(GROUP, dataID)][, .(dataID, N, Complete = N == max(N)), by = GROUP]
  # Merge infProteinCheck$Complete into results as posComplete and negComplete depending on pos/negGroup
  results[infProteinCheck, posComplete := i.Complete, on = .(posGROUP = GROUP, dataID = dataID)]
  results[infProteinCheck, negComplete := i.Complete, on = .(negGROUP = GROUP, dataID = dataID)]
  
  # Separate data by effect: missing, notSig, sigUp, sigDown, infiniteUp, infiniteDown
  results[,effectType := "missing"]
  results[!is.na(log2FC),effectType := "notSig"]
  results[(is.finite(log2FC) & log2FC > fcThresh & adj.pvalue < pThresh), effectType := "sigUp"]
  results[(is.finite(log2FC) & log2FC < -fcThresh & adj.pvalue < pThresh), effectType := "sigDown"]
  results[log2FC == Inf & posComplete == T, effectType := "infUp"]
  results[log2FC == -Inf & negComplete == T, effectType := "infDown"]
  return(results)
}  

filterSigIntensities <- function(intensitiesMat, results){
  return(intensitiesMat[rownames(intensitiesMat) %in%  results[!effectType %in% c("missing","notSig"), dataGene],])
}



###
### Plotting Functions ###
###

plotProteomicsVolcano <- function(results, posColor, negColor, separate, showLabel, titleName, thresholds = NULL){
  
  if (separate == T){
  volcPlotData <- results[, .(negLog10Adj.PValue=-log10(adj.pvalue), log2FC, posGROUP, effectSign, dataGene)]
  } else{
    volcPlotData <- results[, .(negLog10Adj.PValue=-log10(adj.pvalue), log2FC, effectSign, dataGene)]
  }
  
  volcPlotData[, repelLabel := dataGene]
  volcPlotData[effectSign == "notSig", repelLabel := NA ]
  volcPlotData <- volcPlotData[!is.na(log2FC)]
  d <- volcPlotData[is.finite(log2FC),]
  
  v <- ggplot (d) +
    geom_point(mapping = aes( x = log2FC, y = negLog10Adj.PValue, color = effectSign), alpha=0.5)+
    scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 20)) +
    scale_size_manual(values = c("TRUE" = 3, "FALSE" = 2))+
    theme_classic() +
    theme(legend.position = "none") +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=18,face="bold"),
          title =element_text(size=20, face='bold'),
          legend.text = element_text(size = 14),
          strip.text.x = element_text(size = 18, face = "bold"),
          panel.spacing.x = unit(15, "pt")) +
    ylab("-log10( adjusted p value )") +
    ggtitle(paste0("Volcano Plot of ", titleName, "Data"))
  
  if (showLabel == T){
    v <- v + ggrepel::geom_text_repel( aes( x = log2FC, y = negLog10Adj.PValue, color = effectSign, label = repelLabel)) 
  }
  
  if (separate == T){
    v <- v +  facet_wrap(~posGROUP)  
  }
  if (!is.null(thresholds)){     # thresholds is c()
    v <- v + geom_vline(xintercept = thresholds[1], linetype="dashed", 
                        color = posColor) +
             geom_vline(xintercept = -thresholds[1], linetype="dashed", 
                        color = negColor) +
             geom_hline(yintercept = -log10(thresholds[2]), linetype="dashed", 
                        color = "black")
  }
  
  v = v + scale_color_manual( values = c(sigUp = posColor,
                                         sigDown = negColor,
                                         notSig = "grey"))
  
  return(v)
}


plotStackedBarChart <- function(results, dataType, posColor = "red", negColor = "blue", infposColor = "firebrick", infnegColor = "navy", poscond = "Treated", titleName){
  results <- data.table(results)
  # count how many of each effect per condition and make a new DT 
  barChartData <- results[,.(Sum=.N), by=list(effectType, Label)]
  # Exclude irrelevant data
  barChartData <- barChartData[!effectType %in% c("notSig", "missing")]
  # Make legible for ggplot
  setnames(barChartData,  c("Sum","effectType", "Label"), c("Affected_Proteins", "Effect", "Contrast"))
  # Convert to factors and level
  barChartData[,Effect := factor(barChartData$Effect, levels = c("infUp","sigUp","infDown","sigDown"))]
  # Better Labels if Contrast is infected - mock
  levels(barChartData$Effect) <- list(Only_In_Treated = "infUp", Upregulated = "sigUp", Not_In_Treated = "infDown", Downregulated = "sigDown")
  
  xLabel <- "Contrast"
  if (dataType != "Abundance"){
    yLabel <- paste0("Significantly Affected", dataType, "Sites")
  } else {
    yLabel <- "Significantly Affected Proteins"
  }
  
  # Plot as a stacked bar chart w ggplot
  b<-ggplot(data= barChartData, aes(x=Contrast, y=Affected_Proteins,fill=Effect)) +
    geom_bar(stat="identity")+
    theme_classic() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=18,face="bold"),
          title =element_text(size=20, face='bold'),
          legend.text = element_text(size = 14),
          axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) + 
    ylab(yLabel) +
    ggtitle(paste0("Significant Effects in ", titleName, "Data"))
  
  # Recolor
  b = b + scale_fill_manual( values = c(Only_In_Treated = infposColor,
                                        Not_In_Treated = infnegColor,
                                        Upregulated = posColor,
                                        Downregulated = negColor))
  return(b)  
}


plotPCA <- function(intensity.mat, dataType, showEllipse, rmvLabel, rmvLegend, titleName){
  complete.mat <- intensity.mat[complete.cases(intensity.mat),]
  # Run principal component analysis on the transposed (flipped rows and cols) matrix
  pcaOut <- prcomp(t(complete.mat))
  pcaDT <- as.data.table(pcaOut$x, keep.rownames=TRUE)
  setorder(pcaDT, rn)
  pcaDT[, Condition := tstrsplit(rn, "\\.")[[1]]]
  pcaDT[, Replicate := tstrsplit(rn, "\\.")[[2]]]
  
  pcaPercentVar <- round(100 * (pcaOut$sdev^2)/sum(pcaOut$sdev^2), 1)
  
  o <-  ggplot (pcaDT, aes(x=PC1, y=PC2, color = Condition, shape = Replicate, label = rn)) + 
    geom_point(alpha=1.0, size=6) + 
    theme_bw() + 
    xlab (sprintf ("PC1, %.1f%%", pcaPercentVar[1])) + 
    ylab (sprintf ("PC2, %.1f%%", pcaPercentVar[2])) +
    ggtitle(paste0("PCA of Individual Runs in ", titleName, "Data")) +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=18,face="bold"),
          title =element_text(size=20, face='bold'),
          legend.text = element_text(size = 14))#+ 
    #ggtitle (sprintf ("PCA using %d phosphosites (log intensity)", nrow(complete.mat))) 
  
  if (showEllipse){
    o <-  ggplot (pcaDT, aes(x=PC1, y=PC2, color = Condition, label = rn)) + 
      geom_point(alpha=1.0, size=6) + 
      theme_bw() + 
      xlab (sprintf ("PC1, %.1f%%", pcaPercentVar[1])) + 
      ylab (sprintf ("PC2, %.1f%%", pcaPercentVar[2])) +
      ggtitle(paste0("PCA of Individual Runs in", titleName, "Data")) +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=18,face="bold"),
            title =element_text(size=20, face='bold'),
            legend.text = element_text(size = 14))
    o <- o + ggforce::geom_mark_ellipse(aes(fill = Condition, color = Condition, label = Condition), label.fontsize = 14)
  }
  if(!rmvLabel){
    o <- o + ggrepel::geom_text_repel(nudge_x = 4, nudge_y = 1, size = 6)
  }
  if(rmvLegend){
    o <- o + guides(color = "none", shape = "none")
    if (showEllipse){
      o <- o + guides(fill = "none")
    }
  }
  
  
  return(o)
}


rowClusterWithNA <- function(mat, corr = FALSE, na.value = 0, ...){
  # corr bool asserts whether to use euclidean distance or correlation 
  # euclidean distance is more consistent for similar magnitudes, but can overlook patterns in proteins with very different absolute magnitudes but similar kinetics
  if (corr){
    mat[is.na(mat)] <- na.value
    dst <- as.dist(1-cor(t(mat), use = "pairwise"))
    dst[is.na(dst)|!is.finite(dst)] <- na.value
    hclust(dst, method = "ward.D2", ...)
  } else {
    mat[is.na(mat)] <- na.value
    dst<-dist(mat)
    dst[is.na(dst)|!is.finite(dst)] <- na.value
    hclust(dst, method = "ward.D2", ...)
  }
}


plotHeatmap <- function(mat, sampleSize, clusterBool, titleName){
  if (!is.null(sampleSize) && nrow(mat) > sampleSize){
    mat <- mat[sample(nrow(mat), sampleSize, replace = FALSE), ] 
  }
  #title <- "Intensity (Normalized to Ctrl) - Per Run View"
  #filepre <- "Heatmap_SweptIntensity_"
  hclustObject <- rowClusterWithNA(mat)
  
  h <- Heatmap(mat, 
               cluster_rows = hclustObject,
               name = "Log Intensity",
               row_title = NULL,
               column_gap = unit(2, "mm"),
               cluster_columns = clusterBool,
               column_split = colnames(mat),
               column_names_gp = gpar(fontsize = 14),
               show_row_names = FALSE,
               column_title = paste0("Heatmap of Individual Runs in ", titleName, "Data"),
               column_title_gp = gpar(fontsize = 20, fontface = "bold"),
               heatmap_legend_param = list(title = "Log Intensity\n", 
                                           title_gp = gpar(fontsize = 14, fontface = "bold"),
                                           labels_gp = gpar(fontsize = 14),
                                           legend_height = unit(50, "mm"),
                                           param = list(title_gap = unit(10, "mm")) )) 
  return(h)
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
