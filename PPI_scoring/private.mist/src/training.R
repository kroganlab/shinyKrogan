#! /usr/bin/Rscript --vanilla --default-packages=utils

suppressMessages(library(pheatmap))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(MESS))
suppressMessages(library(pROC))
suppressMessages(library(precrec))

#######################################

## parameter grid 
mist.train.getParamGrid = function(steps=0.05) {
  param_weight_detail = seq(from=0,to=1, by=steps)
  param_grid = expand.grid(param_weight_detail, param_weight_detail, param_weight_detail)
  param_grid = param_grid[apply(param_grid[,1:3],1,sum) == 1,]
  colnames(param_grid) = c("Reproducibility","Abundance","Specificity")
  param_grid = data.frame(param_grid, ID=paste(param_grid$Reproducibility, param_grid$Abundance, param_grid$Specificity,sep="|"))
  param_grid
}

mist.train.getPredictionRates = function(predfull, true_negative_size, thresholds = seq(from=0,to=0.99, by=.01)){
  # When the true negative data set is specified
  if(true_negative_size > 0){
    predpos <- predfull[predfull$label==1,]
    predneg <- predfull[predfull$label==0,]
    # CHECKING THAT THERE ARE MORE THAN ZERO TRAINING PROTEINS
    if(dim(predpos)[1] == 0){
      cat("\n\nERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nTHERE IS NO TRUE POSITIVES (IS THE training_file EMPTY or what?\n\n")
      stop("Please, come back when you are ready\n\n")
    }
    
    # CHECKING THAT THE TRUE NEGATIVES SPECIFIED ARE NOT LARGER THAN THE TOTAL NUMBER OF SCORES AVAILABLE..
    if(true_negative_size > dim(predneg)[1]){
      cat("\n\nERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nTHE true_negative_size VALUE SELECTED IS WAY TOO LARGE (LARGER THAN THE TOTAL NUMBER OF SCORES AVAILABLE)\n\n")
      stop("PLEASE, COME BACK WHEN YOU ARE READY\n\n")
    }else{
      set.seed(1)
      sampledNegatives <- predneg[sample(nrow(predneg), true_negative_size),]
      prediction <- rbind(predpos, sampledNegatives)    
    }
    
    cat("\tTRAINING SET, ")
    cat("TRUE POSITIVES:", dim(predpos)[1])
    cat(" - TRUE NEGATIVES:", dim(prediction)[1] - dim(predpos)[1])
    cat(" - TOTAL:", dim(prediction)[1],"\n")
    
  }else if(true_negative_size < 0){
    cat("\n\nERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nSERIOUSLY? A NEGATIVE true_negative_size? ARE YOU KIDDING ME? IF SO... YOU ARE NOT FUNNY\n\n")
    stop("PLEASE, COME BACK WHEN YOU GET SERIOUS ABOUT USING THE MIST PIPELINE\n\n")
  }else{
    prediction <- predfull
  }

  pos = nrow(prediction[prediction$label==1,])
  neg = nrow(prediction[prediction$label==0,])
  pos_false = c()
  pos_true = c()
  neg_true = c()
  neg_false = c()
  ID=unique(prediction$ID)
  for(t in thresholds){
    pos_predicted = prediction[which(prediction$predicted_scores >= t) ,]
    neg_predicted = prediction[which(prediction$predicted_scores < t) ,]
    pf = nrow(pos_predicted[which(pos_predicted$label==0),])  # false positives
    pos_false = c(pos_false, pf)
    pos_true = c(pos_true, nrow(pos_predicted)-pf)  # true positives
    nf = nrow(neg_predicted[which(neg_predicted$label==1),])  # false negaitves
    neg_false = c(neg_false, nf)
    neg_true = c(neg_true, nrow(neg_predicted)-nf)  # true negatives
  }
  
  # Figure out at what point the Specificity=0.95
  x.roc <- roc(prediction$label, prediction$predicted_scores)
  TNR=0.95 #acceptable true negative rate
  idx <- x.roc$sensitivities[ which( abs(x.roc$specificities- TNR)==min(abs(x.roc$specificities - TNR)))  ][1]
  # get the total AUC and the partial AUC where we have a specificity=0.95
  auc.partial <- pROC::auc(x.roc, partial.auc=c(idx, 0), partial.auc.focus="se")[1] 
  auc.all <- pROC::auc(prediction$label, prediction$predicted_scores)[1]
  
  # calculate slope of the first 0.95. If multiple, pick the min
  y= data.frame(sensitivities =rev(x.roc$sensitivities), specificities = rev(1-x.roc$specificities) )
  idx.95 <- which( abs(y$specificities-(1-.95)) == min(abs(y$specificities-(1-.95))) )
  slope.95 <- min( ( y$sensitivities[idx.95] - y$sensitivities[1])/(y$specificities[idx.95] - y$specificities[1]) )
  idx.75 <- which( abs(y$specificities-(1-.75)) == min(abs(y$specificities-(1-.75))) )
  slope.75_95 <- min( ( y$sensitivities[idx.75] - y$sensitivities[idx.95])/(y$specificities[idx.75] - y$specificities[idx.95]) )
  
  # Calculate AUC for precision recall plot
  ## Generate an sscurve object that contains ROC and Precision-Recall curves
  sscurves <- evalmod(scores = prediction$predicted_scores, labels = prediction$label)
  ## Shows AUCs
  pr.auc = precrec::auc(sscurves)
  pr.auc = pr.auc[which(pr.auc$curvetypes == "PRC"),'aucs']
  
  # combine all the parts together
  res = data.frame(R_A_S=ID, threshold=thresholds, tpr=pos_true/pos, fpr=pos_false/neg, specificity=neg_true/neg, precision=pos_true/(pos_true+pos_false), fdr=pos_false/(pos_false+pos_true), acc=(pos_true+neg_true)/(pos+neg), f1=(2*pos_true)/((2*pos_true)+pos_false+neg_false), total_pos=pos_true+pos_false, ROC_auc = auc.all, ROC_auc_partial=auc.partial, slope_TFR95 = slope.95, slope_TFR75_95=slope.75_95, stringsAsFactors=F, PR_auc = pr.auc)    
  res
}


# create a ROC plot based on the MIST scores generated by specific weights
mist.plotROC <- function(output_file, param_set, prediction){
  out_file <- paste0(dirname(output_file),"/plots/ROC_plot_", paste0("R", param_set$Reproducibility, "_A",param_set$Abundance, "_S",param_set$Specificity),".pdf")
  plot_title <- sprintf('ROC based on following weights : R=%s, A=%s, S=%s',param_set$Reproducibility,param_set$Abundance ,param_set$Specificity)
  # plot ROC curve
  pdf(out_file)
  x.roc <- roc(prediction$label, prediction$predicted_scores)
  
  # With more options:  thresh=sum(specificity, sensitivity)
  TNR=0.95 #acceptable true negative rate
  idx <- x.roc$sensitivities[which( abs(x.roc$specificities- TNR)==min(abs(x.roc$specificities - TNR)))   ][1]
  plot.roc(x.roc, print.auc=TRUE, auc.polygon=TRUE, partial.auc=c(idx, 0), 
           partial.auc.focus="se", grid=c(0.1, 0.2), 
           max.auc.polygon=TRUE, print.thres=TRUE, main = plot_title, reuse.auc=FALSE)
  #   plot.roc(x.roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
  #            max.auc.polygon=TRUE, main = plot_title)

  # Add a smoothed ROC:
  plot.roc(smooth(x.roc), add=TRUE, col="blue")
  dev.off()
  
  
  # Create PRECISION-RECALL plot
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  out_file <- paste0(dirname(output_file),"/plots/PR_plot_", paste0("R", param_set$Reproducibility, "_A",param_set$Abundance, "_S",param_set$Specificity),".pdf")
  
  ## Generate an sscurve object that contains ROC and Precision-Recall curves
  sscurves <- evalmod(scores = prediction$predicted_scores, labels = prediction$label)
  ## Shows AUCs
  pr.auc = precrec::auc(sscurves)
  pr.auc = round( pr.auc[which(pr.auc$curvetypes == "PRC"),'aucs'], 4 )
  # create plot title with PR AUC in it
  plot_title <- sprintf('Precision-Recall Plot. Weights : R=%s, A=%s, S=%s, AUC=%s',param_set$Reproducibility,param_set$Abundance ,param_set$Specificity, pr.auc)
  
  # convert to ggplot2 format
  ssdf <- fortify(sscurves)
  ## Plot a ROC curve
  p_roc <- ggplot(subset(ssdf, curvetype == "PRC"), aes(x = x, y = y)) + xlab("Recall") + ylab("Precision")
  p_roc <- p_roc + geom_line() + xlim(0,1) + ylim(0,1) + theme_bw() + ggtitle(plot_title)
  ggsave(out_file, p_roc, width=7, heigh=7)
  

  
  
  # output ROC data to use elsewhere
#   roc_out <- data.frame( specificity = (1-x.roc$specificities), sensitivity=x.roc$sensitivities, RAS=param_set$ID)
#   write.table(roc_out, gsub(".txt","_ROC_metrics.txt", output_file), quote=F, row.names=F, col.names=F, sep='\t', append=T )
}


## generic grid search 
mist.train.gridSearch = function(metrics, param_grid, true_negative_size, predictionFun=mist.train.getPredictionRates, output_file, plot_roc){
  if(true_negative_size == 0 ){
    cat("\n\nUSING THE ENTIRE DATA SET AS TRUE NEGATIVES\nIF YOU WANT TO SPECIFY THE NUMBER OF TRUE NEGATIVES, THEN CHANGE THE <true_negative_size> VALUE\n\n")
  }
  prediction_rates = c()
  
  # creat plot directory
  if(!file.exists( paste0(dirname(output_file),"/plots/") )){
    dir.create(paste0(dirname(output_file),"/plots/"))
  }
  
  for(i in 1:nrow(param_grid)){
    param_set = param_grid[i,]
    print(sprintf('%s/%s : R %s A %s S %s',i,nrow(param_grid),param_set$Reproducibility,param_set$Abundance ,param_set$Specificity))
    predicted_scores = (metrics$Reproducibility * param_set$Reproducibility) + (metrics$Abundance * param_set$Abundance) + (metrics$Specificity * param_set$Specificity)
    prediction = data.frame(metrics, param_set, predicted_scores=predicted_scores)
    prediction_rates = rbind(prediction_rates, predictionFun(prediction, true_negative_size))
    if(plot_roc){
      mist.plotROC(output_file, param_set, prediction)
    }
  }
  prediction_rates
}

mist.train.label = function(metrics, training_set){
  training_set = data.frame(training_set, label=1)
  metrics = merge(metrics, training_set, by=c('Bait','Prey'), all.x=T)
  metrics[which(is.na(metrics$label)),]$label=0
  cat(">>  There were ", sum(metrics$label),"/",dim(training_set)[1], " Gold interactions detected in the data.\n")
  metrics
}

mist.train.main = function(metrics, training_set, output_file, training_steps=0.1, true_negative_size, plot_roc=1){
  metrics = mist.train.label(metrics, training_set)  
  param_grid = mist.train.getParamGrid(training_steps)
  prediction_rates = mist.train.gridSearch(metrics, param_grid, true_negative_size, output_file=output_file, plot_roc=plot_roc)
  #plotROC(prediction_rates)
  write.table(prediction_rates, file=output_file, sep="\t", quote=F, eol="\n", row.names=F)
  prediction_rates = prediction_rates[order(prediction_rates$f1, decreasing=T),]
  
  prediction_weights = as.numeric(unlist(strsplit(as.character(prediction_rates$R_A_S[1]),"\\|")))
  prediction_weights = data.frame(R=prediction_weights[1],A=prediction_weights[2],S=prediction_weights[3], stringsAsFactors=F)
  prediction_weights$threshold = prediction_rates$threshold[1]
  return(prediction_weights)
}

# metrics = read.delim('tests/entero/processed/preprocessed_MAT_MIST_METRICS.txt',stringsAsFactors=F)
# training_set = read.delim('tests/entero/input/PV_goldset.txt',stringsAsFactors=F, header=F)
# colnames(training_set) = c('Bait','Prey')
# mist.train.main(metrics, training_set)
# 
# 
