inRstudio <- "yes"
test_directory <- Root_path

Generate_comparables_inputs_MIST <- function(path_to_object_inputs){
        
        if(inRstudio == "yes"){
                scriptPath <- test_directory 
      #          setwd(scriptPath)
        }
        
        load(file = paste0(path_to_object_inputs, "config" ))
#        source(paste(Code_Path,"mist.R",sep=""))
#        source(paste(Code_Path,'training.R',sep=""))        
        data_file=config$files$data 
        keys_file=config$files$keys
        output_file=paste(config$files$output_dir,'preprocessed.txt',sep='/')
        filter_data=config$preprocess$filter_contaminants
        contaminants_file=config$preprocess$contaminants_file
        rm_co=config$preprocess$remove_carryover
        collapse_file=config$files$collapse
        exclusions_file=config$files$specificity_exclusions
        remove_file=config$files$remove
        id_colname=config$preprocess$id_colname
        prey_colname=config$preprocess$prey_colname
        pepcount_colname=config$preprocess$pepcount_colname
        mw_colname=config$preprocess$mw_colname

        load(file = paste0(path_to_object_inputs, "matrix_file" ))
        
#        if(config$mist$enabled){
                cat(">> MIST\n")
                if(!config$preprocess$enabled){ ## use previous data matrix instead of the one from pre-processing call 
                        matrix_file = config$mist$matrix_file
                }
                mist.results = mist.main(matrix_file=matrix_file, weights=config$mist$weights, w_R=config$mist$reproducibility, w_A=config$mist$abundance, w_S=config$mist$specificity, training_file=config$mist$training_file, training_steps=config$mist$training_steps)
                
## Pick apart Mist.main 

#> mist.main
#function(matrix_file, weights='fixed', w_R=0.30853, w_A=0.00596, w_S=0.68551, training_file, training_steps=0.1){
## create inputs for mist.main
        #weights='fixed'; w_R=0.30853; w_A=0.00596;  w_S=0.68551;  training_file;  training_steps=0.1
        weights=config$mist$weights; w_R=config$mist$reproducibility; w_A=config$mist$abundance; w_S=config$mist$specificity;
        training_file=config$mist$training_file; training_steps=config$mist$training_steps
        dat <- read.delim(matrix_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
        dat1 <- dat
        dat <- mist.processMatrix(dat)
        dat2 <- dat
        m3d_norm <- mist.getM3D_normalized(dat[[1]])
        info = dat[[2]]
        
        dat <- mist.getMetrics(m3d_norm, info)
        R <- mist.vectorize(dat[[1]])
        A <- mist.vectorize(dat[[2]])
        S <- mist.vectorize(dat[[3]])
        
        metrics = data.frame(Bait=A$Bait,Prey=A$Prey,Abundance=A$Xscore,Reproducibility=R$Xscore,Specificity=S$Xscore)
        ## only retain non-zero results
        metrics = metrics[metrics$Abundance>0,]
        ## for debug purposes
        #   output_file <- gsub('.txt', "_MIST_METRICS.txt", matrix_file)
        #   write.table(metrics, output_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t" )
        
        if(weights == 'fixed'){
                mist_scores = metrics$Reproducibility*w_R + metrics$Abundance*w_A + metrics$Specificity*w_S
                results = data.frame(metrics, MIST=mist_scores)  
        }else if(weights == 'PCA'){
                results <- mist.doPCA(R,A,S)
        }else if(weights == 'training'){
                training_set = read.delim(training_file, header=F, stringsAsFactors=F)
                colnames(training_set) = c('Bait','Prey')
                output_file <- paste(dirname(matrix_file), "prediction_rates.txt", sep="/")
                training_weights = mist.train.main(metrics, training_set, output_file, training_steps)
                cat(sprintf("\tWEIGHTS BASED ON TRAINING SET:\n\t  REPRODUCIBILITY: %s\n\t  ABUNDANCE: %s\n\t  SPECIFICITY: %s\n",training_weights$R, training_weights$A, training_weights$S))
                mist_scores = metrics$Reproducibility*training_weights$R + metrics$Abundance*training_weights$A + metrics$Specificity*training_weights$S
                results = data.frame(metrics, MIST=mist_scores)  
        }else{
                print(sprintf('unrecognized MIST option: %s',weights))
        }
        
        ## per bait get all IPs the prey was found in
        ip_occurences = mist.getSampleOccurences(m3d_norm, info)
        
        if(nrow(ip_occurences) != nrow(results)){
                print(" ERROR : INCONSISTENCY BETWEEN SCORE MATRIX AND IP OCCURENCE MATRIX ")
        }
        
        results_with_samples = merge(results, ip_occurences, by=c('Bait','Prey'))
        results_with_samples

#MIST_dfs.list <- list(dat, dat1, dat2, m3d_norm, R, A, S, results, ip_occurences, results_with_samples)
#names(MIST_dfs.list) <- c("dat", "dat1", "dat2", "m3d_norm", "R", "A", "S", "results", 
#                          "ip_occurences", "results_with_samples")

save(list = c("dat", "dat1", "dat2", "m3d_norm", "R", "A", "S", "results", 
              "ip_occurences", "results_with_samples"), file = paste0(path_to_object_inputs, "MIST_dfs") , ascii = F)
}

# path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_1/input_1_objects/"
# Generate_comparables_inputs_MIST(path_to_object_inputs )
# 
# path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_2/input_2_objects/"
# Generate_comparables_inputs_MIST(path_to_object_inputs )
# 
# path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_3/input_3_objects/"
# Generate_comparables_inputs_MIST(path_to_object_inputs )