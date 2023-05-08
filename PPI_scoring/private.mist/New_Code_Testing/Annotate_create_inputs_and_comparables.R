inRstudio <- "yes"
test_directory <- Root_path

Generate_comparables_inputs_Annotate <- function(path_to_object_inputs){
        
        if(inRstudio == "yes"){
                scriptPath <- test_directory 
#                setwd(scriptPath)
        }
        
        load(file = paste0(path_to_object_inputs, "config" ))
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
        load(file = paste0(path_to_object_inputs, "MIST_dfs" ))
        load(file = paste0(path_to_object_inputs, "COMPASS_dfs" ))
 #       source(paste(scriptPath,'/src/annotate.R',sep=""))
 # DBS adde code to set return object within mist testing equal to object name in
# main.R script.  
        if (exists('results_with_samples') ){
                mist.results = results_with_samples
                mist.results1 = mist.results
        }
#####
# main.R script.  
        if (exists('comppass.results') ){
                comppass.results = comppass.results[,c('Bait','Prey','WD')]
                comppass.results1 = comppass.results
        }
#####

        # Combine all results
        if( exists('mist.results') & exists('comppass.results') ){
                results = merge(mist.results, comppass.results[,c('Bait','Prey','WD')], by=c('Bait','Prey'))
        }else if( exists('mist.results') ){
                results = mist.results
        }else if( exists('comppass.results') ){
                results = comppass.results
        }
        
        # If no results calculated this time
        if( !exists('results')  ){  
                output_dir = config$files$output_dir
                if( file.exists(paste(output_dir, "/COMPPASS/", gsub('.txt', '_COMPPASS.txt', basename(config$files$data)),sep='/')) & file.exists(gsub('.txt', "_MIST.txt", matrix_file)) ){
                        cat("\tLOADING MIST SCORES\n")
                        mist.results = read.delim(gsub('.txt', "_MIST.txt", matrix_file), sep = '\t', header=T, stringsAsFactors=F)
                        cat("\tLOADING COMPPASS SCORES\n")
                        comppass.results = read.delim( paste(output_dir, "/COMPPASS/", gsub('.txt', '_COMPPASS.txt', basename(config$files$data)),sep='/'), sep="\t", header=T, stringsAsFactors=F)
                        results = merge(mist.results, comppass.results[,c('Bait','Prey','WD')], by=c('Bait','Prey'))
                }else if( file.exists(paste(output_dir, "/COMPPASS/", gsub('.txt', '_COMPPASS.txt', basename(config$files$data)),sep='/')) ){
                        cat("\tLOADING COMPPASS SCORES\n")
                        results = read.delim(paste(output_dir, "/COMPPASS/", gsub('.txt', '_COMPPASS.txt', basename(config$files$data)),sep='/'), sep="\t", header=T, stringsAsFactors=F)
                }else if(file.exists(gsub('.txt', "_MIST.txt", matrix_file))){
                        cat("\tLOADING MIST SCORES\n")
                        results = mist.results = read.delim(gsub('.txt', "_MIST.txt", matrix_file), sep = '\t', header=T, stringsAsFactors=F)
                }else{
                        stop("NO SCORES FOUND. PLEASE ENABLE ONE OF THE SCORING OPTIONS.")  
                }
        }
 results1 <- results       
        # ~~ Annotations ~~
        if(config$annotate$enabled){
                cat(">> ANNOTATING\n")
                
                results = annotate.queryFile(results, config$annotate$species, config$annotate$uniprot_dir)
        }
        save(list = c("results", "results1"), 
             file = paste0(path_to_object_inputs, "Annotate_dfs") , ascii = F)
}

#path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_1/input_1_objects/"
#Generate_comparables_inputs_Annotate(path_to_object_inputs )

#path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_2/input_2_objects/"
#Generate_comparables_inputs_Annotate(path_to_object_inputs )

#path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_3/input_3_objects/"
#Generate_comparables_inputs_Annotate(path_to_object_inputs )
