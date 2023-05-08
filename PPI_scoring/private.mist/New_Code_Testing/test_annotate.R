
require(testthat)


RUN_TEST_on_annotate <- function(Path_to_input_object ){
        load(file = paste0(Path_to_input_object, "config" ))
        load(file = paste0(Path_to_input_object, "matrix_file" ))
        
        data_file=config$files$data 
        keys_file=config$files$keys
        #        output_file=paste(config$files$output_dir,'preprocessed.txt',sep='/')
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
        
        # Now, start testing functions - they are all loaded up!
        #debug(preprocess.main)
        load(file = paste0(Path_to_input_object, "Annotate_dfs" ))
        source(paste(Path_to_new_Code,'annotate.R',sep=""))
        
        if(config$annotate$enabled){
                cat(">> ANNOTATING\n")
                test_that("Annotate.main", {
                        cat("\t testing function annotate.queryFile \n")
                        expect_that(annotate.queryFile(results1, config$annotate$species, config$annotate$uniprot_dir), equals(results))
                        
                })
                
        }
        
}
