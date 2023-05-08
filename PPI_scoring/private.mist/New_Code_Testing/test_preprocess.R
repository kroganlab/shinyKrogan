
require(testthat)

RUN_TEST_on_preprocess <- function(Path_to_input_object ){
        load(file = paste0(Path_to_input_object, "config" ))
        load(file = paste0(Path_to_input_object, "matrix_file" ))
        load(file = paste0(Path_to_input_object, "preprocessdf" ))
        names(preprocessdf.list)
        
        # create variable names that match function calls from object (this is 
        # simply the initial call of preprocessing.main)
        
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
        
        source(paste(Path_to_new_Code,"preprocess.R",sep=""))
        
        # Now, start testing functions - they are all loaded up!
        #debug(preprocess.main)
        if(FULL_TEST) {
                cat("\tpreprocess.main\n")
                        matrix_file_new = preprocess.main(data_file=config$files$data, keys_file=config$files$keys, 
                                          output_file=paste(config$files$output_dir,'preprocessed.txt',sep='/'), 
                                          filter_data=config$preprocess$filter_contaminants, 
                                          contaminants_file=config$preprocess$contaminants_file, 
                                          rm_co=config$preprocess$remove_carryover, 
                                          collapse_file=config$files$collapse, 
                                          exclusions_file=config$files$specificity_exclusions, 
                                          remove_file=config$files$remove, id_colname=config$preprocess$id_colname, 
                                          prey_colname=config$preprocess$prey_colname, 
                                          pepcount_colname=config$preprocess$pepcount_colname, 
                                          mw_colname=config$preprocess$mw_colname)
        
 #       test_that("Processing.main", {
                cat("\tpreprocess.main.compare\n")
                expect_that(matrix_file_new, equals(matrix_file))
        }
   #             cat("\tpreprocess.remove.Duplicates\n")
        cat("\t testing function preprocess.removeDuplicates\n")
 capture.output(preprocess_remove_duplicates <-preprocess.removeDuplicates(preprocessdf.list$df.removeDuplicates, 
                                                                           id_colname, prey_colname) )
                expect_that(preprocess_remove_duplicates, 
                            equals(preprocessdf.list$df.filterContaminants))
 cat("\t testing function preprocess.filterContaminants\n")
                capture.output(preprocess_filter_contaminants <- preprocess.filterContaminants(contaminants_file, 
                                preprocessdf.list$df.filterContaminants, prey_colname),file = "output")
                expect_that(preprocess_filter_contaminants, 
                            equals(preprocessdf.list$df.mergeData))
 cat("\t testing function preprocess.mergeData\n")
                capture.output <- (preprocess_mergedata <- preprocess.mergeData(preprocessdf.list$df.mergeData, 
                                preprocessdf.list$keys, id_colname) )
                expect_that(preprocess_mergedata, equals(preprocessdf.list$df.orderExperiments))
 cat("\t testing function preprocess.findCarryover\n")
                expect_that(preprocess.findCarryover(preprocessdf.list$df.findCarryover, 
                                                     id_colname, prey_colname, pepcount_colname), 
                            equals(preprocessdf.list$to_remove))
cat("\t testing function preprocess.createMatrix\n")
                capture.output <- (preprocess_createMatrix <- preprocess.createMatrix(preprocessdf.list$df.createMatrix, 
                collapse_file, exclusions_file, remove_file, id_colname, prey_colname, pepcount_colname, mw_colname) )
                expect_that(preprocess_createMatrix, 
                            equals(preprocessdf.list$df_mat))
#        })
        
}

