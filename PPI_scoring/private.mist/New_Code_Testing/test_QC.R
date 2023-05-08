
require(testthat)



RUN_TEST_on_qc <- function(Path_to_input_object ){
        load(file = paste0(Path_to_input_object, "config" ))
        load(file = paste0(Path_to_input_object, "matrix_file" ))
        load(file = paste0(Path_to_input_object, "qc_df" ))
        source(paste(Path_to_new_Code,"qc.R",sep=""))
        
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
        
        # Now, start testing functions - they are all loaded up!
        #debug(preprocess.main)

        test_that("QC.main", {
                if(FULL_TEST) {
                        cat("\t testing qc.read.matrix.file\n")
                         expect_that(qc_dfs.list$ip_matrix, equals(read.delim(matrix_file, stringsAsFactors=F)))
                }
                cat("\t testing function qc.dataMatrix \n")
                expect_that(qc_dfs.list$data_matrix, equals(qc.dataMatrix(qc_dfs.list$ip_matrix)))
                cat("\t testing function qc.getIpToBaits \n")
                expect_that(qc_dfs.list$ip_baits, equals(qc.getIpToBaits(qc_dfs.list$ip_matrix)))
        })
        
}
