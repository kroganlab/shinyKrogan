
require(testthat)

RUN_TEST_on_mist <- function(Path_to_input_object ){
        load(file = paste0(Path_to_input_object, "config" ))
        load(file = paste0(Path_to_input_object, "matrix_file" ))
        load(file = paste0(Path_to_input_object, "MIST_dfs" ))

suppressMessages( source(paste(Path_to_new_Code,"mist.R",sep="")))       
suppressMessages( source(paste(Path_to_new_Code,'training.R',sep="")))
        
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
        
        test_that("mist.main", {
                if(FULL_TEST) {
                cat("\t testing function mist.main \n")
                 expect_that(mist.main(matrix_file=matrix_file, weights=config$mist$weights, w_R=config$mist$reproducibility, 
                                      w_A=config$mist$abundance, w_S=config$mist$specificity, 
                                      training_file=config$mist$training_file, training_steps=config$mist$training_steps), 
                            equals(results_with_samples))
                }
                cat("\t testing function mist.processMatrix \n")
                expect_that(mist.processMatrix(dat1), 
                            equals(dat2))
                cat("\t testing function mist.getM3D_normalized \n")
                expect_that(mist.getM3D_normalized(dat2[[1]]), 
                            equals(m3d_norm))
                cat("\t testing function mist.getMetrics \n")
                expect_that(mist.getMetrics(m3d_norm, dat2[[2]]), 
                            equals(dat))
                cat("\t testing function mist.vectorize.R \n")
                expect_that(mist.vectorize(dat[[1]]), equals(R))
                cat("\t testing function mist.vectorize.A \n")
                expect_that(mist.vectorize(dat[[2]]), equals(A))
                cat("\t testing function mist.vectorize.S \n")
                expect_that(mist.vectorize(dat[[3]]), equals(S))
                cat("\t testing function mist.getSampleOccurences \n")
                expect_that(mist.getSampleOccurences(m3d_norm, dat2[[2]]), 
                            equals(ip_occurences))                                       
        })
        
}
