
require(testthat)

RUN_TEST_on_compass <- function(Path_to_input_object ){
        load(file = paste0(Path_to_input_object, "config" ))
        load(file = paste0(Path_to_input_object, "matrix_file" ))
        load(file = paste0(Path_to_input_object, "COMPASS_dfs" ))
        
        source(paste(Path_to_new_Code,'comppass.R',sep=""))
        
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
        
  #      test_that("compass.main", {
                if(FULL_TEST) {
                cat("\t testing function Comppass.main \n")
                expect_that(Comppass.main(matrix_file, output_file, resampling=F), 
                           equals(summary))
                }
                cat("\t testing function Comppass.cleanMatrix \n")
                expect_that(Comppass.cleanMatrix(data1), equals(data_tmp ))
                cat("\t testing function Comppass.cnames \n")
                expect_that(Comppass.cnames(data1), equals(cnames ))
                cat("\t testing function Comppass.convertMatrix \n")
                expect_that(Comppass.convertMatrix(data_tmp, cnames), equals(data_long ))
                cat("\t testing function Comppass.StatsTable \n")
                expect_that(Comppass.StatsTable(data_long), equals( stats_tab))
                cat("\t testing function Comppass.ReproTable \n")
                expect_that(Comppass.ReproTable(data_long), equals( repro_tab))
                cat("\t testing function Comppass.SpeciTable \n")
                expect_that(Comppass.SpeciTable(stats_tab), equals(speci_tab))
                cat("\t testing function Comppass.Z \n")
                expect_that(Comppass.Z(stats_tab), equals(Z))
                cat("\t testing function Comppass.S \n")
                expect_that(Comppass.S(stats_tab, speci_tab), equals(S))
#                tryCatch(expr, error = cat("error"))
                cat("\t testing function Comppass.D \n")
                expect_that(Comppass.D(stats_tab, speci_tab, repro_tab), equals(D))
                cat("\t testing function Comppass.WD \n")
                expect_that(Comppass.WD(stats_tab, speci_tab, repro_tab), equals(WD))
                cat("\t testing function Comppass.WDv2 \n")
                expect_that(Comppass.WDv2(stats_tab, speci_tab, repro_tab, num_reps), equals(WDv2))
                cat("\t testing function Comppass.Summary \n")
                expect_that(Comppass.Summary(stats_tab, Z, S, D, WD, WDv2), equals(summary1))                
#        })
        
}
