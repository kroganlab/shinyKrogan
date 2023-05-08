inRstudio <- "yes"
test_directory <- Root_path

Generate_comparables_inputs_COMPASS <- function(path_to_object_inputs){
        
        if(inRstudio == "yes"){
                scriptPath <- test_directory 
     #           setwd(scriptPath)
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
        
        # ~~ COMPPASS ~~
#        if(config$comppass$enabled){
                cat(">> COMPPASS\n")
   #             source(paste(scriptPath,'/src/comppass.R',sep=""))
                #data_file = paste(config$files$output_dir,'preprocessed.txt',sep='/')
                output_dir = paste(config$files$output_dir,'COMPPASS/',sep='/')
                dir.create(output_dir, showWarnings = F)  #create Comppass directory    
                output_file = paste(output_dir, gsub('.txt', '_COMPPASS.txt', basename(config$files$data)),sep='/')
                comppass.results = Comppass.main(matrix_file, output_file, resampling=F)
 #               > Comppass.main
 #               function(data_file, output_file, resampling=F){
                        data_file <- matrix_file
                        resampling = F
                        cat("\tREADING\n")
                        data1 = read.delim(data_file,  stringsAsFactors=F, header=T,skip=1, check.names=FALSE)
                        cat("\tCONVERTING\n")
                        ## convert into intermediate format
                        data_tmp = Comppass.cleanMatrix(data1)
                        cnames = Comppass.cnames(data1)
                        data_long = Comppass.convertMatrix(data_tmp, cnames)
                        
                        cat("\tCOMPUTING STATS\n")
                        ## make stats table
                        stats_tab = Comppass.StatsTable(data_long)
                        
                        ## make reproducibility table
                        repro_tab = Comppass.ReproTable(data_long)
                        
                        ## make specificity table
                        speci_tab = Comppass.SpeciTable(stats_tab)
                        
                        cat("\tCOMPUTING SCORES\n")
                        ## compute scores
                        Z = Comppass.Z(stats_tab)
                        S = Comppass.S(stats_tab, speci_tab)
                        D = Comppass.D(stats_tab, speci_tab, repro_tab)
                        WD = Comppass.WD(stats_tab, speci_tab, repro_tab)
                        #get number of replicates
                        num_reps = as.matrix(table(as.character(names(data1)[-c(1:4)])))
                        WDv2 = Comppass.WDv2(stats_tab, speci_tab, repro_tab, num_reps)
                        
                        cat("\tSUMMARIZING\n")
                        ## compile all scores
                        summary = Comppass.Summary(stats_tab, Z, S, D, WD, WDv2)
                        summary1 <- summary
                        if(resampling){
                                print("RESAMPLING SCREEN")
                                ## resample the screen
                                P = Comppass.ResampledPvalues(data_tmp, cnames, summary)
                                
                                summary = cbind(summary, P)
                        }
                        cat("\tWRITING\n")
                        ## write out
                        write.table(summary, file=output_file, row.names=F, col.names=T, eol="\n", sep="\t", quote=F) 
#                        return(summary)
#                }              

                comppass.results = comppass.results[,c('Bait','Prey','WD')]
 #       }
        
#        COMPASS_dfs.list <- list(data, data_tmp, cnames, data_long, stats_tab, repro_tab, speci_tab, Z, S, D, WD, num_reps,
#                                 WDv2, summary1, summary, output_dir, output_file, comppass.results )
#        names(COMPASS_dfs.list) <- c("data1", "data_tmp", "cnames", "data_long", "stats_tab", "repro_tab", "speci_tab", 
#                                     "Z", "S", "D", "WD", "num_reps","WDv2", "summary1", "summary", "output_dir",
#                                     "output_file", "comppass.results")
        
#        save(COMPASS_dfs.list, file = paste0(path_to_object_inputs, "COMPASS_dfs") , ascii = F)
        save(list = c("data1", "data_tmp", "cnames", "data_long", "stats_tab", "repro_tab", "speci_tab", 
                       "Z", "S", "D", "WD", "num_reps","WDv2", "summary1", "summary", "output_dir",
                      "output_file", "comppass.results"), 
             file = paste0(path_to_object_inputs, "COMPASS_dfs") , ascii = F)
}

# path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_1/input_1_objects/"
# Generate_comparables_inputs_COMPASS(path_to_object_inputs )
# 
# path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_2/input_2_objects/"
# Generate_comparables_inputs_COMPASS(path_to_object_inputs )
# 
# path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_3/input_3_objects/"
# Generate_comparables_inputs_COMPASS(path_to_object_inputs )