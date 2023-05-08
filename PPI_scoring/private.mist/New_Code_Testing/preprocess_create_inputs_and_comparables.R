#! /usr/bin/Rscript --vanilla
suppressMessages(library(getopt))

inRstudio <- "yes"

# set these two if doing testing (if inRstudio is set to "no" they are ingored when running via RSCRIPT.
# test_directory <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/"
# yaml_filename <- "Influenza_kshv_weights_1.yml"
# path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_1/input_1_objects/"

test_directory <- Root_path

Generate_comparables_inputs <- function(yaml_filename ,path_to_object_inputs ){

        #########################
        ## MIST MAIN FILE #######
        
        spec = matrix(c(
                'verbose', 'v', 2, "integer", "",
                'help'   , 'h', 0, "logical", "available arguments (this screen)",
                'config'  , 'c', 1, "character", "configuration file in YAML format"),
                byrow=TRUE, ncol=5)
        
        opt = getopt(spec = spec, opt = commandArgs(TRUE), command = get_Rscript_filename(), usage = FALSE, debug = FALSE)
        
        if(inRstudio == "yes"){
                opt$config <- yaml_filename 
                
        }
        
        # if help was asked for print a friendly message
        # and exit with a non-zero error code
        if ( !is.null(opt$help) ) {
                cat(getopt(spec, usage=TRUE));
                q(status=1);
        }
        
        suppressMessages(library(yaml))
        
        PIPELINE=T
        
        # set source directory
        args <- commandArgs(trailingOnly = F) 
        scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
        
        # add this line to run in R studio
        if(inRstudio == "yes"){
                scriptPath <- test_directory 
                setwd(scriptPath)
        }
        
        ## load all externeal files
        if(exists(Code_Path)){
                scriptPath <- Code_Path
        }
  #      source(paste(scriptPath,"/src/preprocess.R",sep=""))
        # source(paste(scriptPath,"/src/qc.R",sep=""))
        # source(paste(scriptPath,"/src/mist.R",sep=""))
        # source(paste(scriptPath,'/src/training.R',sep=""))
        # source(paste(scriptPath,'/src/annotate.R',sep=""))
        
        
        getConfig <- function(config_file){
                x = readLines(config_file)
                x = x[x!=""]  #remove \n\n cases (blank Lines)
                x = gsub(':  ',': ', gsub(":", ': ',x) )   # make sure there is a space between ':' and any character
                x = gsub('\t', '  ', x)
                config = paste(x, collapse='\n')
                config = yaml.load(config)
                return(config)
        }
        
        
        #main <- function(opt){
        config = tryCatch(getConfig(opt$config), error = function(e) { print("!!! Error loading the config file. Please make sure the file follows YAML format."); break} )
        
        # replace any spaces in the colname with a "." to match how R reads in headers with spaces. 
        config$preprocess$ip_colname <- gsub(' ','.',config$preprocess$ip_colname)
        config$preprocess$prey_colname <- gsub(' ','.',config$preprocess$prey_colname)
        config$preprocess$pepcount_colname <- gsub(' ','.',config$preprocess$pepcount_colname)
        config$preprocess$mw_colname <- gsub(' ','.',config$preprocess$mw_colname)
        
        ##  create an outputdir if it doesn't exist 
        #  if(is.null(config$files$output_dir) || config$files$output_dir == '') config$files$output_dir = sprintf('%s/processed/',getwd())
        #    
        #  dir.create(config$files$output_dir, showWarnings = T)
        
        # testing line get input in object form
        save(config, file = paste0(path_to_object_inputs, "config") , ascii = T)
        
        matrix_file = preprocess.main(data_file=config$files$data, keys_file=config$files$keys, 
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
        
        save(matrix_file, file = paste0(path_to_object_inputs, "matrix_file") , ascii = T)
        
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
        
        
        
        #   ## main switches between parts of the pipeline
        #   if(config$preprocess$enabled){
        #     cat(">> PREPROCESSING FILES\n")
        #     matrix_file = preprocess.main(data_file=config$files$data, keys_file=config$files$keys, output_file=paste(config$files$output_dir,'preprocessed.txt',sep='/'), filter_data=config$preprocess$filter_contaminants, contaminants_file=config$preprocess$contaminants_file, rm_co=config$preprocess$remove_carryover, collapse_file=config$files$collapse, exclusions_file=config$files$specificity_exclusions, remove_file=config$files$remove, id_colname=config$preprocess$id_colname, prey_colname=config$preprocess$prey_colname, pepcount_colname=config$preprocess$pepcount_colname, mw_colname=config$preprocess$mw_colname)  
        
        #> preprocess.main
        #function(data_file, keys_file, output_file, filter_data, contaminants_file, rm_co=T, collapse_file, exclusions_file, remove_file, id_colname, prey_colname, pepcount_colname, mw_colname){
        #  cat("\tREADING FILES\n")
        
        #debug(postprocess.main())
        
        keys=tryCatch(read.delim(keys_file, sep="\t", header=F, stringsAsFactors=FALSE), error = function(e) cat(sprintf('\tERROR reading keys from : %s\n',keys_file)))
        df=tryCatch(read.delim(data_file, sep="\t", header=T, stringsAsFactors=FALSE), error = function(e) cat(sprintf('\tERROR reading data from : %s\n',data_file)))
        preprocess.checkNames(df, id_colname, prey_colname, pepcount_colname, mw_colname) 
        df.checknames <- df
        names(keys) = c(id_colname, "BAIT")
        keys$BAIT = gsub(' ', '', keys$BAIT)
        
        # quality control
        ## TO DO GIT ISSUE #1
        if(class(df[,3])=="character"){
                cat("\t!!! CHARACTERS FOUND IN unique_pep COLUMN. CONVERTING TO NUMERIC.\n")
                df[,3] = as.numeric(df[,3])
        }
        df <- df[which(df[,3] > 0 | is.na(df[,3]) | df[,3] == ""),]   # remove ms_unique_pep <= 0
        df.removeDuplicates <- df
        df <- preprocess.removeDuplicates(df, id_colname, prey_colname)
        df.filterContaminants <- df
        #filter contaminants out
        cat("\tFILTERING COMMON CONTAMINANTS\n")
        if(filter_data == 1)
                df <- preprocess.filterContaminants(contaminants_file, df, prey_colname)
        
        df.mergeData <- df
        
        #merge keys with data
        cat("\tMERGING KEYS WITH DATA\n")
        df <- preprocess.mergeData(df, keys, id_colname)
        
        df.orderExperiments <- df
        
        df <- df[preprocess.orderExperiments(df, id_colname),]  #GENERATES WARNINGS WHEN ID# HAS CHARACTERS IN IT: FIXED
        #write.table(df, output_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T, na="")
        df.findCarryover <- df
        # Remove Carryover
        if(rm_co==1){
                to_remove <- preprocess.findCarryover(df, id_colname, prey_colname, pepcount_colname)
                if(length(to_remove)>0){
                        write.table(to_remove, gsub('.txt','_ToRemove.txt',output_file), eol="\n", sep="\t", quote=F, row.names=F, col.names=T, na="")
                        # remove carryover proteins
                        tmp = merge(df, data.frame(to_remove[,c(id_colname, prey_colname)], here=1), by=c(id_colname, prey_colname), all.x=TRUE) #get index of carryover proteins
                        df <- tmp[which(is.na(tmp$here)),-dim(tmp)[2]]  # remove the 'here' column
                        df <- df[preprocess.orderExperiments(df, id_colname),]
                        output_file = gsub('.txt','_NoC.txt',output_file)
                        #    write.table(df, output_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T, na="")
                }
        }
        
        df.createMatrix <- df
        # create matrix to be used by MiST
        cat("\tCONVERTING TO MATRIX\n")
        matrix_output_file = gsub('.txt','_MAT.txt',output_file)
        df_mat <- preprocess.createMatrix(df, collapse_file, exclusions_file, remove_file, id_colname, prey_colname, pepcount_colname, mw_colname) #return a list b/c of space padding
        
        
        
        write.table(df_mat[[1]], matrix_output_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=F, na="")
        write.table(df_mat[[2]], matrix_output_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=F, na="", append=TRUE)
        
        preprocessdf.list <- list(df.checknames, df.createMatrix, df.filterContaminants,
                                  df.findCarryover, df.mergeData, df.removeDuplicates, df.orderExperiments,
                                  df_mat, keys, to_remove)
        names(preprocessdf.list) <- c("df.checknames", "df.createMatrix", "df.filterContaminants",
                                      "df.findCarryover", "df.mergeData", "df.removeDuplicates",
                                      "df.orderExperiments", "df_mat", "keys", "to_remove")
        
        ?list  
        save(preprocessdf.list, file = paste0(path_to_object_inputs, "preprocessdf") , ascii = F)

}

# yaml_filename <- "Influenza_kshv_weights_1_names.yml"
# path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_1/input_1_objects/"
# Generate_comparables_inputs(yaml_filename ,path_to_object_inputs )
# 
# yaml_filename <- "Influenza_mist_weights_2.yml"
# path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_2/input_2_objects/"
# Generate_comparables_inputs(yaml_filename ,path_to_object_inputs )
# 
# yaml_filename <- "Influenza_trained_2_11_3.yml"
# path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_3/input_3_objects/"
# Generate_comparables_inputs(yaml_filename ,path_to_object_inputs )

