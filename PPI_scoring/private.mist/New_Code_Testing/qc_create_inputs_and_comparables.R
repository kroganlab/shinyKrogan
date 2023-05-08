inRstudio <- "yes"
test_directory <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/"
#Path_to_input_object <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_1/input_1_objects/"

Generate_comparables_inputs_qc <- function(path_to_object_inputs){

if(inRstudio == "yes"){
#        scriptPath <- test_directory 
#        setwd(scriptPath)
}


load(file = paste0(path_to_object_inputs, "config" ))
#source(paste(Code_Path,"qc.R",sep=""))

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



if(config$qc$enabled){
        cat(">> QUALITY CONTROL\n")
        if(!config$preprocess$enabled){ ## use previous data matrix instead of the one from pre-processing call 
                matrix_file = config$qc$matrix_file
        }
        load(file = paste0(path_to_object_inputs, "matrix_file" ))
#        qc.main(matrix_file=matrix_file, font_scale=config$qc$cluster_font_scale, cluster=config$qc$cluster, ip_dists=config$qc$ip_distributions)
#> qc.main
#function(matrix_file, font_scale, cluster=T, ip_dists=T){
        font_scale <- config$qc$cluster_font_scale; cluster=T; ip_dists=T
        ip_matrix = read.delim(matrix_file, stringsAsFactors=F)
     
        data_matrix = qc.dataMatrix(ip_matrix)  
        ip_baits = qc.getIpToBaits(ip_matrix)
#         if(cluster){
#                 qc.clusterHeatmap(data_matrix, gsub('.txt','.pdf',matrix_file), ip_baits, font_scale)  
#         }
#         if(ip_dists){
#                 qc.ipDists(data_matrix, ip_baits, matrix_file)
#                 qc.NumUniquePlot(ip_matrix, matrix_file)
#         }
qc_dfs.list <- list(ip_matrix, data_matrix, ip_baits, font_scale)
names(qc_dfs.list) <- c("ip_matrix", "data_matrix", "ip_baits", "font_scale")

save(qc_dfs.list, file = paste0(path_to_object_inputs, "qc_df") , ascii = F)
}
}

# path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_1/input_1_objects/"
# Generate_comparables_inputs_qc(path_to_object_inputs )
# 
# path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_2/input_2_objects/"
# Generate_comparables_inputs_qc(path_to_object_inputs )
# 
# path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_3/input_3_objects/"
# Generate_comparables_inputs_qc(path_to_object_inputs )

