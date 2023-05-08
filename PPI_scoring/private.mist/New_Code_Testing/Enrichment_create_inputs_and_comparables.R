inRstudio <- "yes"
test_directory <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/"
Path_to_input_object <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_1/input_1_objects/"

path_to_object_inputs <- Path_to_input_object
#Generate_comparables_inputs_Enrichment <- function(path_to_object_inputs){
        
        if(inRstudio == "yes"){
                scriptPath <- test_directory 
                setwd(scriptPath)
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
 #       load(file = paste0(path_to_object_inputs, "MIST_dfs" ))
 #       load(file = paste0(path_to_object_inputs, "COMPASS_dfs" ))
 load(file = paste0(path_to_object_inputs, "Annotate_dfs" ))
#        source(paste(scriptPath,'/src/annotate.R',sep=""))
        # Print out all scores
        output_file = gsub('.txt', "_ALLSCORES.txt", matrix_file)
        results = data.frame(results[,2:1], results[,3:dim(results)[2]], stringsAsFactors=F)
        results = results[order(results[,1]),]
        write.table(results, output_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
        
        # ~~ Enrichment (ORIGINAL) ~~
        if(config$enrichment$enabled){  # Enrichment analysis
                cat(">> ENRICHMENT\n")
                output_dir = paste(config$files$output_dir,'/Enrichment/',sep='/')
                
                # create Enrichment directory for files
                if(!file.exists(output_dir))
                        dir.create(output_dir)
                
                # GO HyperG enrichment analysis
                if(config$enrichment$hyperg){
                        cat("    CALCULATING HYPERG PROBABILITIES\n")
                        source(paste(scriptPath,'/src/Enrichment.R',sep=""))    
                        data_file = paste(config$files$output_dir,'preprocessed.txt',sep='/')
                        
#                        Enrichment.main(data_file, output_dir, config$preprocess$prey_colname, grouped=T, enrichment_p_cutoff=config$enrichment$enrichment_p_cutoff, id_type=config$enrichment$id_type)
                enrichment_p_cutoff=config$enrichment$enrichment_p_cutoff
                id_type=config$enrichment$id_type
                grouped=T
                prey_colname <- config$preprocess$prey_colname

## Enrichment.main = function(data_file, output_dir, prey_colname, grouped=T, enrichment_p_cutoff=enrichment_p_cutoff, id_type="UNIPROT"){
        # data_file= "~/projects/mist/tests/small/processed/preprocessed_NoC_MAT_MIST.txt"
        # data_file = "~/HPC/MSPipeline/tests/Benchmark/kinases/data/processed_ok/kinases_data_wKEYS_MAT_ALLSCORES_o.txt"
        base_name = sub("(.+)[.][^.]+$", "\\1", basename(data_file))
        df = read.delim(data_file, stringsAsFactors=F)
        prey_idx = grep(prey_colname, names(df))
        set_idx = grep('BAIT', names(df))
        
        GO_CC_groups = c()
        GO_BP_groups = c()
        GO_MF_groups = c()
        PFAM_groups = c()
        KEGG_groups = c()
        
        sets = unique(df[,set_idx])
        preys = unique(df[,prey_idx])
        
        for(i in 1:length(sets)){
                set = sets[i]
                print(set)
                
                ## GET ALL THE STUFF
                d = df[df[,set_idx]==set,]
                if(id_type == "ENTREZID"){
                        da = Enrichment.annotate(ids=unique(d[,prey_idx]), type="ENTREZID",columns=c("UNIPROT","SYMBOL","ENTREZID"))
                        if(length(is.na(da))==0){next}
                        da = da[!is.na(da$ENTREZID),]
                        da_entrez_ids = unique(da$ENTREZID)
                }else if(id_type == "UNIPROT"){
                        da = Enrichment.annotate(ids=unique(d[,prey_idx]), columns=c("UNIPROT","SYMBOL","ENTREZID"))
                        if(length(is.na(da))==0){next}
                        da = da[!is.na(da$ENTREZID),]
                        da_entrez_ids = unique(da$ENTREZID)
                }
                
                ###################
                ## ENRICHMENT WITH GO
                enrichment_p_cutoff = 1 ## set to 1 here to not get empty set error in getXXXEnrichment methods and filter later before writing out tables
                #print("GO CC enrichment")
                da_GO_CC_enriched = Enrichment.getGOEnrichment(da_entrez_ids, ontology="CC", p_cutoff=enrichment_p_cutoff)
                ## save for testing               
                da_GO_CC_enriched_save1 <- da_GO_CC_enriched
                da_entrez_ids_save <- da_entrez_ids; enrichment_p_cutoff_save <- enrichment_p_cutoff
                ## end save
                
                if(!is.null(da_GO_CC_enriched)){	#corrects for when entrenze ids arent in GO_CC
                        da_GO_CC_enriched$Pvalue = p.adjust(da_GO_CC_enriched$Pvalue, method="fdr")
                        da_GO_CC_groups = sqldf(sprintf("select '%s' as 'setid', DAE.*, UNIPROT as 'uniprot_ac', ENTREZ as 'entrez_id', SYMBOL as 'symbol' from da_GO_CC_enriched DAE join da DA on (DA.ENTREZID = DAE.ENTREZ) order by DAE.Pvalue asc",set))
                        GO_CC_groups = rbind(GO_CC_groups, da_GO_CC_groups)
                }
                ## save for testing
                da_GO_CC_enriched_save2 <- da_GO_CC_enriched
                ## end save.
                #print("GO MF enrichment")    
                da_GO_MF_enriched = Enrichment.getGOEnrichment(da_entrez_ids, ontology="MF", p_cutoff=enrichment_p_cutoff)
                
                da_GO_CC_enriched_save3 <- da_GO_CC_enriched
                
                if(!is.null(da_GO_MF_enriched)){	#corrects for when entrenze ids arent in GO_MF
                        da_GO_MF_enriched$Pvalue = p.adjust(da_GO_MF_enriched$Pvalue, method="fdr")
                        da_GO_MF_groups = sqldf(sprintf("select '%s' as 'setid', DAE.*, UNIPROT as 'uniprot_ac', ENTREZ as 'entrez_id', SYMBOL as 'symbol' from da_GO_MF_enriched DAE join da DA on (DA.ENTREZID = DAE.ENTREZ) order by DAE.Pvalue asc",set))
                        GO_MF_groups = rbind(GO_MF_groups, da_GO_MF_groups)
                }
                da_GO_CC_enriched_save4 <- da_GO_CC_enriched               
                #print("GO BP enrichment")
                da_GO_BP_enriched = Enrichment.getGOEnrichment(da_entrez_ids, ontology="BP", p_cutoff=enrichment_p_cutoff)
                if(!is.null(da_GO_BP_enriched)){	#corrects for when entrenze ids arent in GO_BP
                        da_GO_BP_enriched$Pvalue = p.adjust(da_GO_BP_enriched$Pvalue, method="fdr")
                        da_GO_BP_groups = sqldf(sprintf("select '%s' as 'setid', DAE.*, UNIPROT as 'uniprot_ac', ENTREZ as 'entrez_id', SYMBOL as 'symbol' from da_GO_BP_enriched DAE join da DA on (DA.ENTREZID = DAE.ENTREZ) order by DAE.Pvalue asc",set))
                        GO_BP_groups = rbind(GO_BP_groups, da_GO_BP_groups)
                }
                da_GO_CC_enriched_save5 <- da_GO_CC_enriched                 
                #####################
                ## GET THE PFAM STUFF
                #print("PFAM enrichment")
                
### !!!!!!! seems to break here, no da_PFAM_enriched.                  
                da_PFAM_enriched = Enrichment.getPFAMEnrichment(da_entrez_ids, p_cutoff=enrichment_p_cutoff)
                
                da_PFAM_enriched_save1 <- da_PFAM_enriched 
                
                if(!is.null(da_PFAM_enriched)){	#corrects for when entrenze ids arent in PFAM
                        da_PFAM_enriched$Pvalue = p.adjust(da_PFAM_enriched$Pvalue, method="fdr")
                        da_PFAM_groups = sqldf(sprintf("select '%s' as 'setid', DAE.*, UNIPROT as 'uniprot_ac', ENTREZ as 'entrez_id', SYMBOL as 'symbol' from da_PFAM_enriched DAE join da DA on (DA.ENTREZID = DAE.ENTREZ) order by DAE.Pvalue asc",set))
                        PFAM_groups = rbind(PFAM_groups, da_PFAM_groups)
                }
   
                da_PFAM_enriched_save2 <- da_PFAM_enriched  
                
                #####################
                ## GET THE KEGG STUFF
                #print("KEGG enrichment")
                da_KEGG_enriched = Enrichment.getKEGGEnrichment(da_entrez_ids, p_cutoff=enrichment_p_cutoff)
                if(!is.null(da_KEGG_enriched)){	#corrects for when entrenze ids arent in KEGG
                        da_KEGG_enriched$Pvalue = p.adjust(da_KEGG_enriched$Pvalue, method="fdr")
                        da_KEGG_groups = sqldf(sprintf("select '%s' as 'setid', DAE.*, UNIPROT as 'uniprot_ac', ENTREZ as 'entrez_id', SYMBOL as 'symbol' from da_KEGG_enriched DAE join da DA on (DA.ENTREZID = DAE.ENTREZ) order by DAE.Pvalue asc",set))
                        KEGG_groups = rbind(KEGG_groups, da_KEGG_groups)
                }
        }
        
        res = list(GO_CC_groups=GO_CC_groups, GO_BP_groups=GO_BP_groups, GO_MF_groups=GO_MF_groups, PFAM_groups=PFAM_groups, KEGG_groups=KEGG_groups)
        
        if(grouped){
                res = Enrichment.group(res)
                suffix = "_grouped"
        }else{
                suffix = ""
        }
        
        for(l in names(res)){
                e = res[[l]]
                e = e[e$Pvalue <= Enrichment.CUTOFF,]
                write.table(e, file=sprintf("%s/%s%s%s.txt",output_dir,base_name,l,suffix), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
        }
save(list = c("results", "results1","da_entrez_ids_save", "enrichment_p_cutoff_save"  ,"da_GO_CC_enriched_save1",
              "da_GO_CC_enriched_save2", "da_GO_CC_enriched_save3", "da_GO_CC_enriched_save4", "da_GO_CC_enriched_save5",
              "da_PFAM_enriched_save1", "da_PFAM_enriched_save2",  ), 
     file = paste0(path_to_object_inputs, "Enrichment_dfs") , ascii = F)

}                
# end Enrichment.R file.
                }
                
#                 # Perform over representation analysis based on re-sampling
#                 if(config$enrichment$resampling){
#                         cat("    CALCULATING RESAMPLED PROBABILITIES\n")
#                         source(paste(scriptPath,'/src/overRepresented.R',sep=""))
#                         data_file = output_file
#                         overRepresented.main(data_file, output_dir, config$annotate$uniprot_dir, score_name="MIST", prey_name="Prey", bait_name="Bait", score_threshold=config$enrichment$score_threshold)
#                         
#                 } 
#         }
# }

path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_1/input_1_objects/"
Generate_comparables_inputs_Enrichment(path_to_object_inputs )

path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_2/input_2_objects/"
Generate_comparables_inputs_Enrichment(path_to_object_inputs )

path_to_object_inputs <- "/Users/dsanders/Box Sync/Doug/projects/test_code/test_set_matrix/input_3/input_3_objects/"
Generate_comparables_inputs_Enrichment(path_to_object_inputs )