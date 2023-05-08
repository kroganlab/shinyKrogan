my.load_evidencekey=function(evidence,key)
{
  data=evidence
  colnames(data) <- gsub('\\s','.',colnames(data))
  data <- subset(data, trimws(data$Proteins) != "") # remove white Proteins ids
  colnames(data)[grep(colnames(data), pattern="raw.file", ignore.case = TRUE)] <- "RawFile"
  data$Intensity[data$Intensity<1] <- NA
  if(!'IsotopeLabelType' %in% colnames(data)) data$IsotopeLabelType <- 'L'
  
  keys <- key
  
  # check
  unique_data <- unique(data$RawFile)
  unique_keys <- unique(keys$RawFile)
  keys_not_found = setdiff(unique_keys, unique_data)
  data_not_found = setdiff(unique_data, unique_keys)
  if (length(keys_not_found) > 0) {
    message(sprintf("Keys: %d/%d RawFiles not found: %s", length(keys_not_found), length(unique_keys), paste(keys_not_found, collapse=",")))
  }
  if (length(data_not_found) > 0) {
    message(sprintf("Data: %d/%d RawFiles not found: %s", length(data_not_found), length(unique_data), paste(data_not_found, collapse=",")))
  }
  
  # combine
  data = merge(data, keys, by=c('RawFile','IsotopeLabelType'))
  return(list(data=data,keys=keys))
}

filterMaxqData <- function(data){
  data_selected <- data[grep("CON__|REV__",data$Proteins, invert=T),]
  blank.idx <- which(data_selected$Proteins =="")
  if(length(blank.idx)>0)  data_selected = data_selected[-blank.idx,]
  return(data_selected)
}

filterData <- function(data){
  data_f = data
  #if(config$filters$protein_groups == 'remove') 
  data_f <- data_f[grep(";",data_f$Proteins, invert=T),]
  #if(config$filters$contaminants) 
  data_f <- filterMaxqData(data_f)
  msg <- sprintf("Data Filter: %d/%d (%s%%) records remained", 
                 nrow(data_f), nrow(data), round(100*nrow(data_f)/nrow(data),1))
  message(msg)
  return(data_f)
}

# ******************************************** #
# FUNCTIONS for AP-MS scoring
# ******************************************** #

# Filtering and aggregate intensity/spectral_count for each protein
## rm_itself: remove interactions of bait and itself
## fix_itself: for interactions of bait and itself, fix prey name
my.preprocessAPMS <- function(evidence,key, rm_itself = TRUE, fix_itself = TRUE) {
  # Load data/key
  #config <- my.load_config(dataset_id)
  data <- my.load_evidencekey(evidence,key)$data
  
  # Filter
  data_f <- data
  # if (!is.null(config$filters$resolve_baitgroup)) {
  #   data_f <- resolve_baitgroup(data=data_f, bait_pattern=trimws(config$filters$resolve_baitgroup))
  #   data_f <- remove_bait_contaminant(data_f, bait_pattern=trimws(config$filters$resolve_baitgroup))
  # }
  data_f <- filterData(data_f)
  colnames(data_f)[grep(pattern="ms.ms.count", x = colnames(data_f), ignore.case = TRUE)] <- 'spectral_counts'
  
  if (rm_itself) data_f <- subset(data_f, Condition != Proteins)
  if (!rm_itself & fix_itself) {
    idx <- which(data_f$Condition == data_f$Proteins)
    data_f[idx,"Proteins"] <- paste0(data_f[idx,"Proteins"], "prey")
  }
  
  # Aggregate Intensity
  data_f_agg <- aggregate(Intensity ~ TestControl+BaitName+RawFile+BioReplicate+Run+Condition+Proteins+Sequence+Charge, data=data_f, FUN = max)
  data_f_agg <- aggregate(Intensity ~ TestControl+BaitName+RawFile+BioReplicate+Run+Condition+Proteins, data=data_f_agg, FUN = sum)
  data_f_agg <- subset(data_f_agg, !is.na(Intensity))
  
  # Aggregate SPC
  data_f_spc <- aggregate(spectral_counts ~ TestControl+BaitName+RawFile+BioReplicate+Run+Condition+Proteins+Sequence+Charge,data=data_f,FUN = max)
  data_f_spc <- aggregate(spectral_counts ~ TestControl+BaitName+RawFile+BioReplicate+Run+Condition+Proteins,data=data_f_spc,FUN = sum)
  data_f_spc <- subset(data_f_spc, !is.na(spectral_counts))
  
  return( list(data_f = data_f, agg_intensity = data_f_agg, agg_spc = data_f_spc) )
}

my.MaxQToMIST <- function(evidence,key,ref_proteome_file, outdir = "/bf2/smb-data/tnguyen/projects/fluomics/tempdata", quant_variable="spc") {
  # Load and aggregate data/key
  datalist <- my.preprocessAPMS(evidence,key)
  
  quant_variable <- trimws(tolower(quant_variable))
  if (! quant_variable %in% c("spc","intensity")) stop("Please input quant_variable as 'spc' or 'intensity'")
  data_sel <- datalist$agg_spc
  data_sel$ms_spectral_counts <- data_sel$spectral_counts
  quant_col <- 'ms_spectral_counts'
  if (quant_variable == "intensity") {
    data_sel <- datalist$agg_intensity
    data_sel$ms_intensity <- data_sel$Intensity
    quant_col <- 'ms_intensity'
  }
  keysout <- unique(data_sel[,c("RawFile","Condition")])
  
  # Select columns
  data_sel$ms_unique_pep = ""
  data_sel <- data_sel[,c('RawFile','Proteins','ms_unique_pep', quant_col)]
  colnames(data_sel) = c('id','ms_uniprot_ac','ms_unique_pep', quant_col)
  
  # Uniprot annotate
  #conn <- my.bf2_connect()
  #uniprot <- dbGetQuery(conn, "select distinct Entry as ms_uniprot_ac, Length from view_uniprot")
  #dbDisconnect(conn)
  
  # library(dplyr)
  # uniprot <- distinct(uniprot,Entry,Length)
  # colnames(uniprot)[which(colnames(uniprot)=="Entry")]="ms_uniprot_ac"
  # d <- setdiff(data_sel$ms_uniprot_ac, uniprot$ms_uniprot_ac)
  # if (length(d)>0) {
  #   msg <- sprintf("These proteins are not found in uniprot db, please check: %s", paste(d, collapse=","))
  #   #stop(msg)
  #   message(msg)
  # }
  
  ref_proteome <- read.fasta(
    file = ref_proteome_file,
    seqtype = "AA",
    as.string = TRUE,
    set.attributes = TRUE,
    legacy.mode = TRUE,
    seqonly = FALSE,
    strip.desc = FALSE
  )
  p_lengths <- c()
  p_names <- c()
  for (e in ref_proteome) {
    p_lengths <- c(p_lengths, nchar(e[1]))
    p_names <- c(p_names, attr(e, 'name'))
  }
  ref_table <- data.table(names = p_names, lengths = p_lengths)
  ref_table[, uniprot_ac := gsub('([a-z,0-9,A-Z]+\\|{1})([A-Z,0-9,\\_]+)(\\|[A-Z,a-z,0-9,_]+)',
                                 '\\2',
                                 names)]
  ref_table[, uniprot_id := gsub('([a-z,0-9,A-Z]+\\|{1})([a-z,0-9,A-Z]+\\|{1})([A-Z,a-z,0-9,_]+)',
                                 '\\3',
                                 names)]
  colnames(ref_table)[which(colnames(ref_table)=="uniprot_ac")]="ms_uniprot_ac"
  d <- setdiff(data_sel$ms_uniprot_ac, ref_table$ms_uniprot_ac)
  if (length(d)>0) {
    msg <- sprintf("These proteins are not found in fasta, please check: %s", paste(d, collapse=","))
    #stop(msg)
    message(msg)
  }
  
  # Get mass
  data_sel <- base::merge(data_sel, ref_table, by = "ms_uniprot_ac")
  data_sel$Mass <- 110*data_sel$lengths
  
  # Write out
  ## data
  outdir <- paste0(outdir, "/mist/", quant_variable)
  dir.create(outdir, show=FALSE, recursive = TRUE)
  data_outfile <- sprintf("%s/mist-data.%s.txt", outdir, quant_variable)
  write.table(data_sel, file=data_outfile, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
  
  ## keys
  key_outfile <- sprintf("%s/mist-key.%s.txt", outdir, quant_variable)
  write.table(keysout, file=key_outfile, eol='\n', sep='\t', quote=F, row.names=F, col.names=F)
  
  return(list(data_file=data_outfile, keys_file=key_outfile))
}


evidence=read.table(file="/Users/yzhou/Downloads/private.mist/JyotiBatra/MARV02_03_04/evidence.txt",header=T,sep = "\t",stringsAsFactors=F)
keys=read.table(file="/Users/yzhou/Downloads/private.mist/JyotiBatra/MARV02_03_04/keys-2.txt",header = T,sep = "\t",stringsAsFactors=F)
#uniprot <- read.delim(file="/Users/yzhou/Downloads/view_uniprot.txt",header=T,sep = "\t",stringsAsFactors=F)
library(seqinr)
library(data.table)
my.MaxQToMIST(evidence,keys,"/Users/yzhou/Downloads/private.mist/JyotiBatra/MARV02_03_04/human_proteome_and_MARV.fasta", outdir = "/Users/yzhou/Downloads/private.mist/JyotiBatra/MARV02_03_04", quant_variable="spc")

library(yaml)
mist=read_yaml(file="/Users/yzhou/Downloads/private.mist/JyotiBatra/MARV02_03_04/mist/spc/mist.yaml")
mist$files$data="/Users/yzhou/Downloads/private.mist/JyotiBatra/MARV02_03_04/mist/spc/mist-data.spc.txt"
mist$files$keys="/Users/yzhou/Downloads/private.mist/JyotiBatra/MARV02_03_04/mist/spc/mist-key.spc.txt"
mist$files$remove="/Users/yzhou/Downloads/private.mist/JyotiBatra/MARV02_03_04/mist/spc/a.txt"
mist$files$collapse="/Users/yzhou/Downloads/private.mist/JyotiBatra/MARV02_03_04/mist/spc/a.txt"
mist$files$specificity_exclusions="/Users/yzhou/Downloads/private.mist/JyotiBatra/MARV02_03_04/mist/spc/a.txt"
mist$files$output_dir="/Users/yzhou/Downloads/private.mist/JyotiBatra/MARV02_03_04/mist/spc/result"

mist$preprocess$filter_contaminants=0
mist$preprocess$contaminants_file=NULL
mist$preprocess$id_colname="id"
mist$preprocess$prey_colname="ms_uniprot_ac"
mist$preprocess$pepcount_colname="ms_spectral_counts"
mist$preprocess$mw_colname="Mass"

mist$qc$matrix_file=NULL

mist$mist$matrix_file=NULL
mist$mist$weights="fixed"

mist$annotate$enabled=1
mist$annotate$species="HUMAN"
mist$annotate$uniprot_dir="/Users/yzhou/Downloads/uniprot-human-filtered-organism__Homo.tab"

mist$enrichment$enabled=0

write_yaml(mist,file="/Users/yzhou/Downloads/private.mist/JyotiBatra/MARV02_03_04/mist/spc/mist.yaml")



Rscript /Users/yzhou/Downloads/private.mist/main.R -c /Users/yzhou/Downloads/private.mist/JyotiBatra/MARV02_03_04/mist/spc/mist.yaml







evidence_190=read.table(file="/Users/yzhou/Downloads/190.evidence.txt",header=T,sep = "\t",stringsAsFactors=F)
keys_190=read.table(file="/Users/yzhou/Downloads/190.keys.txt",header = T,sep = "\t",stringsAsFactors=F)
uniprot <- read.delim(file="/Users/yzhou/Downloads/view_uniprot.txt",header=T,sep = "\t",stringsAsFactors=F)
my.MaxQToMIST(evidence_190,keys_190, uniprot,outdir = "/Users/yzhou/Downloads/private.mist/JyotiBatra", quant_variable="spc")


