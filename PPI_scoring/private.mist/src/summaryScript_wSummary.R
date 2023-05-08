#! /usr/bin/Rscript
# inputs: -d data_file, -o output_file, -k keys_file: Use MT_keys
####################################################
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))

# get the summarie values of each variable
getSums <- function(y,baits){
	temp <- c(bait_name=unique(y[,2]), bait_uniprot_ac=unique(y[,1]), ID=unique(y$ID), num_prey=length(unique(y[,4])))
  if(length(unique(y[,2]))>1){
    stop(paste("MULTIPLE BAITS LISTED FOR A SINGLE ID: ", unique(y$ID), sep="" ) )
  }
	temp <- c(temp, num_ranks=length(y[,4][-grep("-", y[,4])]))
	
	top_ranked_protein=y[which(y$ms_rank=="[1]"),]$ms_uniprot_ac
	idx <- which(y[,6]==temp[2]) #location of bait in the prey list
	if(length(idx)>0){
	  if(length(idx)>1){   # if the bait has the same rank in more than one sample
      temp1 = y[idx,c(4,7,8)]
      names(temp1) = c('bait_rank','pep_count','cov')
      temp2 = as.data.frame(matrix(rep(temp, length(idx)), nrow=length(idx), byrow=T), stringsAsFactors=F)
      names(temp2) = names(temp)
      temp <- cbind(temp2, temp1, top_ranked_protein=top_ranked_protein)
	  }else{
	    temp <- c(temp, bait_rank=y[idx,4], pep_count=y[idx,7], cov=y[idx,8], top_ranked_protein=top_ranked_protein)
	  }
	}else{
		temp <- c(temp, bait_rank="0", pep_count=0, cov=0, top_ranked_protein=top_ranked_protein)
	}
	
  temp = as.data.frame(t(temp), stringsAsFactors=F)
  #temp <- cbind(as.data.frame(t(temp), stringsAsFactors=F), top_ranked_protein)  # allow for multiple equally ranked top ranks
  
  # find highest ranking bait
  idx = which(y$ms_uniprot_ac %in% baits)
  if(length(idx)>0){
    idx = idx[sortRank(y[idx,'ms_rank'])][1]
    temp <- cbind(temp, top_bait=y[idx,'ms_uniprot_ac'], top_bait_pep_count=y[idx,7], top_bait_cov=y[idx,'ms_x_cov'], stringsAsFactors=F)
  }else{
    temp <- cbind(temp, top_bait=NA, top_bait_pep_count=0, top_bait_cov=0, stringsAsFactors=F)
  }
	return(temp)
}

# return index of sorted ranks (in vector form)
sortRank <- function(x){
  # remove brackets
  x = gsub("\\[|\\]","",x)
  x = as.numeric(gsub("-",".", x))
  idx = order(x)
  return(idx)
}

# finds which proteins are unique for that experiment
findProteins <- function(y, temp){
	prots <- which(table(y[,'ms_uniprot_ac'])<2)
	idx <- match(names(prots), y[,'ms_uniprot_ac'])
	tmp <- cbind(names(prots), y$ID[idx])
	tmp <- table(tmp[,2])
	tmp <- as.data.frame(cbind(names(tmp),tmp), stringsAsFactors=FALSE)
	names(tmp) <- c("ID", "original_peptides")
	tmp <- merge(temp, tmp, by="ID", all.x=TRUE)
  # Re-order the columns
  tmp = tmp[,c('ID', 'bait_name', 'bait_uniprot_ac', 'num_prey', 'num_ranks', 'original_peptides', 'bait_rank', 'pep_count', 'cov', 'top_ranked_protein', 'top_bait', 'top_bait_pep_count', 'top_bait_cov')]
	return(tmp)
}


# wrapper to summarize each bait ID
summarize <- function(x){
	tmp <- c()
  baits = as.character(unique(x$bait_uniprot_ac))
	for(i in unique(x$ID)){
		idx <- which(x$ID == i)
		tmp <- rbind(tmp, getSums(x[idx,], baits))
	}
	tmp <- findProteins(x, tmp)
	# remove NA's
	tmp$original_peptides <- as.character(tmp$original_peptides)
	tmp$original_peptides[is.na(tmp$original_peptides)] <- 0
	tmp <- as.data.frame(tmp, stringsAsFactors=FALSE)
	tmp <- tmp[order(tmp$bait_rank),]
  
	return(tmp)
}

# get bait names from summary file and merge them with the Prospector data file
mergeWkeys <- function(datafile, summaryfile){
	dat <- read.delim(datafile, header=TRUE, stringsAsFactors=FALSE, sep="\t")
	keys <- read.delim(summaryfile, header=TRUE, stringsAsFactors=FALSE, sep="\t")
 
  # rename the columns named id
  names(dat)[grep('^id$|^ID$', names(dat))] = names(keys)[grep('^id$|^ID$', names(keys))] = 'ID'
	# use consistent column name for accession numbers
  if(!any(grep('^bait_uniprot_ac$', names(keys))))  #check if 'bait_uniprot_ac' already exists
    names(keys)[grep('^bait_uniprot_id$', names(keys))] = 'bait_uniprot_ac'
	names(dat)[grep('^Acc\\.\\.$', names(dat))] = 'ms_uniprot_ac'
	names(dat)[grep('^Rank$', names(dat))] = 'ms_rank'
	names(dat)[grep('^Num\\.Unique$', names(dat))] = 'ms_num_unique_peptide'
	names(dat)[grep('cov|Cov', names(dat))] = 'ms_x_cov'
	names(dat)[grep('^Peptide\\.Count$', names(dat))] = 'top_bait_pep_count'
	names(dat)[grep('Spectral|spectral|Count|count', names(dat))] = 'ms_spectral_counts'
	
  
  dat$ms_uniprot_ac = toupper(dat$ms_uniprot_ac)
  keys$bait_uniprot_ac = toupper(keys$bait_uniprot_ac)
	
  dat <- merge(dat, keys[,c('ID','bait_name', 'bait_uniprot_ac')], by="ID", all.x=TRUE)
	n <- dim(dat)[2]
	x1 <- dat[,c(n,n-1)]
	x2 <- dat[,-c(n,n-1)]
	dat <- cbind(x1, x2)
	return(unique(dat))
}

# list all occurances of a bait protein being pulled down
listBaits <- function(dat, outdir){
  baits <- unique(dat[,'bait_uniprot_ac'])
  baits = baits[!is.na(baits)]  # remove NA's
  idx <- which(!is.na(match(dat[,'ms_uniprot_ac'], baits)))
  if(length(idx) == 0){
    cat("  !!  No Baits found in the Preys!!")
    break
  }
  listOfBaits <- dat[idx,c(3,1,2,4:dim(dat)[2])]
  # order them by experiment
  idx <- orderExperiments(listOfBaits)
  listOfBaits <- listOfBaits[idx,]
  
  # get bait-bait pairs
  temp <- unique(listOfBaits[listOfBaits$bait_uniprot_ac == listOfBaits$ms_uniprot_ac, c("bait_uniprot_ac","bait_name")])
  names(temp)[which("bait_name" == names(temp))] = "ms_uniprot_name"
  listOfBaits = merge(listOfBaits, temp, by.x='ms_uniprot_ac', by.y='bait_uniprot_ac', all.x=T)
  n = dim(listOfBaits)[2]
  # reorder columns
  listOfBaits = listOfBaits[,c(2:4,1,n,5:(n-1))]  
  # reorder by rank & id
  listOfBaits = listOfBaits[sortRank(listOfBaits$ms_rank),]
  listOfBaits = listOfBaits[orderExperiments(listOfBaits),]
  
  write.table(listOfBaits, paste(outdir,"baits_found.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  return(listOfBaits)
}

# plot heatmapt
carryoverPlot <- function(listOfBaits, outfile, sc){
  x = listOfBaits
  x$pairs = paste(x$ID, x$bait_name, sep="-")
  x$pairs <- factor(x$pairs, levels = unique(x$pairs))
  
  x$ms_rank = gsub('\\[|\\]|-.*', '', x$ms_rank)
  
  #breaks = unique(c( seq(min(x$ms_spectral_counts),quantile(x$ms_spectral_counts,.25),length.out=2), seq(quantile(x$ms_spectral_counts,.26),quantile(x$ms_spectral_counts,.50),length.out=2), seq(quantile(x$ms_spectral_counts,.51),quantile(x$ms_spectral_counts,.75),length.out=2), seq(quantile(x$ms_spectral_counts,.76),quantile(x$ms_spectral_counts,1),length.out=2) ))
  
  
  height = max(7, length(unique(x$ms_uniprot_name))/3.9)  #  78/20
  width = max(7, length(unique(x$pairs))/5)               # 273/55
  #x$ms_num_unique_peptide = log(x$ms_num_unique_peptide)
  
  
  if(sc==1){
    n = x$ms_spectral_counts
    ggplot(x, aes(pairs, ms_uniprot_name)) + geom_tile(aes(fill = ms_spectral_counts), colour = "white") + 
      scale_fill_gradient(limits=c(min(n),as.vector(quantile(n,.95))), low = "steelblue", high = "red", breaks=c( as.vector(quantile(n, c(.25,.5,.75,.95)))), name="Peptide Count" ) + 
      geom_text(aes(fill = ms_spectral_counts, label = ms_spectral_counts), size = 3) + 
      labs(title="Baits in Experiments as they appear.") + xlab("ID - bait name") + ylab("Bait uniprot ac") + theme(axis.text.x = element_text(angle=45, hjust=1))
  }else{    
    n = x$ms_num_unique_peptide
    ggplot(x, aes(pairs, ms_uniprot_name)) + geom_tile(aes(fill = ms_num_unique_peptide), colour = "white") + 
      scale_fill_gradient(limits=c(min(n),as.vector(quantile(n,.95))), low = "steelblue", high = "red", breaks=c( as.vector(quantile(n, c(.25,.5,.75,.95)))), name="Peptide Count" ) + 
      geom_text(aes(fill = ms_num_unique_peptide, label = ms_num_unique_peptide), size = 3) + 
      labs(title="Baits in Experiments as they appear. (Numbers = rank, colors = peptide count)") + xlab("ID - bait name") + ylab("Bait uniprot ac") + theme(axis.text.x = element_text(angle=45, hjust=1))
  }
  ggsave(outfile, height=height, width=(width*1.75), limitsize=FALSE)
  
}

# order the data based on it's id
orderExperiments <- function(y){
  idnum <- sub("(^[A-Za-z]+|^[A-Za-z]+)-","",y$ID)	#strip out numbers
  idname <- sub("([0-9].*|-[0-9].*)","",y$ID)			#strip out id characters
  tmp = strsplit(idnum,"-")
  tmp = matrix( unlist(lapply(tmp, function(x) if(length(x)>1){c(x[1], x[2])}else{c(x[1],0)} )), ncol=2, byrow=T)
  tmp = cbind(idname, tmp)
  tmp[,2] = sub("[^0-9]+","",tmp[,2])  #added to remove extra characters/puctuation here after number like "rerun"
  idx = order(tmp[,1], as.numeric(tmp[,2]), as.numeric(tmp[,3]))
  return(idx)
}

# begin summarizing
main <- function(datafile, summaryfile, outdir, sc){
	dat <- mergeWkeys(datafile, summaryfile)
	x <- summarize(dat)
  x <- x[orderExperiments(x),]
  outfile = paste(outdir,"data_summarized.txt",sep="")
  write.table(x, outfile , sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  # save matrix of all bait proteins found in all samples
  x <- listBaits(dat, outdir)
 
  # plot a heatmap of the frequenct of baits
	x <- x[orderExperiments(x),]
  outfile = gsub('data_summarized.txt','carryover_plot_ordered_by_ID.pdf',outfile)
	carryoverPlot(x, outfile, sc)
  
	outfile = gsub('carryover_plot_ordered_by_ID', 'carryover_plot_ordered_by_BAIT',outfile)
  x <- x[order(x$bait_name),] 
	carryoverPlot(x, outfile, sc)
  
#   if(sc==1){
#     x1 = dcast(x, ms_uniprot_ac~ID+bait_uniprot_ac, value.var="ms_spectral_counts",fill=0)
#     write.table(x1, paste(outdir,"baits_matrix_spectral_counts.txt",sep=""), , sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)    
#   }else{
#     x1 = dcast(x, ms_uniprot_ac~ID+bait_uniprot_ac, value.var="ms_num_unique_peptide",fill=0)
#     write.table(x1, paste(outdir,"baits_matrix_ms_num_peptide.txt",sep=""), , sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
#   }
# 	x1 = dcast(x, ms_uniprot_ac~ID+bait_uniprot_ac, value.var="ms_rank",fill=0)
# 	write.table(x1, paste(outdir,"baits_matrix_ms_rank.txt",sep=""), , sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
	
	return(x)
}


option_list <- list(
  make_option(c("-d", "--data_file"), 
              help="data file containing Prospector output"),
  make_option(c("-o", "--output_dir"),
              help="output file for summary matrix"),
  make_option(c("-s", "--spectral_counts", default=0), 
              help="Use spectral counts instead of num_unique_pep"),
  make_option(c("-k", "--summary_file"), 
              help="summary file (txt)") #should be from MT          
)

parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
summarizedData <- main(parsedArgs$data_file, parsedArgs$summary_file, parsedArgs$output_dir, parsedArgs$spectral_counts)

#write.table(summarizedData, parsedArgs$output_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


###################################
# datafile= "~/Box Sync/Projects/Dengue/20140904 Mosquito Rerun/20140904 PS-327-341-proteins.txt"
# summaryfile = "~/projects/mist/tests/krogan/dengue/summary/summary.txt"
# #keysfile= "~/Desktop/MT_keys_merged.tsv"
# outdir= "~/projects/mist/tests/krogan/dengue/summary/"
# sc=1
# x= main(datafile, summaryfile, outdir, sc)


# datafile= "~/projects/mist/tests/engel/elwell/input/elwell_data.txt"
# summaryfile = "~/projects/mist/tests/engel/elwell/input/elwell_summary.txt"
# #keysfile= "~/Desktop/MT_keys_merged.tsv"
# outdir= "~/projects/mist/tests/engel/elwell/summary/"
# spectral_counts=sc=0
# x= main(datafile, summaryfile, outdir, sc)



# datafile= "~/projects/mist/tests/NV/input/data.txt"
# summaryfile = "~/projects/mist/tests/NV/input/summary.txt"
# #keysfile= "~/Desktop/MT_keys_merged.tsv"
# outdir= "~/projects/mist/tests/NV/summary/"
# spectral_counts=sc=0
# x= main(datafile, summaryfile, outdir, sc)

# datafile= "~/Box Sync/Projects/WNV/data/base/WNV_data.txt"
# summaryfile = "~/Box Sync/Projects/WNV/summary.txt"
# #keysfile= "~/Desktop/MT_keys_merged.tsv"
# outdir= "~/Desktop/"
# spectral_counts=sc=0
# x= main(datafile, summaryfile, outdir, sc)





