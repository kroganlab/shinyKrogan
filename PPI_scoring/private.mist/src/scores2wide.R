suppressMessages(library(dplyr))


# takes a Mist scores (results) file and converts it to "wide" format.
scores2wide <-function(dat, outfile){
  # dat <- unique(read.delim(scores_file, stringsAsFactors = F))
  # save protein annotations for later
  dat.labels <- unique(dat[,c('Prey', 'Entry.name','Protein.names', 'Gene.names')])
  
  # get percent ranking of interactions globally based on MIST scores
  dat <- dat %>% mutate(mist_percent_ranking=percent_rank(MIST)) %>% mutate(mist_ranking=rank(desc(MIST), ties.method = 'min'))
  # get percent ranking of interactions within each bait based on MIST scores
  dat <- dat %>% group_by(Bait) %>% mutate(mist_percent_ranking_perBait=percent_rank(MIST)) %>% mutate(mist_ranking_perBait=rank(desc(MIST), ties.method = 'min'))
  
  # remove the "group by" 
  dat <- as.data.frame(dat, stringsAsFactors=F)
  
  # get percent ranking of interactions globally based on COMPPASS scores
  dat <- dat %>% mutate(wd_percent_ranking=percent_rank(WD)) %>% mutate(WD_ranking=rank(desc(WD), ties.method = 'min'))
  # get percent ranking of interactions within each bait based on MIST scores
  dat <- dat %>% group_by(Bait) %>% mutate(wd_percent_ranking_perBait=percent_rank(WD)) %>% mutate(wd_ranking_perBait=rank(desc(WD), ties.method = 'min'))
  
  head(dat[,c("Bait","mist_percent_ranking", "mist_ranking", "mist_percent_ranking_perBait","mist_ranking_perBait")])
  head(dat[,c("Bait","wd_percent_ranking", "WD_ranking", "wd_percent_ranking_perBait","wd_ranking_perBait")])
  
  
  # write Long version results out. Includes the percentages
  outfile <- gsub(".txt", "-long_percentiles.txt", outfile)
  write.table(dat, outfile, quote=F, sep='\t', row.names=F)
  
  # create wide format of the data
  # Reproducibility wide format
  R <- dcast(dat, Prey~Bait, value.var='Reproducibility')
  names(R)[-1] <- paste0(names(R)[-1], "_Reproducibility")
  # Abundance wide format
  A <- dcast(dat, Prey~Bait, value.var='Abundance')
  names(A)[-1] <- paste0(names(A)[-1], "_Abundance")
  # Specificity wide format
  S <- dcast(dat, Prey~Bait, value.var='Specificity')
  names(S)[-1] <- paste0(names(S)[-1], "_Specificty")
  
  # mist wide format
  mist <- dcast(dat, Prey~Bait, value.var='MIST')
  names(mist)[-1] <- paste0(names(mist)[-1], "_mist")
  # comppass wide format
  wd <- dcast(dat, Prey~Bait, value.var='WD')
  names(wd)[-1] <- paste0(names(wd)[-1], "_wd")
  # mist as a global percentage wide format
  mist_percent_global <- dcast(dat, Prey~Bait, value.var="mist_percent_ranking")
  names(mist_percent_global)[-1] <- paste0(names(mist_percent_global)[-1], "_mist_percent_global")
  # wd as a global percentage wide format
  wd_percent_global <- dcast(dat, Prey~Bait, value.var="wd_percent_ranking")
  names(wd_percent_global)[-1] <- paste0(names(wd_percent_global)[-1], "_wd_percent_global")
  # wd as a percentage per bait wide format
  wd_percent_perBait <- dcast(dat, Prey~Bait, value.var="wd_percent_ranking_perBait")
  names(wd_percent_perBait)[-1] <- paste0(names(wd_percent_perBait)[-1], "_wd_percent_perBait")
  
  # combine all the percentages into a single file
  dat.wide <- data.frame(stringsAsFactors = F, R, A[,-1], S[,-1], mist[,-1], wd[,-1], mist_percent_global[,-1], wd_percent_global[,-1], wd_percent_perBait[,-1])
  
  # Annotate the data
  dat.wide <- merge(dat.wide, dat.labels, by='Prey')
  # write results out
  outfile <- gsub(".txt", "-wide.txt", outfile)
  write.table(dat.wide, outfile, quote=F, sep='\t', row.names=F)
}



