prepIBS <- function(IBS, lcut1, lcut2, icut, meta, missing = -1){
  #filter poor loci
  IBS[IBS == missing] <- NA
  IBS <- IBS[-which(rowSums(is.na(IBS)) >= (ncol(IBS))*(1-lcut1)),]
  
  #remove poorly sequenced individuals
  IBS <- t(IBS)
  IBS <- cbind(meta, IBS)
  IBS <- IBS[-which(rowSums(is.na(IBS)) >= (ncol(IBS) - ncol(meta))*(1-icut)),]
  
  #refilter poor loci
  IBS <- IBS[,-which(colSums(is.na(IBS)) >= nrow(IBS)*(1-lcut2))]
  
  miss <- rowSums(is.na(IBS))
  
  #interpolate missing--see snpR::format_snps for the source of this code
  dat <- as.matrix(IBS[,-c(1:ncol(meta))])
  afs <- colMeans(dat, TRUE)
  temp <- which(is.na(dat))/nrow(dat)
  fill_element <- floor(temp) + 1 #get the column for each missing data point
  fill_element[which(temp %% 1 == 0)] <- fill_element[which(temp %% 1 == 0)] - 1 #correct for anything in the last row
  dat[which(is.na(dat))] <- afs[fill_element] #fill with the appropriate allele frequency.
  
  IBS <- cbind(as.data.frame(cbind(IBS[,1:3], miss = miss)), dat)
  return(IBS)
}