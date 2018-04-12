library(snpR); library(ggplot2)
setwd("../tSNE_data/steelhead/snps/RAPTURE/")
#==========================Import metadata====================

bamlist <- read.table("bamlist.txt")
bamlist <- as.vector(t(bamlist))

#get ind IDs
samps <- gsub(".+bam/", "", bamlist)
samps <- gsub("_.+$", "", samps)

#get pops
pops <- gsub(".+[0-9]{4}_", "", bamlist)
pops <- gsub("_.+", "", pops)

#get families
fams <- gsub(".+_", "", bamlist)
fams <- gsub("\\.bam", "", fams)

#===========================IBS====================
setwd("new/")

##########
#IBS mat, doesn't work very well. Ask Mike?
# IBS <- "rapIBS_clean.ibsMat"
# IBS <- read.table(IBS)
# dIBS <- as.dist(IBS)
# #set colnames
# colnames(IBS) <- (samps)
# 
# #append
# IBS <- cbind(data.frame(pop = pops, fam = fams, samp = samps), IBS)
# 
# #remove nas
# minmiss <- rowSums(ifelse(as.matrix(is.na(IBS)), 1, 0))
# IBS <- IBS[-which(minmiss > min(minmiss)),]
# IBS <- IBS[,-(which(minmiss > min(minmiss)) + 3)]
# 
# # mds <- cmdscale(as.dist(as.matrix(IBS[,-c(1:3)])))
# # mds <- as.data.frame(mds)
# # mds <- cbind(IBS[,1:3], mds)
# 
# ggplot(mds, aes(V1, V2, color = fam)) + geom_point() + guides(color = FALSE)
# 
# pca <- PCAfromPA(IBS, 3, "fam")
# pca$plot + guides(color = FALSE)
# 
# tsne <- tSNEfromPA(IBS, 3, "fam")
# tsne$plot + guides(color = FALSE)

#############
#IBS raw after filtering, works well. Hard to say if better than genos or not!
IBS <- read.table("rapIBS_clean.ibs", header = T)


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
# 
# #remove poorly sequenced loci
# cuttoff <- 0.5
# IBS[IBS == -1] <- NA
# IBS <- IBS[-which(rowSums(is.na(IBS)) >= (ncol(IBS)- 4)*(1-cuttoff)),]
# 
# #remove poorly sequenced individuals
# icut <- 0.5
# IBS <- t(IBS[,-c(1:4)])
# IBS <- cbind(data.frame(pop = pops, fam = fams, samp = samps), IBS)
# IBS <- IBS[-which(rowSums(is.na(IBS)) >= (ncol(IBS) - 3)*(1-icut)),]
# 
# #refilter poor loci
# cuttoff2 <- 0.75
# IBS <- IBS[,-which(colSums(is.na(IBS)) >= nrow(IBS)*(1-cuttoff2))]
# 
# miss <- rowSums(is.na(IBS))
# 
# #interpolate missing--see snpR::format_snps for the source of this code
# dat <- as.matrix(IBS[,4:ncol(IBS)])
# afs <- colMeans(dat, TRUE)
# temp <- which(is.na(dat))/nrow(dat)
# fill_element <- floor(temp) + 1 #get the column for each missing data point
# fill_element[which(temp %% 1 == 0)] <- fill_element[which(temp %% 1 == 0)] - 1 #correct for anything in the last row
# dat[which(is.na(dat))] <- afs[fill_element] #fill with the appropriate allele frequency.
# 
# IBS <- cbind(IBS[,1:3], dat)

IBS <- prepIBS(IBS[,-c(1:4)], .5, .75, .5, data.frame(samp = samps, pop = pops, fam = fams))

#get the pca and tSNE
pca <- PCAfromPA(IBS, 4, "fam")
pca$plot + guides(color = FALSE)

tSNE <- tSNEfromPA(IBS, 4, "fam")
tSNE$plot + guides(color = FALSE)
dat <- tSNE$plot$data

ggplot(dat, aes(x = V1, V2, color = miss)) + geom_point()


#==============================Posts===============
Posts <- read.table("rapPosts_clean.geno")
Posts <- Posts[,-c(1:2)]
meta <- data.frame(samp = samps, pop = pops, fam = fams)

#remove crappy samples

#get the expected number of each allele
PG1 <- t(Posts[,seq(1,ncol(Posts), 3)])
PG2 <- t(Posts[,seq(2,ncol(Posts), 3)])
PG3 <- t(Posts[,seq(3,ncol(Posts), 3)])

#figure out likely missing data
missing <- ifelse(round(PG1, 3) == 0.333 & round(PG2, 3) == 0.333 & round(PG3, 3) == 0.333, 1, 0)

#get the expected number of each allele
PMaj <- PG1*2 + PG2
PMin <- PG3*2 + PG2
remove(PG2)

#remove crappy alleles
cuttoff <- 0.5
cl <- which(colSums(missing)/nrow(missing) > (1-cuttoff))
if(length(cl) > 0){
  PMaj <- PMaj[,-cl]
  PMin <- PMin[,-cl]
  missing <- missing[,-cl]
}

#remove crappy samples
cuttoff <- 0.5
cl <- which(rowSums(missing)/ncol(missing) > (1-cuttoff))
if(length(cl) > 0){
  PMaj <- PMaj[-cl,]
  PMin <- PMin[-cl,]
  meta <- meta[-cl,]
  missing <- missing[-cl,]
}

#refilter loci
cuttoff <- 0.75
cl <- which(colSums(missing)/nrow(missing) > (1-cuttoff))
if(length(cl) > 0){
  PMaj <- PMaj[,-cl]
  PMin <- PMin[,-cl]
  missing <- missing[,-cl]
}

PGs <- cbind(meta, miss = rowSums(missing), as.data.frame(cbind(PMaj, PMin)))

tSNEp <- tSNEfromPA(PGs, 4, "fam", initial_dims = floor(length(unique(PGs$fam))*2/3), perplex = 30)
tSNEp$plot + guides(color = FALSE)

dat <- tSNEp$plot$data

ggplot(dat, aes(V1, V2, color = miss)) + geom_point(size = 6)
ggplot(dat, aes(V1, V2, color = fam)) + geom_point(size = 6) + guides(color = FALSE)
