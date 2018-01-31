library(snpR)
##############################
#get data
setwd("~/GitHub/tSNE_data/stickleback")


#format
genos <- read.table("snps_numeric_filt.txt", header = T, stringsAsFactors = F, colClasses = "character")


##full
#pa_genos <- format_snps(genos, 3, 7, "0000", "00")
#pa_genos$pop <- gsub("\\.\\w+", "", pa_genos$samp)
#pa_genos <- pa_genos[,c(1, ncol(pa_genos), 2:ncol(pa_genos))]
#pa_genos <- pa_genos[,-ncol(pa_genos)]
#saveRDS(pa_genos, "pa_genos.RDS")

##noXIX
#pa_nXIX <- format_snps(genos[genos$group != "groupXIX",], 3, 7, "0000", "00")
#pa_nXIX$pop <- gsub("\\.\\w+", "", pa_nXIX$samp)
#pa_nXIX <- pa_nXIX[,c(1, ncol(pa_nXIX), 2:ncol(pa_nXIX))]
#pa_nXIX <- pa_nXIX[,-ncol(pa_nXIX)]
#saveRDS(pa_nXIX, "pa_nXIX.RDS")

##XIX only
#pa_XIX <- format_snps(genos[genos$group == "groupXIX",], 3, 7, "0000", "00")
#pa_XIX$pop <- gsub("\\.\\w+", "", pa_XIX$samp)
#pa_XIX <- pa_XIX[,c(1, ncol(pa_XIX), 2:ncol(pa_XIX))]
#pa_XIX <- pa_XIX[,-ncol(pa_XIX)]
#saveRDS(pa_XIX, "pa_XIX.RDS")

##inversion only
#rFST <- read.table("~/Stickleback/Full Data/2017_reruns/PlotData/FST/rawfst.txt", header = T)
#isnps <- rFST[rFST$group == "groupIX" & rFST$position >= 13*1000000& rFST$position <= 18.5*1000000 & rFST$comp == "ASP_OPL",]
#isnps <- isnps[isnps$Fst >= 0.15, 1]
#pa_IX <- format_snps(genos[genos$snp %in% isnps,], 3, 7, "0000", "00")
#pa_IX$pop <- gsub("\\.\\w+", "", pa_IX$samp)
#pa_IX <- pa_IX[,c(1, ncol(pa_IX), 2:ncol(pa_IX))]
#pa_IX <- pa_IX[,-ncol(pa_IX)]
#saveRDS(pa_IX, "pa_IX.RDS")

#read in
pa_genos <- readRDS("pa_genos.RDS")
pa_XIX <- readRDS("pa_XIX.RDS")
pa_nXIX <- readRDS("pa_nXIX.RDS")
pa_IX <- readRDS("pa_IX.RDS")


###############################
#run PCA and tSNE on the full data set

fpca <- PCAfromPA(pa_genos, 2)
ftSNE <- tSNEfromPA(pa_genos, 2)

###############################
#with XIX only!

XIXpca <- PCAfromPA(pa_XIX, 2)
XIXtSNE <- tSNEfromPA(pa_XIX, 2)

###############################
#with noXIX

nXIXpca <- PCAfromPA(pa_nXIX, 2)
nXIXtSNE <- tSNEfromPA(pa_nXIX, 2)

###############################
#inversion
IXpca <- PCAfromPA(pa_IX, 2)
IXtSNE <- tSNEfromPA(pa_IX, 2)

