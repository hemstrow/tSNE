library(ggplot2)
library(snpR)
########################################################################################################
#set up data
setwd("~/GitHub/tSNE_data/monarch/")
#genotypes <- read.table("pa_snps.rds", header = T, stringsAsFactors = F)

#pa_snps <- format_snps(genotypes, 3, 7)
# pops <- pa_snps$samp
# pops <- as.character(pops)
# pops <- substr(pops, 1, 3)
# pa_snps$Population <- pops
# pa_snps <- pa_snps[,c(1, ncol(pa_snps), 2:(ncol(pa_snps) - 1))]

#saveRDS(pa_snps, "pa_snps.rds")

#read in formated data
pa_snps <- readRDS("pa_snps.rds")
colnames(pa_snps)[2] <- "Population"
#################################################################
#run the PCA and tSNE
#PCA:
pca <- PCAfromPA(pa_snps, 2, do.plot = "Population")
pca$plot

#tSNE
tSNE <- tSNEfromPA(pa_snps, 2, do.plot = "Population")
tSNE$plot



