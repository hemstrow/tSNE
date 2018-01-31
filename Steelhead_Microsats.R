library(ggplot2)
###########################################
#get data and format.
setwd("~/GitHub/tSNE_data/steelhead/msats/")

#genos <- read.table("steelhead_genotypes.txt", sep='', header=FALSE, colClasses = "character")
#genos <- t(genos)
#pops <- genos[1,]
#genos <- genos[-1,]
#colnames(genos) <- pops
#lnam <- c("Oki23", "Ssa407", "mSsa408", "Ots209", "OtsG249b", "OtsG85", "Omy27", "Omy1001", "Ots243", "Ots409", "OtsG3", "Ots212", "Omm1087")
#genos <- cbind(locus = lnam, lnum = 1:13, genos)

#pa_genos <- format_snps(genos, 2, 7, "msat_2", "00", lnames = lnam)
#pa_genos <- cbind(samp = 1:nrow(pa_genos), pop = pa_genos$samp, pa_genos[,-1])

#saveRDS(pa_genos, "pa_genos.RDS")

#read data in
pa_genos <- readRDS("pa_genos.RDS")

#######################################
#run PCA and tSNE

pca <- PCAfromPA(pa_genos, 2, TRUE)
tSNE <- tSNEfromPA(pa_genos, 2, c.dup = TRUE)
