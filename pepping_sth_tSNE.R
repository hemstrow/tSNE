setwd("../tSNE_data/steelhead/snps/")
library(snpR)
###########################
#import and prepare data
# sth.geno <- read.table("pepping_sth_snps.geno", sep = "\t", header = F, stringsAsFactors = F)
# sth.geno <- sth.geno[,-ncol(sth.geno)]
# 
# #filter for non poly, non bi, hf hets only
# flt.geno <- filter_snps(sth.geno, 2, non_poly = T, hf_hets = 0.55)
# 
# pa.geno <- format_snps(flt.geno, 2, 7, lnames = paste0(flt.geno[,1], "_", flt.geno[,2]))
# 
# #add parentage info
# fams <- read.table("pepping_sth_families.txt", header = T)
# fams <- fams[fams$Run != "L",]
# pa.geno <- cbind(fams, pa.geno)
# 
# saveRDS(pa.geno, "pa.geno.RDS")
pa.geno <- readRDS("pa.geno.RDS")

##########################
#run tSNE and PCA

pca <- PCAfromPA(pa.geno, 6, "Family")
pca$plot <- pca$plot + ggplot2::guides(color = FALSE)
pca$plot

tSNE <- tSNEfromPA(pa.geno, 6, "Family")
tSNE$plot <- tSNE$plot + ggplot2::guides(color = FALSE)
tSNE$plot

#interesting, great with some crappy with others... try coloring points by inclusive family prob?