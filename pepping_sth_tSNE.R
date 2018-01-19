###########################
#import and prepare data
sth.geno <- read.table("pepping_sth_snps.geno", sep = "\t", header = F, stringsAsFactors = F)
sth.geno <- sth.geno[,-ncol(sth.geno)]

#filter for non poly, non bi, hf hets only
flt.geno <- filter_snps(sth.geno, 2, T, T, F, F, 0.6)

pa.geno <- format_snps(flt.geno, 2, 8, l_names = paste0(flt.geno[,1], "_", flt.geno[,2]))

##########################
#run tSNE and PCA

library(Rtsne)
library(mmtsne)

hout <- hbeta(pa.geno[,-1])

out2 <- Rtsne(as.matrix(pa.geno[-1]), initial_dims = 50, dims = 2, 
             perplexity = hout$H, theta = 0, max_iter = 2000, verbose = TRUE,
             check_duplicates = FALSE)

dists <- as.data.frame(out2$Y)

library(ggplot2)
ggplot(dists, aes(x = V1, y = V2)) + geom_point()

#what else floats on water?

pca_out <- prcomp(pa.geno[,-1])

pca <- as.data.frame(pca_out$x)

ggplot(pca, aes(x = PC1, y = PC2)) + geom_point()



