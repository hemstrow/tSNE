##############################
#whole shbang, minus XIX
genotypes <- read.table("../../../Stickleback/Full Data/2017_reruns/snps_numeric_filt.txt", header = T,
                        stringsAsFactors = F, colClasses = "character")

test <- genotypes[1:100,]

pa_genotypes <- format_snps(genotypes[genotypes$group != "groupXIX",], 3, 8, "0000", "00", l_names = genotypes[genotypes$group != "groupXIX",]$snp, interp_miss = T)



library(Rtsne)
library(mmtsne)

hout <- hbeta(pa_genotypes[,-1])

out <- Rtsne(as.matrix(pa_genotypes[-1]), initial_dims = 50, dims = 2, perplexity = hout$H, theta = 0, max_iter = 1000, verbose = TRUE)
dists <- as.data.frame(out$Y)
dists$pop <- substr(pa_genotypes$samp,1,3)

library(ggplot2)
ggplot(dists, aes(x = V1, y = V2, color = pop)) + geom_point() + scale_color_brewer(palette = "Set1")

pca_out <- prcomp(pa_genotypes[,-1])

pca <- as.data.frame(pca_out$x)
pca$pop <- substr(pa_genotypes$samp,1,3)

ggplot(pca, aes(x = PC1, y = PC2, color = pop)) + geom_point() + scale_color_brewer(palette = "Set1")

#############################3
#with XIX only!
paxix_genotypes <- format_snps(genotypes[genotypes$group == "groupXIX",], 3, 8, "0000", "00", l_names = genotypes[genotypes$group == "groupXIX",]$snp, interp_miss = T)



library(Rtsne)
library(mmtsne)

hout <- hbeta(paxix_genotypes[,-1])

out <- Rtsne(as.matrix(paxix_genotypes[-1]), initial_dims = 50, dims = 2, perplexity = hout$H, theta = 0, max_iter = 1000, verbose = TRUE)
dists <- as.data.frame(out$Y)
dists$pop <- substr(pa_genotypes$samp,1,3)

library(ggplot2)
ggplot(dists, aes(x = V1, y = V2, color = pop)) + geom_point() + scale_color_brewer(palette = "Set1")

pca_out <- prcomp(paxix_genotypes[,-1])

pca <- as.data.frame(pca_out$x)
pca$pop <- substr(paxix_genotypes$samp,1,3)

ggplot(pca, aes(x = PC1, y = PC2, color = pop)) + geom_point() + scale_color_brewer(palette = "Set1")



#################################################################################################
#inversion
FST <- read.table("../../../Stickleback/Full Data/2017_reruns/PlotData/FST/LGsnoXIX.txt", header = T)
rFST <- read.table("../../../Stickleback/Full Data/2017_reruns/PlotData/FST/rawfst.txt", header = T)
isnps <- rFST[rFST$group == "groupIX" & rFST$position >= 13*1000000& rFST$position <= 18.5*1000000 & rFST$comp == "ASP_OPL",]
divsnps <- isnps[isnps$Fst >= 0.15, 1]

iv_genotypes <- format_snps(genotypes[genotypes$snp %in% divsnps,], 3, 8, "0000", "00", l_names = divsnps, interp_miss = T)


hout <- hbeta(iv_genotypes[,-1])
dups <- duplicated.matrix(as.matrix(iv_genotypes[,-1]))
iv_genotypes <- iv_genotypes[!dups,]

out <- Rtsne(as.matrix(iv_genotypes[-1]), initial_dims = 50, dims = 2, perplexity = hout$H, theta = 0, max_iter = 1000, verbose = TRUE)
dists <- as.data.frame(out$Y)
dists$pop <- substr(iv_genotypes[!dups,]$samp,1,3)
dists$hap <- haplotypes

library(ggplot2)
ggplot(dists, aes(x = V1, y = V2, color = pop, shape = hap)) + geom_point() + scale_color_brewer(palette = "Set1") +
  theme_bw()

pca_out <- prcomp(t(iv_genotypes[,-1]))

pca <- as.data.frame(pca_out$rotation)
pca$pop <- substr(iv_genotypes$samp,1,3)
pca$hap <- haplotypes

ggplot(pca, aes(x = PC1, y = PC2, color = pop, shape = hap)) + geom_point() + scale_color_brewer(palette = "Set1") +
  theme_bw()


