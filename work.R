genotypes <- read.table("../../../Stickleback/Full Data/2017_reruns/snps_numeric_filt.txt", header = T,
                        stringsAsFactors = F, colClasses = "character")


sim_snps <- format_snps(genotypes, 3, 7, "0000", "00")








########################################
#all except group XIX
tsne_snps <- sim_snps[sim_snps$group != "groupXIX",]
tsne_snps <- tsne_snps[,-c(1:3)]
tsne_snps <- t(tsne_snps)

#fill missing with locus averages
lavs <- colMeans(tsne_snps, na.rm = T)
test <- tsne_snps
for (i in 1:ncol(test)){
  test[is.na(test[,i]),i] <- lavs[i]
}

#tsna
library(Rtsne)
library(mmtsne)
hbout <- hbeta(test)
out <- Rtsne(test, initial_dims = 50, dims = 2, perplexity = hbout$H, theta = 0, max_iter = 1000, verbose = TRUE)

dists <- as.data.frame(out$Y)
dists$pop <- substr(colnames(genotypes[,4:ncol(genotypes)]),1,3)
ggplot(dists, aes(x = V1, y = V2, color = pop)) + geom_point() + scale_color_brewer(palette = "Set1")

#pca
pca_out <- prcomp(t(test))
summary(pca_out)
pcdat <- as.data.frame(pca_out$rotation)

pcdat <- cbind(pop = substr(colnames(genotypes[,4:ncol(genotypes)]),1,3), pcdat)

ggplot(pcdat, aes(x = PC1, y = PC2, color = pop)) + geom_point() + scale_color_brewer(palette = "Set1")
ggplot(pcdat, aes(x = PC2, y = PC3, color = pop)) + geom_point() + scale_color_brewer(palette = "Set1")









##########################################################################################
#divergent inversion snps only
#first, get list of high FST group IX snps
FST <- read.table("../../../Stickleback/Full Data/2017_reruns/PlotData/FST/LGsnoXIX.txt", header = T)
rFST <- read.table("../../../Stickleback/Full Data/2017_reruns/PlotData/FST/rawfst.txt", header = T)
ggplot(rFST[rFST$group == "groupIX" & rFST$comp == "ASP_OPL",], 
       aes(x = position, y = Fst)) + geom_point()
isnps <- rFST[rFST$group == "groupIX" & rFST$position >= 13*1000000& rFST$position <= 18.5*1000000 & rFST$comp == "ASP_OPL",]
divsnps <- isnps[isnps$Fst >= 0.15, 1]
divpos <- isnps[isnps$Fst >= 0.15, 2]

div_tsne <- sim_snps[sim_snps$snp %in% divsnps, 4:ncol(sim_snps)]
div_tsne <- t(div_tsne)


#fill missing with locus averages
lavs <- colMeans(div_tsne, na.rm = T)
test2 <- div_tsne
for (i in 1:ncol(test)){
  test2[is.na(test2[,i]),i] <- lavs[i]
}

#tsna
library(Rtsne)
library(mmtsne)
dups <- !duplicated.matrix(as.matrix(test2))
test2 <- test2[dups,]
hbout <- hbeta(test2)
out <- Rtsne(test2, initial_dims = 50, dims = 2, perplexity = hbout$H*.5, theta = 0, max_iter = 1000, verbose = TRUE)

dists <- as.data.frame(out$Y)
pops <- substr(colnames(genotypes[,4:ncol(genotypes)]),1,3)
dists$pop <- pops[dups]
ggplot(dists, aes(x = V1, y = V2, color = pop)) + geom_point() + scale_color_brewer(palette = "Set1")


#get individual genotypes for inversion
haplotypes <- numeric(nrow(div_tsne))
commons <- numeric(ncol(div_tsne))
for(i in 1:ncol(div_tsne)){
  commons[i] <- as.numeric(names(which.max(table(div_tsne[,i]))))
}

for(i in 1:nrow(div_tsne)){
  cgts <- commons[!is.na(div_tsne[i,])]
  calls <- div_tsne[i,]
  calls <- calls[!is.na(calls)]
  hets <- sum(calls == 2)
  choms <- sum(calls == cgts)
  rhoms <- sum(calls != 2 & calls != cgts)
  haplotypes[i] <- names(which.max(c(het = hets,chom = choms,rhom = rhoms)))
}

haps2 <- haplotypes[dups]

dists$hap <- haps2

ggplot(dists, aes(x = V1, y = V2, color = pop, shape = hap)) + geom_point() + scale_color_brewer(palette = "Set1")

#get categories for now
ndups <- div_tsne[dups,]
c1 <- row.names(ndups[dists$V2 < 0 & dists$V1 < -40,])
c2 <- row.names(ndups[dists$V2 < 0 & dists$V1 > -40 & dists$V1 < 20,])
c3 <- row.names(ndups[dists$V2 > -10 & dists$V1 > 20,])
c4 <- row.names(ndups[dists$V2 > 20 & dists$V1 < 0,])

###################
#get consensus sequences for each category.
divgenos <- genotypes[genotypes$snp %in% divsnps,]
divgenos <- format_snps(divgenos, 3, 6, "0000", "00")

#function to get sequence distributions. Takes NN format
summerize_genotypes <- function(x, ecs){
  meta <- x[,1:ecs]
  x <- x[,-c(1:ecs)]
  out <- data.frame(A = numeric(nrow(x)), T = numeric(nrow(x)), 
                    G = numeric(nrow(x)), C = numeric(nrow(x)))
  srow <- function(y){
    as <- c(substr(y, 1, 1), substr(y, 2, 2))
    as <- as[as != "N"]
    tab <- table(as)
    tab <- tab/sum(tab)
    return(tab)
  }
  for (i in 1:nrow(x)){
    tab <- srow(x[i,])
    if("A" %in% names(tab)){out[i,"A"] <- tab[names(tab) == "A"]}
    if("T" %in% names(tab)){out[i,"T"] <- tab[names(tab) == "T"]}
    if("G" %in% names(tab)){out[i,"G"] <- tab[names(tab) == "G"]}
    if("C" %in% names(tab)){out[i,"C"] <- tab[names(tab) == "C"]}
  }
  out <- cbind(meta, out)
  return(out)
}

#run on the correct groups
c1gs <- cbind(cluster = "c1", summerize_genotypes(cbind(divgenos[,1:3], divgenos[,colnames(divgenos) %in% c1,]), 3))
c2gs <- cbind(cluster = "c2", summerize_genotypes(cbind(divgenos[,1:3], divgenos[,colnames(divgenos) %in% c2,]), 3))
c3gs <- cbind(cluster = "c3", summerize_genotypes(cbind(divgenos[,1:3], divgenos[,colnames(divgenos) %in% c3,]), 3))
c4gs <- cbind(cluster = "c4", summerize_genotypes(cbind(divgenos[,1:3], divgenos[,colnames(divgenos) %in% c4,]), 3))
#melt together
gs <- rbind(c1gs, c2gs, c3gs, c4gs)
library(reshape2)
gs_melt <- melt(gs, id.vars = c("cluster", "snp", "position", "group"))

ggplot(gs_melt, aes(x = snp, y = value, color = variable)) +geom_point() + facet_wrap(~cluster, ncol = 1) +
  scale_y_continuous(limits = c(0.0001, 1))




############
#remove last three snps
test3 <- div_tsne
for (i in 1:ncol(test)){
  test3[is.na(test3[,i]),i] <- lavs[i]
}
test3 <- test3[,-c((ncol(test3)-2):ncol(test2))]
dups3 <- !duplicated.matrix(as.matrix(test3))
hbout <- hbeta(test3)
out2 <- Rtsne(test3[dups3,], initial_dims = 50, dims = 2, perplexity = hbout$H*.5, theta = 0, max_iter = 1000, verbose = TRUE)
haps3 <- haplotypes[dups3]
dists2 <- as.data.frame(out2$Y)
dists2$haps <- haps3

ggplot(as.data.frame(dists2), aes(V1, V2, color = haps)) + geom_point()




#pca
pca_out <- prcomp(t(test2))
summary(pca_out)
pcdat <- as.data.frame(pca_out$rotation)

pcdat <- cbind(pop = dists$pop, pcdat)
pcdat$hap <- dists$hap
ggplot(pcdat, aes(x = PC1, y = PC2, color = hap)) + geom_point() + scale_color_brewer(palette = "Set1")

