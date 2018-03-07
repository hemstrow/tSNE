library(snpR)
setwd("../tSNE_data/steelhead/snps/RAPTURE/")
############################################
#read in data
# IBS <- read.table("rapIBS.covMat")
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
# ##############################################
# #run pca
# pca <- PCAfromPA(IBS, 3, c("fam", "pop"))
# pca$plot <- pca$plot + guides(color = FALSE, fill = FALSE)
# pca$plot
# tSNE <- tSNEfromPA(IBS, 3, c("fam", "pop"))
# tSNE$plot <- tSNE$plot + guides(color = FALSE, fill = FALSE)
# tSNE$plot
# 
# 
# ##############################################
# #try raw, with and without removing missings
# IBSr <- read.table("rapIBS.ibs", header = T, stringsAsFactors = F)
# IBSr <- IBSr[,-c(3,4)]
# cnames <- paste0(IBSr[,1], "_", IBSr[,2])
# IBSr <- t(IBSr[,-c(1:2)])
# colnames(IBSr) <- cnames
# IBSr <- cbind(data.frame(pop = pops, fam = fams, samp = samps, stringsAsFactors = F), IBSr, stringsAsFactors = F)
# tSNEr1 <- tSNEfromPA(IBSr, 3, "fam")
# tSNEr1$plot <- tSNEr1$plot + guides(color = FALSE)
# tSNEr1$plot

##############################################
#from genotypes
genos <- read.table("rapIBS.geno")
colnames(genos) <- c("group", "position", paste0(samps, "_", pops, "_", fams))
flt_genos <- filter_snps(genos, 2, 0.05, 0.55, 500, .25)
pa_genos <- format_snps(flt_genos, 2, 7)
pa_genos$samp <- as.character(pa_genos$samp)
meta <- matrix(unlist(strsplit(pa_genos$samp, "_")), nrow(pa_genos), 3, T)
pa_genos <- cbind(data.frame(pop = meta[,2], fam = meta[,3], stringsAsFactors = F), pa_genos, stringsAsFactors = F)

saveRDS(pa_genos, "pa_genos_flex.RDS")

tSNE <- tSNEfromPA(pa_genos, 3, "fam", initial_dims = floor(length(unique(pa_genos$fam))/1.5))
tSNE$plot <- tSNE$plot #+ guides(color = FALSE) 
tSNE$plot

pca <- PCAfromPA(pa_genos, 3, "pop")
pca$plot <- pca$plot + guides(color = FALSE)
pca$plot

#more stringent filtering
flt_genos <- filter_snps(genos, 2, 0.05, 0.55, 500, .5)
pa_genos <- format_snps(flt_genos, 2, 7)
pa_genos$samp <- as.character(pa_genos$samp)
meta <- matrix(unlist(strsplit(pa_genos$samp, "_")), nrow(pa_genos), 3, T)
pa_genos <- cbind(data.frame(pop = meta[,2], fam = meta[,3], stringsAsFactors = F), pa_genos, stringsAsFactors = F)

saveRDS(pa_genos, "pa_genos_strin.RDS")
tSNE <- tSNEfromPA(pa_genos, 3, c("fam"), initial_dims = floor(length(unique(pa_genos$fam))/1.5))
tSNE$plot <- tSNE$plot + guides(color = FALSE)
tSNE$plot

pca <- PCAfromPA(pa_genos, 3, "pop")
pca$plot <- pca$plot
pca$plot
