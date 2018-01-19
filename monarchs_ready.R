#import bamlist and annotations
bamlist <- read.table("bamlist.txt")
DNA274 <- "ACAGTG"
DNA275 <- "CAGATC"
plate2 <- read.csv("plate2_MFMB1.csv", header = T, skip = 1)
plate3 <- read.csv("plate3_MFMB2.csv", header = T, skip = 1)

plate2 <- rbind(plate2[,1:3], setNames(plate2[,4:6], colnames(plate2)[1:3]))
plate3 <- rbind(plate3[,1:3], setNames(plate3[,4:6], colnames(plate3)[1:3]))

barcodes <- read.table("RAD barcodes.txt", header = T)

#get sample info

plate2 <- merge(plate2, barcodes, by = "Well")
plate3 <- merge(plate3, barcodes, by = "Well")
sampinfo <- data.frame(plate = substr(bamlist[,1], 8, 13), ind = substr(bamlist[,1], 20, 27))

plate2 <- merge(sampinfo[sampinfo$plate == DNA274,], plate2, by.x = "ind", by.y = "Index")
plate3 <- merge(sampinfo[sampinfo$plate == DNA275,], plate3, by.x = "ind", by.y = "Index")

combplates <- rbind(plate2, plate3)

#import, name, and sort raw genotypes.
raw_genos <- read.table("genotypes.geno", header = F, stringsAsFactors = F)
colnames(raw_genos) <- c("group", "position", paste0(combplates$Pop, 1:nrow(combplates)))
raw_genos <- cbind(snp = 1:nrow(raw_genos), raw_genos)

#filter raw genotypes.
#flt_genos <- filter_snps(raw_genos, 3, T, T, 0.05, l, 0.55, 75)
rg1 <- raw_genos[1:936600,]
rg2 <- raw_genos[936601:nrow(raw_genos),]
remove(raw_genos)

flt_genos_p1 <- filter_snps(rg1, 3, T, T, 0.05, hf_hets =  0.55, min_ind =  100)
flt_genos_p2 <- filter_snps(rg2, 3, T, T, 0.05, hf_hets = 0.55, min_ind = 100)
flt_genos <- rbind(flt_genos_p1, flt_genos_p2)
write.table(flt_genos, "flt_genotypes2.txt", col.names = T, quote = F, sep = "\t")


sub_snps <- sample(1:nrow(flt_genos), 8000)
pa_snps <- format_snps(flt_genos[sub_snps,], 3, 8, l_names = flt_genos$snp[sub_snps])



#pca
pca_raw <- prcomp(pa_snps[,-1])
pca <- as.data.frame(pca_raw$x)
pca$pop <- combplates$Pop
pca$pop <- substr(pca$pop, 1, 3)
pca[pca$pop == "NZ2",]$pop <- "NZ"
pca[pca$pop == "NZR",]$pop <- "NZ"

#which samples are crappy?
counts <- ifelse(flt_genos[sub_snps, 4:ncol(flt_genos)] == "NN", FALSE, TRUE)
counts <- colSums(counts)

gpa_snps <- pa_snps[counts >= 4000,]

library(ggplot2)
ggplot(pca, aes(PC1, PC2, color = pop)) + geom_point() + scale_color_brewer(palette = "Set3")

library(Rtsne)
hbout <- hbeta(gpa_snps[,-1])
out <- Rtsne(gpa_snps[,-1], initial_dims = 50, dims = 2, perplexity = hbout$H, theta = 0, max_iter = 5000, verbose = TRUE)

dists <- as.data.frame(out$Y)
dists$pop <- combplates$Pop[counts >= 4000]
ggplot(dists, aes(x = V1, y = V2, color = pop)) + geom_point() + theme_bw()


#PCA with the filtered data
gpca_raw <- prcomp(gpa_snps[,-1])
gpca <- as.data.frame(gpca_raw$x)
gpca$pop <- combplates$Pop[counts >= 4000]
gpca$pop <- substr(gpca$pop, 1, 3)
ggplot(gpca, aes(PC1, PC2, color = pop)) + geom_point() + theme_bw()
