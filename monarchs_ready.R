setwd("~/2017-2018/Monarchs/")
genotypes <- read.table("flt_genotypes.txt", header = T, stringsAsFactors = F)
genotypes_meta <- genotypes[,1:3]
genotypes_data <- genotypes[,4:ncol(genotypes)]

#order the data by pop name and perpare l for formatting.
abnames <- substr(colnames(genotypes_data), 1, 3)
names.tab <- table(abnames)
colnames(genotypes_data) <- abnames
colnames(genotypes_data)[colnames(genotypes_data) == "NZ2" | colnames(genotypes_data) == "NZR"] <- "NZ"
genotypes_data <- genotypes_data[,order(colnames(genotypes_data))]
genotypes <- cbind(genotypes_meta, genotypes_data)
remove(genotypes_data, genotypes_meta)

names.tab <- c(names.tab, NZ = sum(names.tab["NZR"] + names.tab["NZ2"]))
names.tab <- names.tab[-which(names(names.tab) %in% c("NZR", "NZ2"))]
names.tab <- names.tab[order(names(names.tab))]

l <- list(names(names.tab), as.numeric(names.tab))

#sub_snps <- sample(nrow(genotypes), 10000, replace = F)
#gsub <- genotypes[sub_snps,]

#pa_snps <- format_snps(gsub, 3, 7, l_names = gsub[,1])
#writeLines(as.character(sub_snps), "subset_indices.txt", sep = "\t")
#write.table(pa_snps, "pa_subset.txt", col.names = T, quote = F, sep = "\t")


pa_snps <- read.table("pa_subset.txt", header = T, stringsAsFactors = F)
sub_snps <- readLines("subset_indices.txt")
sub_snps <- as.numeric(unlist(strsplit(sub_snps, "\t")))
pops <- gsub("\\.[0-9]+", "", pa_snps$samp)

#
mc <- 8000

#pca
pca_raw <- prcomp(pa_snps[,-1])
pca <- as.data.frame(pca_raw$x)
pca$pop <- pops

#which samples are crappy?
counts <- ifelse(genotypes[sub_snps, 4:ncol(genotypes)] == "NN", FALSE, TRUE)
counts <- colSums(counts)

gpa_snps <- pa_snps[counts >= mc,]

library(ggplot2)
ggplot(pca, aes(PC1, PC2, color = pop)) + geom_point() + scale_color_brewer(palette = "Set3")

library(Rtsne)
library(mmtsne)
hbout <- hbeta(gpa_snps[,-1])
out <- Rtsne(gpa_snps[,-1], initial_dims = 50, dims = 2, perplexity = hbout$H, theta = 0, max_iter = 5000, verbose = TRUE)

dists <- as.data.frame(out$Y)
dists$pop <- pops[counts >= mc]
ggplot(dists, aes(x = V1, y = V2, color = pop)) + geom_point() + theme_bw()


#PCA with the filtered data
gpca_raw <- prcomp(gpa_snps[,-1])
gpca <- as.data.frame(gpca_raw$x)
gpca$pop <- pops[counts >= mc]
gpca$pop <- substr(gpca$pop, 1, 3)
ggplot(gpca, aes(PC1, PC2, color = pop)) + geom_point() + theme_bw()
