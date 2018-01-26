library(ggplot2)
library(Rtsne)
library(mmtsne)

########################################################################################################
#set up data


#get whole genome info for picking out crappy samples.
setwd("~/GitHub/tSNE_data/monarch/")
genotypes <- read.table("flt_genotypes.txt", header = T, stringsAsFactors = F)

#order the data by pop name get pop info. Do this only once.
#genotypes_meta <- genotypes[,1:3]
#genotypes_data <- genotypes[,4:ncol(genotypes)]
#abnames <- substr(colnames(genotypes_data), 1, 3)
#names.tab <- table(abnames)
#colnames(genotypes_data) <- abnames
#colnames(genotypes_data)[colnames(genotypes_data) == "NZ2" | colnames(genotypes_data) == "NZR"] <- "NZ"
#genotypes_data <- genotypes_data[,order(colnames(genotypes_data))]
#genotypes <- cbind(genotypes_meta, genotypes_data)
#remove(genotypes_data, genotypes_meta)
#write.table(genotypes, "flt_genotypes.txt", quote = F, col.names = T, sep = "\t")



#subset and reformat. Don't need to do this more than once.
#names.tab <- c(names.tab, NZ = sum(names.tab["NZR"] + names.tab["NZ2"]))
#names.tab <- names.tab[-which(names(names.tab) %in% c("NZR", "NZ2"))]
#names.tab <- names.tab[order(names(names.tab))]
#l <- list(names(names.tab), as.numeric(names.tab))
#sub_snps <- sample(nrow(genotypes), 10000, replace = F)
#gsub <- genotypes[sub_snps,]

#pa_snps <- format_snps(gsub, 3, 7, l_names = gsub[,1])
#writeLines(as.character(sub_snps), "subset_indices.txt", sep = "\t")
#write.table(pa_snps, "pa_subset.txt", col.names = T, quote = F, sep = "\t")

#read in formated data and attach metadata.
pa_snps <- read.table("pa_subset.txt", header = T, stringsAsFactors = F)
sub_snps <- readLines("subset_indices.txt")
sub_snps <- as.numeric(unlist(strsplit(sub_snps, "\t")))
pops <- gsub("\\.[0-9]+", "", pa_snps$samp)
pa_snps$pop <- pops
pa_snps <- pa_snps[,c(1,ncol(pa_snps), 2:(ncol(pa_snps) - 1))]

#generate counts of the number of sequenced loci per individual for this subset of data.
counts <- ifelse(genotypes[sub_snps, 4:ncol(genotypes)] == "NN", FALSE, TRUE)
counts <- colSums(counts)


#################################################################
#run the PCA and tSNE
#PCA:
pca <- PCAfromPA(pa_snps, 2, 8000, counts)

#tSNE
tSNE <- tsnefromPA(pa_snps, 2, 8000, counts)

  


