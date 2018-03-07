setwd("~/GitHub/tSNE_data/sturgeon/")
library(readr)
library(snpR)
###############################################
#import and prepare data

# Fr.genos <- read_delim("Fraser River data for t-SNE test.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# 
# #need to make a nested list containing locus information
# lnames <- unlist(strsplit(colnames(Fr.genos)[2:ncol(Fr.genos)], "_"))
# lnames <- lnames[nchar(lnames) > 3]
# lnames <- unique(lnames)
# lopts <- list()
# for(i in 1:length(lnames)){
#   tcols <- grep(paste0(lnames[i], "_"), colnames(Fr.genos))
#   tcols <- sort(c(which(colnames(Fr.genos) == lnames[i]),tcols))
#   topts <- as.numeric(Fr.genos[1,tcols])
#   lopts[[i]] <- list(locus = lnames[i], alleles = topts, cols = tcols)
# }
# 
# #function to get this data into presence absence. Returns this plus a count of missing data per ind
# #Arguments:
# # x: input data. make sure that the FIRST row is the first row of DATA (not header info)!
# #    Note that the program expects absences to be noted by NA. Anything at all in the cell will be noted as a presence.
# # ecs: number of header COLUMNS before the first column of data.
# # linfo: locus info. A nested list. At the outer level, each element is a locus. Each inner level list
# #        contains locus information, in a list as given by:
# #        list(locus = "locus name, alleles = c("a1", "a2", "a3", ...), cols = numeric vector of the columns containing the allele presence absences for the locus))
# format_mp_pa <- function(x, ecs, linfo){
#   #grab header info and convert the NAs to zeros and the others to 1
#   x[,(ecs+1):ncol(x)] <- ifelse(is.na(x[,(ecs+1):ncol(x)]), 0, 1)
# 
#   #initialize missing count
#   miss <- numeric(length(nrow(x)))
# 
#   #get observed allele frequencies
#   for(i in 1:length(linfo)){
#     cat("Formating locus:", linfo[[i]]$locus, "\n")
#     xi <- x[,linfo[[i]]$cols]
#     xm <- colMeans(xi)
#     miss <- miss + (rowSums(xi) == 0)
#     xi[rowSums(xi) == 0,] <- xm
#     x[,linfo[[i]]$cols] <- xi
#     colnames(x)[linfo[[i]]$cols] <- paste0(linfo[[i]]$locus, "_", linfo[[i]]$alleles)
#   }
#   return(list(pam = x, miss = miss))
# }
# 
# pa.geno <- format_mp_pa(Fr.genos[2:nrow(Fr.genos),], 1, lopts)
# saveRDS(pa.geno, "formatted_frasier_pa.RDS")

pa.geno <- readRDS(pa.geno, "formatted_frasier_pa.RDS")

mcount <- pa.geno$miss
pa.geno <- pa.geno$pam



##########################################
#run tSNE and PCA

minmiss <- 0
rdat <- as.data.frame(pa.geno[mcount <= minmiss,])
rdat <- cbind(Population = gsub("[0-9]+", "", rdat$X1), rdat)

pca <- PCAfromPA(rdat, 2, do.plot = "Population", c.dup = T)
tSNE <- tSNEfromPA(rdat, 2, do.plot = "Population", c.dup = T)
