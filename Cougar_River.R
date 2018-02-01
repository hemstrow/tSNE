library(ggplot2)
library(snpR)
setwd("../tSNE_data/chinook/")

##########################################3
#prepare data


# sample_table <- read.table("14JUL2015_cougar_juvenile_master.txt",sep='', header=TRUE, colClasses = "character", stringsAsFactors = F)
# samps <- sample_table[,1]
# 
# # Transposing the sample table and cutting out metadata
# sample_table <- t(sample_table)
# sample_table <- sample_table[-1,]
# sample_table <- sample_table[1:22,]
# colnames(sample_table) <- samps
# 
# # empty matrix to fill in the reformat function's input data
# format_input <- matrix(NA, nrow = nrow(sample_table)/2, ncol = ncol(sample_table))
# 
# 
# # loop to turn msat data from one allele per cell to both alleles per cell (diploid data)
# j = 1
# i = 1
# while (i < nrow(sample_table)){
#   format_input[j,] <- paste0(sample_table[i,],sample_table[i+1,])
#   i <- i + 2
#   j = j + 1
# }
# 
# # Locus names
# lnam <- c("Ots201b", "Ots209", "Ots249", "Ots253b", "Ots215", "OtsG311", "OtsG409", "Ots211", "Ots208b", "Ots212", "Ots515")
# 
# # Throwing loci and sample names in
# colnames(format_input) <- samps
# format_input <- cbind(loci = 1:nrow(format_input), marker = lnam, format_input)
# 
# # Will's reformatting function
# pa_genos <- format_snps(format_input, 2, output = 7, 
#                                input_form = "msat_3", miss = "000", lnames = lnam)
# 
# # Creating a new table of those 2011 juvs from the master table
# single_year_class <- pa_genos[grep("MR11",pa_genos[,1]),]
# 
# 
# #merge with sib data.
# sibdata <- read.table("aj.pedigree.10.mendel.checked.txt", sep = "\t", header=T)
# colnames(single_year_class)[1] <- "off"
# comb <- merge(x = sibdata, single_year_class, by = "off")
# 
# #save
# saveRDS(comb, "pa11_genos.RDS")


pa11_genos <- readRDS("pa11_genos.RDS")
#########################################################
#plot the PCA

pca <- PCAfromPA(pa11_genos, 4, c("mom", "dad"), TRUE)
pca$plot <- pca$plot + guides(color = FALSE, fill = FALSE)

tSNE <- tSNEfromPA(pa11_genos, 4, c("mom", "dad"), c.dup = TRUE, iter = 1000)
tSNE$plot <- tSNE$plot + guides(color = FALSE, fill = FALSE)
