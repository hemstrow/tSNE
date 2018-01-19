library(Rtsne)
library(mmtsne)
library(ggplot2)

setwd("C:/Users/mattt/Dropbox/Cougar River Family Data")

sample_table <- read.table("14JUL2015_cougar_juvenile_master.txt",sep='', header=TRUE, colClasses = "character")
samps <- sample_table[,1]

# Transposing the sample table and cutting out metadata
sample_table <- t(sample_table)
sample_table <- sample_table[-1,]
sample_table <- sample_table[1:22,]
colnames(sample_table) <- samps

# empty matrix to fill in the reformat function's input data
format_input <- matrix(NA, nrow = nrow(sample_table)/2, ncol = ncol(sample_table))


# loop to turn msat data from one allele per cell to both alleles per cell (diploid data)
j = 1
i = 1
while (i < nrow(sample_table)){
  format_input[j,] <- paste0(sample_table[i,],sample_table[i+1,])
  i <- i + 2
  j = j + 1
}

# Locus names
lnam <- c("Ots201b", "Ots209", "Ots249", "Ots253b", "Ots215", "OtsG311", "OtsG409", "Ots211", "Ots208b", "Ots212", "Ots515")

# Throwing loci and sample names in for funsies
colnames(format_input) <- samps
rownames(format_input) <- lnam

# Will's reformatting function
formatted_table <- format_snps(format_input, 0, output = 8, 
                               input_form = "msat_3", miss = "000", l_names = lnam, 
                               interp_miss = TRUE)

# Pulling out juveniles from year class 2011, searching for "MR11"
juvs <- grep("MR11",formatted_table[,1])
# Creating a new table of those 2010 juvs from the master table
single_year_class <- formatted_table[juvs,]

# Running the t-SNE on this subset data
perplexity <- hbeta(single_year_class[,-1], beta = 1)

tsne <- Rtsne(single_year_class[,-1], dims = 2, initial_dims = 50, perplexity = perplexity$H, theta = 0.5, 
              check_duplicates = FALSE, pca = TRUE, max_iter = 1000, verbose = TRUE, 
              is_distance = FALSE)

plot(tsne$Y)

sibdata <- read.table("aj.pedigree.10.mendel.checked.txt", sep = "\t", header=T)

plot_data <- data.frame(tsne$Y)
rownames(plot_data) <- samps[juvs]
# adding a new column to plot_data to make merging easier
plot_data$off <- row.names(plot_data)
# combining plot data and sibdata to have parentage with t-SNE coordinates
comb <- merge(x = sibdata, plot_data, by = "off")


# plotting the t-SNE. Change color = dad or mom as needed
# geompoint fill makes the inside of the points the mom's color, and the outline the dad's
# theme legend position removes the imcomprehensible legend
ggplot(comb, aes(comb$X1, comb$X2, color = dad)) + geom_point(aes(fill=mom), pch=21, size=2.5) + theme(legend.position="none")


write.table(comb, file = "2011_juv_tsne.csv")

pca <- prcomp(t(formatted_table[,-1]))
plot(pca$rotation)
p_dat <- as.data.frame(pca$rotation)
p_dat$pop <- metadata

metadata <- metadata[,1]
#metadata <- metadata[dupes]

#dupes <- !duplicated.matrix(as.matrix(sample_table))

plot_data$pop <- formatted_table["samp"]


# plotting PCA data
ggplot(p_dat, aes(PC1, PC2, color = pop)) + geom_point()

multiple <- mmtsne(formatted_table[,-1])
