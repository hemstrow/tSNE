library(Rtsne)
library(mmtsne)
library(ggplot2)

sample_table <- read.table("groupedfp90_genepop.txt",sep='', header=FALSE, colClasses = "character")

sample_table <- t(sample_table)
pops <- sample_table[1,]
sample_table <- sample_table[-1,]

sample_table <- sample_table[-c(1:3),]

lnam <- c("Oki23", "Ssa407", "mSsa408", "Ots209", "OtsG249b", "OtsG85", "Omy27", "Omy1001", "Ots243", "Ots409", "OtsG3", "Ots212", "Omm1087")
formatted_table <- format_snps(sample_table, 0, output = 8, 
                               input_form = "msat_2", miss = "00", l_names = lnam, 
                               interp_miss = TRUE)


hbeta(formatted_table[,-1], beta = 1)

tsne <- Rtsne(formatted_table[,-1], dims = 2, initial_dims = 50, perplexity = 25, theta = 0, 
              check_duplicates = FALSE, pca = TRUE, max_iter = 1000, verbose = TRUE, 
              is_distance = FALSE)
plot(tsne$Y)

plot_data <- data.frame(tsne$Y)

pca <- prcomp(t(formatted_table[,-1]))
plot(pca$rotation)
p_dat <- as.data.frame(pca$rotation)
p_dat$pop <- metadata

metadata <- read.table("metadata.txt",sep='\t', header=FALSE)
metadata <- metadata[,1]
#metadata <- metadata[dupes]

#dupes <- !duplicated.matrix(as.matrix(sample_table))

plot_data$pop <- metadata

# plotting the t-SNE
ggplot(plot_data, aes(plot_data[,1], plot_data[,2], color = pop)) + geom_point()

# plotting PCA data
ggplot(p_dat, aes(PC1, PC2, color = pop)) + geom_point()

multiple <- mmtsne(formatted_table[,-1])

