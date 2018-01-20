# little script to download the files of a package. Untar unzips it, download downloads it
untar(download.packages(pkgs = "Rtsne", 
                  destdir = ".",
                  type = "source")[,2])


library(tsne)
library(ggplot2)
library(ggfortify)
library(gplots)
library(mmtsne)
library(FNN)
library(Rtsne)
library(proxy)
library(kernlab)
library(KRLS)
library(randomForest)

# Reading in my Colony data (1 present 0 absent) and turning it into a matrix
sample_table = read.table("Training_Genotypes_tsne.txt",sep='\t', header=FALSE)
sturg_matrix = as.matrix(sample_table, dimnames)
names <- read.table("Names.txt", sep='\t', header = FALSE)
plot_colors <- read.table("colors.txt", sep='\t', header = FALSE)
# Making a 2d TSNE plot

# tsne fxn using Rtsne, which uses Barnes-Hut simulation at theta>0 and returns more data than tsne()
# Rtsne() is also much faster than tsne()!
tsne_fxn <- function(dims, initial_dims, perplex, iter, pca = TRUE){
  tsne.out = Rtsne(sturg_matrix, dims = dims, perplexity = perplex, theta = 0, max_iter = iter, verbose=TRUE)
  saved_tsne2 <<- tsne.out
  tsne_plot <- as.data.frame(tsne.out$Y)
  tsne_colors <- plot_colors[,1]
  tsne_plot$colors <- tsne_colors
  require(ggplot2)
  ggplot(tsne_plot, aes(tsne_plot[,1],tsne_plot[,2], color = colors))  + geom_point() + geom_text(aes(label=names))
}

# running the t-SNE function
tsne_fxn(2,100,8.77,1000)

# Playing with a k nearest neighbors algorithm
cl <- factor(c(rep("fullsib",25), rep("halfsib",25), rep("nosib",4)))
knn(saved_tsne$Y, saved_tsne2$Y, cl, k=5, prob=TRUE)

saved_tsne <- Rtsne(sturg_matrix, dims = 2, perplexity = 8.77, theta = 0.5, max_iter = 1000, verbose=TRUE)

View(saved_tsne$Y)
testframe <- data.frame(saved_tsne$Y)
tsne_colors <- plot_colors[,1]
testframe$colors <- tsne_colors
require(ggplot2)
ggplot(testframe, aes(testframe[,1],testframe[,2], color = colors))  + geom_point() + geom_text(aes(label=names))

write.table(saved_tsne$Y, "tsne.csv")

# hbeta estimates perplexity and probability values for my data, from the mmtsne package
hbeta(sturg_matrix, beta = 1)
num_maps <- 15

multiple_maps <- mmtsne(sturg_matrix, no_maps = num_maps, no_dims = 4, perplexity = 8.77, max_iter = 2000)
par(mfrow=c(1,2))
plot(multiple_maps$Y[,,1])
plot(multiple_maps$Y[,,2])

multiple_plot <- as.data.frame(multiple_maps)
par(mfrow=c(1,2))
plot(multiple_plot$Y.1, multiple_plot$Y.2)
ggplot(multiple_plot, aes(multiple_plot$Y.7, multiple_plot$Y.8)) + geom_point(aes(multiple_plot$Y.7, multiple_plot$Y.8)) + geom_text(aes(label=names))
ggplot(multiple_plot, aes(multiple_plot$Y.3, multiple_plot$Y.4)) + geom_point(aes(multiple_plot$Y.3, multiple_plot$Y.4)) + geom_text(aes(label=names))


# calculating a similarity matrix based on the student-t kernel method with 1 DF
simil_matrix <- function(X)
  {
  return(1/(1+as.matrix(dist(X))^2))
  }

kullback_fxn <- function(original_data, tsne_data){
  out_matrix <- matrix(0, nrow(original_data), nrow(original_data))
  for(i in 1:num_maps){
    student <- simil_matrix(tsne_data$Y[,,i])
    weight <- as.matrix(tsne_data$weights[,i])
    weighted <- outer(student, weight, FUN = "*")
    weightmatrix <- matrix(weighted, nrow = nrow(weighted), ncol = ncol(weighted))
    out_matrix <- out_matrix + weightmatrix
  }
  out_matrix
  diag(out_matrix) <- 0
  qmatrix <- out_matrix/sum(out_matrix)
  pmatrix <- as.matrix(simil(sturg_matrix, diag = T, upper = T))
  diag(pmatrix) <- 0
  kullback <- KL.divergence(qmatrix, pmatrix, k = 7)
  return(kullback)
}  

kullback_fxn(sturg_matrix, multiple_maps)


nearest_fxn <- function(nearest, q, p){
  # Finding K nearest neighbors for the neighborhood preservation ratio
  knn_output <- get.knn(q, k=nearest, algo="kd_tree")
  knn_input <- get.knn(p, k=nearest, algo="kd_tree")
  # these are just the individuals and their nearest neighbors for the NPR
  outputs <- as.matrix(knn_output$nn.index)
  inputs <- as.matrix(knn_input$nn.index)
  intersection <- 0
  # beginning of a function to iterate through each matrix, comparing rows
  for(i in 1:nrow(outputs)){
      intersection <- intersection + length(intersect(outputs[i,], inputs[i,]))
  }
  # finding the neighborhood preservation ratio
  npr <- intersection/length(outputs)
  return(npr)
}

nearest_fxn(5, qmatrix, pmatrix)

intersect_list <- list()
for(i in 1:nrow(outputs)){
  intersect_list[[i]] <- intersect(outputs[i,], inputs[i,])
}
intersect_matrix <- as.matrix(intersect_list)

# Making a 3d TSNE
third_tsne.out = tsne(sturg_matrix, initial_config = NULL, k = 3)
write.csv(third_tsne.out, file="3d_tsne.txt", append = FALSE, quote = TRUE, sep = "\t")

sample_names <- read.table("BLCJ2016_Repat.txt", sep='\t', header=FALSE)
third_tsne.out = tsne(sturg_matrix, initial_config = NULL, k = 3, perplexity = 8.77)



# creating a Euclidean distance matrix between samples from the 3d TSNE
distance_matrix <- dist(third_tsne.out, method = "euclidean", diag = T)

# Writing that matrix out to a CSV
matrix <- as.matrix(distance_matrix)
write.csv(matrix, file = "3d_dist_matrix.csv", append = TRUE, quote = TRUE, sep = "\t")

# Plotting that TSNE in 3d
library(scatterplot3d)
scatterplot3d(third_tsne.out)

# plotting a quick PCA
# genetic distance from GenALex?
pca_matrix <- sturg_matrix[,-c(55,74)]
pca <- principal(pca_matrix)
loadings <- as.data.frame(pca$loadings[,1:2])
ggplot(loadings, aes(loadings[,1],loadings[,2]))  + geom_point(aes(loadings[,1],loadings[,2])) + geom_text(aes(label=names))

# Playing with a k nearest neighbors algorithm
cl <- factor(c(rep("fullsib",25), rep("halfsib",25), rep("nosib",4)))
knn(saved_tsne$Y, pca$rotation[,1:2], cl, k=5, prob=TRUE)


pca_dist <- dist(loadings, method = "euclidean", diag = T)
pca_matrix <- as.matrix(pca_dist)
write.csv(pca_matrix, file = "2d_PCA_dist.csv")



# playing with the mmtsne and a known data set
data("iris")
# Estimate a mmtsne model with 2 maps, 2 dimensions each
model <- mmtsne(iris[,1:4], no_maps=2, max_iter=100)
# Plot the results side-by-side for inspection
# Points scaled by map proportion weights plus constant factor
par(mfrow=c(1,2))
plot(model$Y[,,1], col=iris$Species, cex=model$proportions[,1] + .2)
plot(model$Y[,,2], col=iris$Species, cex=model$proportions[,2] + .2)
par(mfrow=c(1,1))

# gen_distance
GD_table = read.table("Sterling_GD.txt",sep='\t', header=FALSE)
GD_matrix = data.matrix(GD_table)
GD_pca <- prcomp(GD_matrix)
GD_loadings <- as.data.frame(GD_pca$rotation[,1:2])
ggplot(GD_loadings, aes(GD_loadings[,1],GD_loadings[,2]))  + geom_point(aes(GD_loadings[,1],GD_loadings[,2])) + geom_text(aes(label=names))
