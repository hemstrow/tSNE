library(tsne)
library(ggplot2)
library(mmtsne)
library(FNN)
library(Rtsne)


# Reading in my Colony data (1 present 0 absent) and turning it into a matrix
sample_table = read.table("Training_Genotypes_tsne.txt",sep='\t', header=FALSE)
sturg_matrix = as.matrix(sample_table, dimnames = NULL)
names <- read.table("Names.txt", sep='\t', header = FALSE)
plot_colors <- read.table("colors.txt", sep='\t', header = FALSE)


# tsne fxn using Rtsne, which uses Barnes-Hut simulation at theta>0 and returns more data than tsne()
# Rtsne() is also much faster than tsne()!
tsne_fxn <- function(data, dims = FALSE, initial_dims, perplex = FALSE, gravity = FALSE, iter = FALSE){
  if(!dims){
      dims <- 2
  }
  if(!perplex){
    perp <- hbeta(sturg_matrix, beta = 1)
    perplex <- perp$H
  }
  if(!iter){
    iter <- 5000
  }
  if(!gravity){
    gravity <- 0
  }
  tsne.out = Rtsne(data, dims = dims, perplexity = perplex, theta = gravity, max_iter = iter, verbose=TRUE)
  saved_tsne2 <<- tsne.out
  tsne_plot <- as.data.frame(tsne.out$Y)
  tsne_colors <- plot_colors[,1]
  tsne_plot$colors <- tsne_colors
  require(ggplot2)
  ggplot(tsne_plot, aes(tsne_plot[,1],tsne_plot[,2], color = colors))  + geom_point() + geom_text(aes(label=names))
}

# running the t-SNE function
tsne_fxn(data = sturg_matrix)


# hbeta estimates perplexity and probability values for my data, from the mmtsne package
hbeta(sturg_matrix, beta = 1)


# plotting a quick PCA
# genetic distance from GenALex?
pca_matrix <- sturg_matrix[,-c(55,74)]
pca <- principal(pca_matrix)
loadings <- as.data.frame(pca$loadings[,1:2])
ggplot(loadings, aes(loadings[,1],loadings[,2]))  + geom_point(aes(loadings[,1],loadings[,2])) + geom_text(aes(label=names))

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
