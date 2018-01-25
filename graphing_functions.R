#'PCAs from genetic data.
#'
#'\code{PCAfromPA} creates a ggplot object PCA from presence/absense allelic data. If individuals which were sequenced at too few loci are to be filtered out, both the mc and counts argument must be provided.
#'
#'Description of x:
#'    SNP or other allelic data in presence/absence format, as given by \code{\link{format_snps}} option 7. An additional column of population IDs titled "pop" must also be provided.
#'
#' @param x Input presence/absence data, as described in details.
#' @param ecs The number of extra metadata columns at the start of the input. Must be more that two to avoid errors. I really should fix that at some point. Includes the "pop" collumn.
#' @param mc If the data is filtered of poorly sequenced individuals, how many loci must the individuals have to be kept?
#' @param counts Numeric vector containing the number of loci sequenced per individual.
#'
#' @return A list containing the raw PCA output and a PCA plot in the form of a ggplot graphical object. The plot can be changed as usual with ggplot objects.
#'
#' @examples
#' counts <- colSums(ifelse(stickSNPS[, 4:ncol(stickSNPs)] == "NN", FALSE, TRUE))
#' PCAfromPA(stickPA, 2, 2000, counts)
PCAfromPA <- function(x, ecs, mc = FALSE, counts = FALSE){
  #grab metadata and data
  meta <- x[,1:ecs]
  x <- x[,(ecs+1):ncol(x)]
  
  #filter if requested
  if(mc & counts){
    keeps <- which(counts >= mc)
    x <- x[keeps,]
    meta <- meta[keeps,]
  }
  else if ((mc & !counts) | (!mc & counts)){
    stop("Counts and mc must either both be defined or neither must be.")
  }
  
  pca_r <- prcomp(x)
  pca <- as.data.frame(pca_r$x)
  pca$pop <- meta$pop
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  
  out <- ggplot(pca, aes(PC1, PC2, color = pop)) + geom_point() + theme_bw() + 
    scale_color_manual(values = cbbPalette) + guides(color = guide_legend(title="Population"))
  
  return(list(raw = pca_r, plot = out))
}


# tsne fxn using Rtsne, which uses Barnes-Hut simulation at theta>0 and returns more data than tsne()
tsnefromPA <- function(data, dims = FALSE, initial_dims, perplex = FALSE, gravity = FALSE, iter = FALSE, ecs , mc = FALSE, counts = FALSE){
  #grab metadata and data
  meta <- data[,1:ecs]
  data <- data[,(ecs+1):ncol(data)]
  
  #filter if requested
  if(mc & counts){
    keeps <- which(counts >= mc)
    data <- data[keeps,]
    meta <- meta[keeps,]
  }
  else if ((mc & !counts) | (!mc & counts)){
    stop("Counts and mc must either both be defined or neither must be.")
  }
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
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  ggplot(tsne_plot, aes(tsne_plot[,1],tsne_plot[,2], color = colors))  + geom_point() + theme_bw()
}

