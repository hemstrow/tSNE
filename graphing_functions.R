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
  
  #check for any duplicates, which need to be removed!
  dups <- which(duplicated(x) | duplicated(x, fromLast=TRUE))
  if(length(dups) > 0){
    warning("Duplicates detected, indices:", dups, "\nRemoving all of these!")
    x <- x[-dups,]
    meta <- meta[-dups,]
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

#'tSNE from genetic data.
#'
#'\code{tSNEfromPA} creates a ggplot object tSNE from presence/absense allelic data using the Barnes-Hut simulation at theta>0 implemented in \code{\link[Rtsne]{Rtsne}}. Returns more data than \code{\link{tsne}}. If individuals which were sequenced at too few loci are to be filtered out, both the mc and counts argument must be provided.
#'
#'See the documentaion for \code{\link[Rtsne]{Rtsne}} for details. Defaults match those of \code{\link[Rtsne]{Rtsne}}. Argument documentation taken from that function.
#'
#'Description of x:
#'    SNP or other allelic data in presence/absence format, as given by \code{\link{format_snps}} option 7. An additional column of population IDs titled "pop" must also be provided.
#'
#' @param x Input presence/absence data, as described in details.
#' @param ecs The number of extra metadata columns at the start of the input. Must be more that two to avoid errors. I really should fix that at some point. Includes the "pop" collumn.
#' @param mc If the data is filtered of poorly sequenced individuals, how many loci must the individuals have to be kept?
#' @param counts Numeric vector containing the number of loci sequenced per individual.
#' @param dims integer, output dimensionality
#' @param initial_dims integer, the number of dimensions retained in the initial PCA step.
#' @param perplex Perplexity parameter, by default found by \code{\link[mmtsne]{hbeta}}, with beta = 1.
#' @param gravity Theta parameter from \code{\link[Rtsne]{Rtsne}}.
#' @param iter Integer. Number of tSNE iterations to perform.
#' @param ... Other arguments, passed to \code{\link[Rtsne]{Rtsne}}.
#'
#' @return A list containing the raw tSNE output and a tSNE plot in the form of a ggplot graphical object. The plot can be changed as usual with ggplot objects.
#'
#' @examples
#' counts <- colSums(ifelse(stickSNPS[, 4:ncol(stickSNPs)] == "NN", FALSE, TRUE))
#' PCAfromPA(stickPA, 2, 2000, counts)
tSNEfromPA <- function(x, ecs, mc = FALSE, counts = FALSE, dims = 2, initial_dims = 50, perplex = FALSE, gravity = 0, iter = 5000, ...){
  #grab metadata and data
  meta <- x[,1:ecs]
  x <- x[,(ecs+1):ncol(x)]
  x <- as.matrix(x)
  
  #filter if requested
  if(mc & counts){
    keeps <- which(counts >= mc)
    x <- x[keeps,]
    meta <- meta[keeps,]
  }
  else if ((mc & !counts) | (!mc & counts)){
    stop("Counts and mc must either both be defined or neither must be.")
  }
  
  #check for any duplicates, which need to be removed!
  dups <- which(duplicated(x) | duplicated(x, fromLast=TRUE))
  if(length(dups) > 0){
    warning("Duplicates detected, indices:", dups, "\nRemoving all of these!")
    x <- x[-dups,]
    meta <- meta[-dups,]
  }
  
  #get perplexity if not provided
  if(!perplex){
    perp <- hbeta(x, beta = 1)
    perplex <- perp$H
  }
  
  #run the tSNE
  tsne.out = Rtsne(x, dims, initial_dims, perplex, gravity, iter, check_duplicates = FALSE, verbose=TRUE, ...)
  
  #saved_tsne2 <- tsne.out
  tsne_plot <- as.data.frame(tsne.out$Y)
  tsne_plot$pop <- meta$pop
  #tsne_colors <- plot_colors[,1]
  #tsne_plot$colors <- tsne_colors
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  
  out <- ggplot(tsne_plot, aes(V1, V2, color = pop))  + geom_point() + theme_bw() +
    theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(),
          panel.grid = element_blank()) +
    guides(color = guide_legend(title="Population"))
  return(list(tSNE = tsne.out, plot = out))
}

