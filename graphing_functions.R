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
PCAfromPA <- function(x, ecs, do.plot = "pop", c.dup = FALSE, mc = FALSE, counts = FALSE){
  library(ggplot2)
  
  #grab metadata and data
  meta <- x[,1:ecs]
  x <- x[,(ecs+1):ncol(x)]
  
  ##############################
  #sanity checks...
  if((mc != FALSE & !counts != FALSE) | (!mc != FALSE & counts != FALSE)){
    stop("Counts and mc must either both be defined or neither must be.")
  }
  
  if(is.character(do.plot)){
    if(length(do.plot) > 2){
      stop("Only two plotting variables supported. For more, set do.plot to FALSE and plot manually.\n")
    }
  }
  else if (do.plot != FALSE){
    stop("do.plot must be either FALSE or between 1 and 2 variables to plot by.")
  }
  else if (!all(do.plot) %in% colnames(meta)){
    stop("Plotting variables specified in do.plot must match column names present in the metadata of x.\n")
  }
  
  ############################################
  #filter if requested
  if(mc != FALSE & counts != FALSE){
    keeps <- which(counts >= mc)
    x <- x[keeps,]
    meta <- meta[keeps,]
  }
  
  #check for any duplicates, which need to be removed!
  if(c.dup){
    cat("Checking for duplicates...\n")
    dups <- which(duplicated(x) | duplicated(x, fromLast=TRUE))
    if(length(dups) > 0){
      cat("Duplicates detected, indices:", dups, "\nRemoving all of these!\n")
      x <- x[-dups,]
      meta <- meta[-dups,]
    }
  }
  
  cat("Preparing pca...\n")
  pca_r <- prcomp(as.matrix(x))
  pca <- as.data.frame(pca_r$x) #grab the PCA vectors.
  pca <- cbind(meta, pca)  #add metadata that is present in the input.
  
  
  #return just the pca data if the plot isn't requested, which is useful for complex plots.
  if(!is.character(do.plot)){
    return(list(raw = pca_r, pca = pca))
  }
  
  ################################################
  #construct plot.
  cat("Preparing plot...\n")
  
  #Categories (pops, fathers, mothers, ect.) are given in do.plot argument. Supports up to two!
  #make the base plot, then add categories as color and fill.
  out <- ggplot(pca, aes(PC1, PC2)) + theme_bw() #initialize plot
  long <- FALSE #are any of the sets of categories really long?
  
  
  #add variables.
  if(length(do.plot) == 1){
    v1 <- pca[,which(colnames(pca) == do.plot[1])] #get the factors
    v1u <- length(unique(v1)) #number of categories
    
    out <- out + geom_point(aes(color = v1)) #add the factor
    
    if(v1u >= 8){long <- TRUE} #are there too many categories to color with the cbb palette?
  }
  
  if(length(do.plot) == 2){
    v1 <- pca[,which(colnames(pca) == do.plot[1])]
    v2 <- pca[,which(colnames(pca) == do.plot[2])]
    v1u <- length(unique(v1))
    v2u <- length(unique(v2))
    
    out <- out + geom_point(aes(color = v1, fill = v2), pch = 21, size = 2.5, stroke = 1.25)
    
    if(v1u >= 8 | v2u >= 8){long <- TRUE}
  }
  
  #use the color blind friendly palette if possible!
  if(!long){
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    out <- out + scale_color_manual(values = cbbPalette)
  }
  
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
tSNEfromPA <- function(x, ecs, do.plot = "pop", dims = 2, initial_dims = 50, 
                       perplex = FALSE, gravity = 0, iter = 5000, 
                       c.dup = FALSE, mc = FALSE, counts = FALSE, ...){
  library(ggplot2)
  #grab metadata and data
  meta <- x[,1:ecs]
  x <- x[,(ecs+1):ncol(x)]
  x <- as.matrix(x)
  
  ##############################
  #sanity checks...
  if((mc != FALSE & !counts != FALSE) | (!mc != FALSE & counts != FALSE)){
    stop("Counts and mc must either both be defined or neither must be.")
  }
  
  if(c.dup == FALSE){
    warning("If there are duplicates in x, expect wierd results! Set c.dup to TRUE to check.\n")
  }
  
  if(is.character(do.plot)){
    if(length(do.plot) > 2){
      stop("Only two plotting variables supported. For more, set do.plot to FALSE and plot manually.\n")
    }
  }
  else if (do.plot != FALSE){
    stop("do.plot must be either FALSE or between 1 and 2 variables to plot by.")
  }
  else if (!all(do.plot) %in% colnames(meta)){
    stop("Plotting variables specified in do.plot must match column names present in the metadata of x.\n")
  }
  
  ##############################
  #filter if requested
  if(mc & counts){
    keeps <- which(counts >= mc)
    x <- x[keeps,]
    meta <- meta[keeps,]
  }
  
  
  if(c.dup){
    #check for any duplicates, which need to be removed!
    cat("Checking for duplicates...\n")
    dups <- which(duplicated(x) | duplicated(x, fromLast=TRUE))
    if(length(dups) > 0){
      cat("Duplicates detected, indices:", dups, "\nRemoving all of these!\n")
      x <- x[-dups,]
      meta <- meta[-dups,]
    }
  }
  
  #get perplexity if not provided
  if(!perplex){
    cat("Estimating perplexity...\n")
    perp <- mmtsne::hbeta(x, beta = 1)
    perplex <- perp$H
  }
  
  #run the tSNE
  cat("Running tSNE...\n")
  tsne.out <- Rtsne::Rtsne(x, dims, initial_dims, perplex, 
                           gravity, iter, check_duplicates = FALSE, 
                           verbose=TRUE, ...)
  #saved_tsne2 <- tsne.out
  tsne_plot <- cbind(meta, as.data.frame(tsne.out$Y))
  
  if(!is.character(do.plot)){
    return(list(raw = tsne.out, pca = tsne_plot))
  }
  
  #############################
  #plot the result
  cat("Preparing plot...\n")
  
  #Categories (pops, fathers, mothers, ect.) are given in do.plot argument. Supports up to two!
  #make the base plot, then add categories as color and fill.
  out <- ggplot(tsne_plot, aes(V1, V2)) + theme_bw() + 
    theme(axis.ticks = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(),
          panel.grid = element_blank()) #initialize plot
  
  
  long <- FALSE #are any of the sets of categories really long?
  
  #add variables.
  if(length(do.plot) == 1){
    v1 <- tsne_plot[,which(colnames(tsne_plot) == do.plot[1])] #get the factors
    v1u <- length(unique(v1)) #number of categories
    
    out <- out + geom_point(aes(color = v1)) #add the factor
    
    if(v1u >= 8){long <- TRUE} #are there too many categories to color with the cbb palette?
  }
  
  if(length(do.plot) == 2){
    v1 <- tsne_plot[,which(colnames(tsne_plot) == do.plot[1])]
    v2 <- tsne_plot[,which(colnames(tsne_plot) == do.plot[2])]
    v1u <- length(unique(v1))
    v2u <- length(unique(v2))
    
    out <- out + geom_point(aes(color = v1, fill = v2), pch = 21, size = 2.5, stroke = 1.25)
    
    if(v1u >= 8 | v2u >= 8){long <- TRUE}
  }
  
  #use the color blind friendly palette if possible!
  if(!long){
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    out <- out + scale_color_manual(values = cbbPalette)
  }
  
  return(list(tSNE = tsne.out, plot = out))
}

