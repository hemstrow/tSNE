summerize_genotypes <- function(x, ecs){
  meta <- x[,1:ecs]
  x <- x[,-c(1:ecs)]
  out <- data.frame(A = numeric(nrow(x)), T = numeric(nrow(x)), 
                    G = numeric(nrow(x)), C = numeric(nrow(x)))
  srow <- function(y){
    as <- c(substr(y, 1, 1), substr(y, 2, 2))
    as <- as[as != "N"]
    tab <- table(as)
    tab <- tab/sum(tab)
    return(tab)
  }
  for (i in 1:nrow(x)){
    tab <- srow(x[i,])
    if("A" %in% names(tab)){out[i,"A"] <- tab[names(tab) == "A"]}
    if("T" %in% names(tab)){out[i,"T"] <- tab[names(tab) == "T"]}
    if("G" %in% names(tab)){out[i,"G"] <- tab[names(tab) == "G"]}
    if("C" %in% names(tab)){out[i,"C"] <- tab[names(tab) == "C"]}
  }
  out <- cbind(meta, out)
  return(out)
}