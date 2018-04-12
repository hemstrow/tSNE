library(snpR); library(ggplot2); library(gridExtra);
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#Figures:
#1: SNP/Microsat pop structure and pedigrees in PCA and tSNE
#   a. Steelhead SNP pedigree
#   b. Chinook microsat pedigree
#   c. Monarch SNP pop structure
#   d. Steelhead microsat pop structure
#2: Stickleback
#   a: Whole genome
#   b: XIX
#   c: just inversion

#================Figure 1=============================

#############a#####
#SNPped
SNPped <- readRDS("../tSNE_data/steelhead/snps/RAPTURE/new/pa_genos.RDS")

SPP<- PCAfromPA(SNPped, 3, "fam")
SPP$plot <- SPP$plot + ggplot2::guides(color = FALSE)

SPT <- tSNEfromPA(SNPped, 3, "fam")
SPT$plot <- SPT$plot + ggplot2::guides(color = FALSE)

#############b#####
#Mped
Mped <- readRDS("../tSNE_data/chinook/pa11_genos.RDS")

MPP <- PCAfromPA(Mped, 4, c("mom", "dad"), TRUE)
MPP$plot <- MPP$plot + guides(color = FALSE, fill = FALSE)

MPT <- tSNEfromPA(Mped, 4, c("mom", "dad"), c.dup = TRUE, iter = 1000)
MPT$plot <- MPT$plot + guides(color = FALSE, fill = FALSE)

#############c####
SNPpop <- readRDS("../tSNE_data/monarch/pa_genotypes.RDS")

colnames(SNPpop)[2] <- "Population"

SPoP <- PCAfromPA(SNPpop, 2, "Population")
SPoP$plot
SPoT <- tSNEfromPA(SNPpop, 2, "Population")
SPoT$plot

#############d####
Mpop <- readRDS("../tSNE_data/steelhead/msats/pa_genos.RDS")

MPoP <- PCAfromPA(Mpop, 2, c.dup = TRUE)
MPoT <- tSNEfromPA(Mpop, 2, c.dup = TRUE)

#############combined####
#combine
grid.arrange(SPP$plot + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()), 
             SPT$plot, 
             MPP$plot + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()), 
             MPT$plot, 
             SPoP$plot + ggplot2::guides(color = FALSE) + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()),
             SPoT$plot + ggplot2::guides(color = FALSE),
             MPoP$plot + ggplot2::guides(color = FALSE) + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()),
             MPoT$plot + ggplot2::guides(color = FALSE),
             ncol = 2, layout_matrix = cbind(c(2,4,6,8), c(1,3,5,7)))

#why MPoP is so wierd: six very influential alleles! Maybe cutthroat associted or something?
ggplot(as.data.frame(MPoP$raw$rotation), aes(PC1, PC2, label = rownames(MPoP$raw$rotation))) + geom_text() + 
  theme_bw() + scale_x_continuous(limits = c(-.9, .9))


#================Figure S1==============
Polyped <- readRDS("../")

#================Figure 2===============
#reaq in data
Swhole <- readRDS("../tSNE_data/stickleback/pa_genos.RDS")
SXIX <- readRDS("../tSNE_data/stickleback/pa_XIX.RDS")
Sinv <- readRDS("../tSNE_data/stickleback/pa_IX.RDS")

#run tSNE at a range of perplexities.
##whole
colnames(Swhole)[2] <- "Population"
base_p <- mmtsne::hbeta(Swhole[,3:ncol(Swhole)], beta = 1)$H
swtSNElp <- tSNEfromPA(Swhole, 2, "Population", perplex = base_p/2)
swtSNEbp <- tSNEfromPA(Swhole, 2, "Population", perplex = base_p)
swtSNEhp <- tSNEfromPA(Swhole, 2, "Population", perplex = base_p*2)
swPCA <- PCAfromPA(Swhole, 2, "Population")
##XIX
colnames(SXIX)[2] <- "Population"
base_p <- mmtsne::hbeta(SXIX[,3:ncol(SXIX)], beta = 1)$H
sXIXtSNElp <- tSNEfromPA(SXIX, 2, "Population", perplex = base_p/2)
sXIXtSNEbp <- tSNEfromPA(SXIX, 2, "Population", perplex = base_p)
sXIXtSNEhp <- tSNEfromPA(SXIX, 2, "Population", perplex = base_p*2)
sXIXPCA <- PCAfromPA(SXIX, 2, "Population")
##only inversion
colnames(Sinv)[2:3] <- c("Population", "Haplotype")
base_p <- mmtsne::hbeta(Sinv[,4:ncol(Sinv)], beta = 1)$H
SinvtSNElp <- tSNEfromPA(Sinv, 3, "Population", perplex = base_p/2)
SinvtSNEbp <- tSNEfromPA(Sinv, 3, "Population", perplex = base_p)
SinvtSNEhp <- tSNEfromPA(Sinv, 3, "Population", perplex = base_p*2)
SinvPCA <- PCAfromPA(Sinv, 3, "Population")
###fix these plots
SinvtSNElp$plot <- ggplot(SinvtSNElp$plot$data, aes(V1, V2, color = Population, shape = Haplotype)) +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
  geom_point() + guides(color = FALSE, shape = FALSE) + scale_shape_manual(values = c(3,16,4)) +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"))

SinvtSNEbp$plot <- ggplot(SinvtSNEbp$plot$data, aes(V1, V2, color = Population, shape = Haplotype)) +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
  geom_point() + guides(color = FALSE, shape = FALSE) + scale_shape_manual(values = c(3,16,4)) +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"))

SinvtSNEhp$plot <- ggplot(SinvtSNEhp$plot$data, aes(V1, V2, color = Population, shape = Haplotype)) +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
  geom_point() + guides(color = FALSE, shape = FALSE) + scale_shape_manual(values = c(3,16,4)) +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"))

SinvPCA$plot <- ggplot(SinvPCA$plot$data, aes(PC1, PC2, color = Population, shape = Haplotype)) +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
  geom_point() + scale_shape_manual(values = c(3,16,4)) +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"))



#get legend
f2_legend <- get_legend(SinvPCA$plot)

#combine plots
grid.arrange(
  swtSNElp$plot+guides(color = FALSE),
  swtSNEbp$plot+guides(color = FALSE),
  swtSNEhp$plot+guides(color = FALSE),
  swPCA$plot+ theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())+guides(color = FALSE),
  sXIXtSNElp$plot+guides(color = FALSE),
  sXIXtSNEbp$plot+guides(color = FALSE),
  sXIXtSNEhp$plot+guides(color = FALSE),
  sXIXPCA$plot+ theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())+guides(color = FALSE),
  SinvtSNElp$plot,
  SinvtSNEbp$plot,
  SinvtSNEhp$plot,
  SinvPCA$plot + guides(color = FALSE, shape = FALSE),
  f2_legend,
  ncol = 5, 
  layout_matrix = cbind(c(1,5,9), c(2,6,10), c(3,7,11), c(4,8,12), c(13,13,13)),
  widths = c(1,1,1,1,0.5)
)

#================Figure S2==============
#cougar at multiple parameter sets
perps <- seq(.1, 1.9, by = 1.8/4) #which perplexities?
idims <- seq(.1, 1.9, by = 1.8/4) #which number of intial dimensions?

compDat <- Mped[1:500,]

base_hbeta <- mmtsne::hbeta(compDat[,5:ncol(compDat)])$H #base perplexity
base_dims <- length(unique(paste0(compDat$mom, compDat$dad))) #base number of families (dimensions)

#run each tSNE and save. Probably should have shunted this straight to a data frame, but oh well.
##run first
sth_parms <- tSNEfromPA(compDat, 4, c("mom", "dad"), initial_dims = floor(base_dims*idims[1]), perplex = base_hbeta*perps[1], gravity = 0.25)$plot$data
sth_parms <- cbind(perplexity = perps[1], initial_dims = idims[1], sth_parms)

for(i in 1:length(perps)){
  if(i == 1){j <- 2}
  else{j <- 1}
  for(j in j:length(idims)){
    out <- tSNEfromPA(compDat, 4, c("mom", "dad"), initial_dims = floor(base_dims*idims[j]), perplex = base_hbeta*perps[i], gravity = 0.25)$plot$data
    out <- cbind(perplexity = perps[i], initial_dims = idims[j], out)
    sth_parms <- rbind(sth_parms, out)
  }
}

#save this so we don't have to run it again
saveRDS(sth_parms, "../tSNE_data/chinook/multi_parameter_tSNE_figureS2.RDS")

ggplot(sth_parms, aes(V1, V2, color = mom, fill = dad)) + 
  facet_grid(perplexity ~ initial_dims, switch = "both") +
  geom_point(pch = 21, stroke = 1.25, size = 2.5) +
  theme_bw() +
  guides(color = FALSE, fill = FALSE) +
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        strip.text = element_text(size = 18),
        title = element_text(size = 22),
        axis.ticks = element_blank()) +
  ylab("Perplexity") + xlab("Initial Dimensions")

#looks great!
#================Figure S3===========
#IBS and post plots

#import data
IBS <- readRDS("../tSNE_data/steelhead/snps/RAPTURE/new/pa_IBS.RDS")
Posts <- readRDS("../tSNE_data/steelhead/snps/RAPTURE/new/pa_Posts.RDS")


#run PCA ant tSNE
idims <- length(unique(IBS$fam))

colnames(IBS)[3] <- "Family"
colnames(Posts)[3] <- "Family"

IBS.tSNE <- tSNEfromPA(IBS, 4, "Family")
IBS.PCA <- PCAfromPA(IBS, 4, "Family")

Posts.tSNE <- tSNEfromPA(Posts, 4, "Family")
Posts.PCA <- PCAfromPA(Posts, 4, "Family")

colnames(IBS.PCA$plot$data)[5:6] <- c("V1", "V2")
colnames(Posts.PCA$plot$data)[5:6] <- c("V1", "V2")

#combine and plot
IPcomb <- rbind(cbind(data = "IBS", method = "tSNE", IBS.tSNE$plot$data),
                cbind(data = "IBS", method = "PCA", IBS.PCA$plot$data[,1:6]),
                cbind(data = "Posterior", method = "tSNE", Posts.tSNE$plot$data),
                cbind(data = "Posterior", method = "PCA", Posts.PCA$plot$data[,1:6]))

ggplot(IPcomb, aes(V1, V2, fill = Family, color = miss)) +
  geom_point(pch = 21, size = 2.5, stroke = 1.5) +
  facet_wrap(data~method, scales = "free") + guides(fill = FALSE, color = guide_colorbar(title = "# Missing Genotypes")) +
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(size = 18),
        title = element_text(size = 22),
        axis.ticks = element_blank())
