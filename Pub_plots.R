library(snpR); library(ggplot2); library(gridExtra); library(grid); library(viridis)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
tb <- ttheme_minimal(core = list(fg_params = list(fontsize = 20)))
tbr <- ttheme_minimal(core=list(fg_params=list(rot=90, fontsize = 20)))
tr <- ttheme_minimal(core = list(fg_params = list(rot = 90)))

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
SNPped <- readRDS("../tSNE_data/chinook/snps/RAPTURE/new/pa_genos.RDS")

SPP<- PCAfromPA(SNPped, 3, "fam")
SPP$plot <- SPP$plot + ggplot2::guides(color = FALSE) + scale_color_viridis(discrete = T)

SPT <- tSNEfromPA(SNPped, 3, "fam")
SPT$plot <- SPT$plot + ggplot2::guides(color = FALSE) + scale_color_viridis(discrete = T)

#############b#####
#Mped
Mped <- readRDS("../tSNE_data/chinook/pa11_genos.RDS")

MPP <- PCAfromPA(Mped, 4, c("mom", "dad"), TRUE)
MPP$plot <- MPP$plot + guides(color = FALSE, fill = FALSE) + scale_color_viridis(discrete = T) + scale_fill_viridis(discrete = T)

MPT <- tSNEfromPA(Mped, 4, c("mom", "dad"), c.dup = TRUE, iter = 1000)
MPT$plot <- MPT$plot + guides(color = FALSE, fill = FALSE) + scale_color_viridis(discrete = T) + scale_fill_viridis(discrete = T)

#############c####
SNPpop <- readRDS("../tSNE_data/monarch/pa_genotypes.RDS")

colnames(SNPpop)[2] <- "Population"

SPoP <- PCAfromPA(SNPpop, 2, "Population")
SPoP$plot <- SPoP$plot + guides(color = FALSE) + scale_color_viridis(discrete = T)
SPoT <- tSNEfromPA(SNPpop, 2, "Population")
SPoT$plot <-SPoT$plot + guides(color = FALSE) + scale_color_viridis(discrete = T)

#############d####
Mpop <- readRDS("../tSNE_data/steelhead/msats/pa_genos.RDS")

MPoP <- PCAfromPA(Mpop, 2, c.dup = TRUE)
MPoP$plot <- MPoP$plot + guides(color = FALSE)
MPoT <- tSNEfromPA(Mpop, 2, c.dup = TRUE)
MPoT$plot <- MPoT$plot + guides(color = FALSE)

#############combined####
#combine

p1r1 <- grid.arrange(tableGrob("SNP", theme = tr),
                     SPP$plot + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()), 
                     SPT$plot,
                     nrow = 1, widths = c(0.05, 1, 1))
p1r2 <- grid.arrange(tableGrob("Microsatellite", theme = tr),
                     MPP$plot + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()), 
                     MPT$plot,
                     nrow = 1, widths = c(0.05, 1, 1))
p1r3 <- grid.arrange(tableGrob("SNP", theme = tr),
                     SPoP$plot + ggplot2::theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()),
                     SPoT$plot,
                     nrow = 1, widths = c(0.05, 1, 1))
p1r4 <- grid.arrange(tableGrob("Microsatellite", theme = tr),
                     MPoP$plot + ggplot2::theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()), 
                     MPoT$plot,
                     nrow = 1, widths = c(0.05, 1, 1))

p1t <- grid.arrange(p1r1, p1r2)
p1t <- grid.arrange(tableGrob("Pedigree", theme = tbr), p1t, widths = c(0.05, 2.05))
p1b <- grid.arrange(p1r3, p1r4)
p1b <- grid.arrange(tableGrob("Population Structure", theme = tbr), p1b, widths = c(0.05, 2.05))

p1x <- grid.arrange(tableGrob("", theme = ttheme_minimal()),
                    tableGrob("", theme = ttheme_minimal()),
                    tableGrob("PCA", theme = ttheme_minimal()), 
                    tableGrob("t-SNE", theme = ttheme_minimal()), nrow = 1,
                    widths = c(0.05, 0.05, 1, 1))
p1xb <- grid.arrange(tableGrob("", theme = ttheme_minimal()),
                     tableGrob("", theme = ttheme_minimal()),
                     tableGrob("Method", theme = tb),
                     nrow = 1, widths = c(0.05, 0.05, 2))

p1c <- grid.arrange(p1t, p1b, p1x, p1xb,
                    nrow = 4, heights = c(1,1,0.07, 0.07))

#why MPoP is so wierd: six very influential alleles! Maybe cutthroat associted or something?
ggplot(as.data.frame(MPoP$raw$rotation), aes(PC1, PC2, label = rownames(MPoP$raw$rotation))) + geom_text() + 
  theme_bw() + scale_x_continuous(limits = c(-.9, .9))


#================Figure S1==============
#pop data
sturpop <- readRDS("../tSNE_data/sturgeon/formatted_frasier_pa.RDS")

stuPotsne <- tSNEfromPA(Polyped, 2, "Population")
stuPoPCA <- PCAfromPA(sturpop, 2, "Population")

#remove guides
stuPoPCA$plot <- stuPoPCA$plot + guides(color = FALSE)
stuPotsne$plot <- stuPotsne$plot + guides(color = FALSE)


#pedigree data
sturped <- readRDS("../tSNE_data/sturgeon/sturg_ped_pa.RDS")
stuPetsne <- tSNEfromPA(sturped, 2, c("dad", "mom"))
stuPePCA <- PCAfromPA(sturped, 2, c("dad", "mom"))

#fix color scheme in ped
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
stuPePCA$plot <- stuPePCA$plot + scale_fill_manual(values = cbbPalette) + scale_color_viridis(discrete = T) + guides(color = FALSE, fill = FALSE)
stuPetsne$plot <- stuPetsne$plot + scale_fill_manual(values = cbbPalette) + scale_color_viridis(discrete = T) + guides(color = FALSE, fill = FALSE)

#remove axis stuff from PCAs
stuPePCA$plot <- sturPePCA$plot + theme(axis.title = element_blank(), 
                                        axis.text = element_blank(), 
                                        axis.ticks = element_blank(), 
                                        panel.grid = element_blank())


stuPoPCA$plot <- stuPoPCA$plot + theme(axis.title = element_blank(), 
                                        axis.text = element_blank(), 
                                        axis.ticks = element_blank(), 
                                        panel.grid = element_blank())


#combine
fS1plots <- grid.arrange(stuPoPCA$plot, stuPotsne$plot, stuPePCA$plot, stuPetsne$plot,
                        layout_matrix = cbind(c(1,3), c(2,4)))

fS1x <- grid.arrange(tableGrob("", theme = ttheme_minimal()),
                    tableGrob("PCA", theme = ttheme_minimal()), 
                    tableGrob("t-SNE", theme = ttheme_minimal()), nrow = 1,
                    widths = c(0.1, 1, 1))

fS1xb <- grid.arrange(tableGrob("", theme = ttheme_minimal()),
                     tableGrob("Method", theme = tb),
                     nrow = 1, widths = c(0.1, 2))

fS1y <- grid.arrange(tableGrob("Population Structure", theme = tbr),
                     tableGrob("Pedigree", theme = tbr),
                     ncol = 1)

fS1 <- grid.arrange(fS1y, fS1plots, nrow = 1, widths = c(0.1, 2))
fS1 <- grid.arrange(fS1, fS1x, fS1xb, ncol = 1, heights = c(2, 0.075, 0.1))

#================Figure 2===============
#reaq in data
Swhole <- readRDS("../tSNE_data/stickleback/pa_genos.RDS")
SXIX <- readRDS("../tSNE_data/stickleback/pa_nXIX.RDS")
Sinv <- readRDS("../tSNE_data/stickleback/pa_IX.RDS")

#run tSNE at a range of perplexities.
##whole
colnames(Swhole)[2] <- "Population"
swtSNElp <- tSNEfromPA(Swhole, 2, "Population", perplex = 2)
swtSNEbp <- tSNEfromPA(Swhole, 2, "Population", perplex = 15)
swtSNEhp <- tSNEfromPA(Swhole, 2, "Population", perplex = 50)
swPCA <- PCAfromPA(Swhole, 2, "Population")
##no XIX
colnames(SXIX)[2] <- "Population"
sXIXtSNElp <- tSNEfromPA(SXIX, 2, "Population", perplex = 2)
sXIXtSNEbp <- tSNEfromPA(SXIX, 2, "Population", perplex = 15)
sXIXtSNEhp <- tSNEfromPA(SXIX, 2, "Population", perplex = 50)
sXIXPCA <- PCAfromPA(SXIX, 2, "Population")
##only inversion
colnames(Sinv)[2:3] <- c("Population", "Haplotype")
SinvtSNElp <- tSNEfromPA(Sinv, 3, "Population", perplex = 2)
SinvtSNEbp <- tSNEfromPA(Sinv, 3, "Population", perplex = 15)
SinvtSNEhp <- tSNEfromPA(Sinv, 3, "Population", perplex = 50)
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


F2xh <- grid.arrange(tableGrob("2", theme = ttheme_minimal()), 
                    tableGrob("15", theme = ttheme_minimal()),
                    tableGrob("50", theme = ttheme_minimal()),
                    tableGrob("", theme = ttheme_minimal()),
                    tableGrob("NA", theme = ttheme_minimal()), nrow = 1, widths = c(1,1,1,0.1,1),
                    left = " ")
F2xl <- grid.arrange(tableGrob("t-SNE", theme = tb),
                     tableGrob("", theme = tb),
                     tableGrob("PCA", theme = tb),
                     nrow = 1, widths = c(3, .1, 1),
                     left = " ")
F2r1 <- grid.arrange(swtSNElp$plot+guides(color = FALSE),
                     swtSNEbp$plot+guides(color = FALSE),
                     swtSNEhp$plot+guides(color = FALSE),
                     tableGrob("", theme = ttheme_minimal()),
                     swPCA$plot+ theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) + guides(color = FALSE),
                     nrow = 1, widths = c(1,1,1,0.1,1),
                     left = "All chromosomes")
F2r2 <- grid.arrange(sXIXtSNElp$plot+guides(color = FALSE),
                     sXIXtSNEbp$plot+guides(color = FALSE),
                     sXIXtSNEhp$plot+guides(color = FALSE),
                     tableGrob("", theme = ttheme_minimal()),
                     sXIXPCA$plot+ theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())+guides(color = FALSE),
                     nrow = 1, widths = c(1,1,1,0.1,1),
                     left = "No sex chromosome")
F2r3 <- grid.arrange(SinvtSNElp$plot,
                     SinvtSNEbp$plot,
                     SinvtSNEhp$plot,
                     tableGrob("", theme = ttheme_minimal()),
                     SinvPCA$plot + guides(color = FALSE, shape = FALSE),
                     nrow = 1, widths = c(1,1,1,0.1,1),
                     left = "Inversion")


F2 <- grid.arrange(F2r1, F2r2, F2r3, F2xh, F2xl, ncol = 1, heights = c(1,1,1,.15,.16))

F2 <- grid.arrange(F2, f2_legend, ncol = 2, widths = c(1, 0.15))


#================Figure S2==============
#cougar at multiple parameter sets
perps <- seq(2, 50, length = 4) #which perplexities?
thetas <- seq(0, 1, length = 4) #which thetas?

compDat <- Mped[1:500,]


#run each tSNE and save. Probably should have shunted this straight to a data frame, but oh well.
##run first
sth_parms <- tSNEfromPA(compDat, 4, c("mom", "dad"), perplex = perps[1], gravity = thetas[1])$plot$data
sth_parms <- cbind(perplexity = perps[1], theta = thetas[1], sth_parms)

for(i in 1:length(perps)){
  if(i == 1){j <- 2}
  else{j <- 1}
  for(j in j:length(thetas)){
    out <- tSNEfromPA(compDat, 4, c("mom", "dad"), perplex = perps[i], gravity = thetas[i])$plot$data
    out <- cbind(perplexity = perps[i], theta = thetas[j], out)
    sth_parms <- rbind(sth_parms, out)
  }
}

#save this so we don't have to run it again
saveRDS(sth_parms, "../tSNE_data/chinook/multi_parameter_tSNE_figureS2.RDS")
sth_parms <- readRDS("../tSNE_data/chinook/multi_parameter_tSNE_figureS2.RDS")

FS2 <- ggplot(sth_parms, aes(V1, V2, color = mom, fill = dad)) + 
  facet_wrap(perplexity ~ theta, scales = "free") +
  geom_point(pch = 21, stroke = 1.25, size = 2.5) +
  theme_bw() +
  guides(color = FALSE, fill = FALSE) +
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        strip.text = element_blank(),
        title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  scale_color_viridis(discrete = T)

FS2xh <- grid.arrange(tableGrob("", theme = ttheme_minimal()),
                     tableGrob("", theme = ttheme_minimal()),
                     tableGrob("0", theme = ttheme_minimal()), 
                     tableGrob("1/3", theme = ttheme_minimal()),
                     tableGrob("2/3", theme = ttheme_minimal()),
                     tableGrob("1", theme = ttheme_minimal()), 
                     widths = c(0.075, 0.075, 1, 1, 1, 1), nrow = 1)

FS2xl <- grid.arrange(tableGrob("", theme = ttheme_minimal()),
                     tableGrob("", theme = ttheme_minimal()),
                     tableGrob("Theta", theme = tb), nrow = 1,
                     widths = c(0.075, 0.075, 4))

FS2y <- grid.arrange(tableGrob("2", theme = ttheme_minimal()), 
                     tableGrob("18", theme = ttheme_minimal()),
                     tableGrob("34", theme = ttheme_minimal()),
                     tableGrob("50", theme = ttheme_minimal()), ncol = 1)

FS2c <- grid.arrange(tableGrob("Perplexity", theme = tbr), FS2y, FS2, nrow = 1, widths = c(0.075, 0.075, 4))
FS2c <- grid.arrange(FS2c, FS2xh, FS2xl, ncol = 1, heights = c(1, 0.05, 0.05))


#looks great!
#================Figure S3===========
#IBS and post plots

#import data
IBS <- readRDS("../tSNE_data/chinook/snps/RAPTURE/new/pa_IBS.RDS")
Posts <- readRDS("../tSNE_data/chinook/snps/RAPTURE/new/pa_Posts.RDS")


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

f3sl <- get_legend(ggplot(IPcomb, aes(V1, V2, fill = Family, color = miss)) +
  geom_point(pch = 21, size = 2.5, stroke = 1.5) +
  facet_wrap(data~method, scales = "free") +
  theme_bw() +
  guides(fill = FALSE, color = guide_colorbar(title = "# Missing Genotypes", title.theme = element_text(size = 12, angle = 0))) +
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text = element_blank(),
        title = element_text(size = 22),
        axis.ticks = element_blank()))

f3sm <- ggplot(IPcomb, aes(V1, V2, fill = Family, color = miss)) +
  geom_point(pch = 21, size = 2.5, stroke = 1.5) +
  facet_wrap(data~method, scales = "free") +
  theme_bw() +
  guides(fill = FALSE, color = FALSE) +
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text = element_blank(),
        title = element_text(size = 22),
        axis.ticks = element_blank())

f3sx <- grid.arrange(tableGrob("", theme = tb),
                     tableGrob("t-SNE", theme = tb), 
                     tableGrob("PCA", theme = tb), nrow = 1,
                     widths = c(0.1, 1, 1))

f3sy <- grid.arrange(tableGrob("IBS", theme = tbr), 
                     tableGrob("Posteriors", theme = tbr), ncol = 1,
                     heights = c(1, 1))

f3st <- grid.arrange(f3sy, f3sm, nrow = 1, widths = c(0.1, 2))
f3s <- grid.arrange(f3st, f3sx, ncol = 1, heights = c(2, 0.1))
f3s <- grid.arrange(f3s, f3sl, nrow = 1, widths = c(2.1, .5))

#================Figure S4=================
#bootstrap the inversion data a few times and plot
colnames(Sinv)[2:3] <- c("Population", "Haplotype")
bdat <- Sinv[,4:ncol(Sinv)]
boot1 <- cbind(Sinv[,1:3], bdat[,sample(168, 186, TRUE)])
btsne1 <- tSNEfromPA(boot1, 3, "Population")

boot2 <- cbind(Sinv[,1:3], bdat[,sample(168, 186, TRUE)])
btsne2 <- tSNEfromPA(boot2, 3, "Population")

boot3 <- cbind(Sinv[,1:3], bdat[,sample(168, 186, TRUE)])
btsne3 <- tSNEfromPA(boot3, 3, "Population")

boot4 <- cbind(Sinv[,1:3], bdat[,sample(168, 186, TRUE)])
btsne4 <- tSNEfromPA(boot4, 3, "Population")

bootc <- rbind(cbind(btsne1$plot$data, boot = 1),
               cbind(btsne2$plot$data, boot = 2),
               cbind(btsne3$plot$data, boot = 3),
               cbind(btsne4$plot$data, boot = 4))

f4sl <- get_legend(ggplot(bootc, aes(V1, V2, color = Population, shape = Haplotype)) + facet_wrap(~boot) +
  theme_bw() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  geom_point() + scale_shape_manual(values = c(3,16,4)) +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")))

f4so <- tSNEfromPA(Sinv, 3, "Population")
f4so<- ggplot(f4so$plot$data, aes(V1, V2, color = Population, shape = Haplotype)) +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
  geom_point() + guides(color = FALSE, shape = FALSE) + scale_shape_manual(values = c(3,16,4)) +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
  ggtitle("Data:")

f4sb <- ggplot(bootc, aes(V1, V2, color = Population, shape = Haplotype)) + facet_wrap(~boot) +
  theme_bw() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  geom_point() + scale_shape_manual(values = c(3,16,4)) +
  guides(shape = FALSE, color = FALSE) +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) +
  ggtitle("Bootstraps:")

f4S <- grid.arrange(f4so, f4sb, f4sl, layout_matrix = cbind(c(1,2), c(3,3)), widths = c(1,0.2))


#================Table S1============
#table showing what allele freqnecy differences the tSNE run is working with for the stickle inversion.

#get clusters, note that this code needs to be adjusted every time!
r <- tSNEfromPA(Sinv, 3, c("Population", "Haplotype"))
dat2 <- r$plot$data

dat2$clust <- NA

dat2$clust[dat2$V1 > 15] <- 1
dat2$clust[dat2$V1 < -30] <- 2
dat2$clust[dat2$V2 > 0 & dat2$V1 < 30 & dat2$V1 > -30] <- 3
dat2$clust[dat2$V2 < 0 & dat2$V1 > -40 & dat2$V1 < 20] <- 4

ggplot(dat2, aes(V1, V2, color = as.factor(clust))) + geom_point() #check

#add genotypes.
sdat <- read.table("../tSNE_data/stickleback/snps_numeric_filt.txt", colClasses = "character", header = T)
sdat$position <- as.numeric(sdat$position)

colnames(sdat)[4:ncol(sdat)] <- paste0(dat2$clust, "_", dat2$samp)

meta <- sdat[,1:3]
sdat <- sdat[,4:ncol(sdat)]
sdat <- sdat[,order(colnames(sdat))]
sdat <- cbind(meta, sdat)

cl <- table(substr(colnames(sdat)[4:ncol(sdat)], 1, 1))
cl <- list(paste0("Cluster_", names(cl)), as.numeric(cl))

#subset inversion snps
rFST <- read.table("~/Stickleback/Full Data/2017_reruns/PlotData/FST/rawfst.txt", header = T)
isnps <- rFST[rFST$group == "groupIX" & rFST$position >= 13*1000000& rFST$position <= 18.5*1000000 & rFST$comp == "ASP_OPL",]
isnps <- isnps[isnps$Fst >= 0.15, 1]

sdat <- sdat[sdat$snp %in% isnps,]


#ac format
ac_sdat <- format_snps(sdat, 3, "ac", input_form = "0000", miss = "00", pop = cl)
ac_sdat$fp <- ac_sdat$ni1/ac_sdat$n_total

#get allele frequency info, subset out anything that is just different on the inversion
mac <- reshape2::dcast(ac_sdat[,c(3,4,9)], pop ~ position)
mac <- t(mac)
mac <- mac[-1,]
mac <- as.numeric(mac)
mac <- matrix(mac, ncol = 4)
colnames(mac) <- paste0("Cluster ", 1:4)
rownames(mac) <- ac_sdat[ac_sdat$pop == "Cluster_1",]$position
mac <- mac[-which(rowSums(ifelse(mac == 1, 1, 0)) == 3),]
mac <- mac[-which(rowSums(ifelse(mac == 0, 1, 0)) == 3),]

#polish
mac <- as.data.frame(mac)
mac <- cbind('position(mb)' = as.numeric(rownames(mac))/1000000, mac)
mac[,2:5] <- round(mac[,2:5], 2)
