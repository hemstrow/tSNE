#This script prepares presence absence data from the microsatellite and SNP data as described in
#Thorstensen and Hemstrom 2018. Note that it will not run unless the raw data and metadata for each set is provided!


#order:
# 1. Monarch SNP data
# 2. Steelhead microsats.
# 3. Steelhead SNP data.
# 4. Chinook microsats
# 5. Stickle SNP data
# 6. Sturgeon pop structure
# 7. Sturgeon family data
# 8. Steelhead IBS and Posteriors.



library(snpR); library(ggplot2); library(readr)

#===================================================Monarch SNP data=================

#import data

setwd("~/2017-2018/Monarchs/")
#import bamlist and annotations
bamlist <- read.table("gbamlist.txt")
plate1 <- read.csv("plate1.csv", header = T)
plate2 <- read.csv("plate2_MFMB1.csv", header = T, skip = 1)
plate3 <- read.csv("plate3_MFMB2.csv", header = T, skip = 1)

DNA274 <- "ACAGTG"
DNA275 <- "CAGATC"

plate2 <- rbind(plate2[,1:3], setNames(plate2[,4:6], colnames(plate2)[1:3]))
plate3 <- rbind(plate3[,1:3], setNames(plate3[,4:6], colnames(plate3)[1:3]))


#add index
barcodes <- read.table("RAD barcodes.txt", header = T)
plate2 <- merge(plate2, barcodes, by = "Well")
plate3 <- merge(plate3, barcodes, by = "Well")
plate1 <- merge(plate1, barcodes, by = "Well")


#get sample info
sampinfo <- data.frame(plate = substr(bamlist[,1], 8, 13), ind = substr(bamlist[,1], 20, 27), stringsAsFactors = FALSE)
sampinfo$plate[sampinfo$plate == "_split"] <- "plate1"



#add sample info to plates
plate2 <- merge(sampinfo[sampinfo$plate == DNA274,], plate2, by.x = "ind", by.y = "Index")
plate3 <- merge(sampinfo[sampinfo$plate == DNA275,], plate3, by.x = "ind", by.y = "Index")
plate1 <- merge(sampinfo[sampinfo$plate == "plate1",], plate1, by.x = "ind", by.y = "Index")

#combine plates in the correct order
combplates <- rbind(plate2, plate3, plate1)

#fix up some pop IDs
pops <- substr(combplates$Pop, 1, 3)
table(pops)
pops[pops == "Gua"] <- "GUA"
pops[pops == "New"] <- "NCA"
pops[pops == "Nor"] <- "NOR"
pops[pops == "NZ2"] <- "NZL"
pops[pops == "NZR"] <- "NZL"
pops[pops == "Rot"] <- "ROT"
pops[pops == "Sai"] <- "SAI"
pops[pops == "Sam"] <- "SAM"
table(pops)
combplates$Pop <- pops


#import raw genotypes.
raw_genos <- read.table("genotypes.geno", header = F, stringsAsFactors = F)

#associate sample IDs to genotypes. Works as long as the combplates are in the same order as the bamfile.
colnames(raw_genos) <- c("group", "position", paste0(combplates$Pop, "_", combplates$ID))

#add snp numbers
raw_genos <- cbind(snp = 1:nrow(raw_genos), raw_genos)

#sort by sample ID (and thus by population)
temp <- raw_genos[,4:ncol(raw_genos)]
raw_genos <- cbind(raw_genos[,1:3], temp[,order(colnames(temp))])

#grab pop info
pops <- substr(colnames(raw_genos)[4:ncol(raw_genos)], 1, 3)
pops <- table(pops)
pops <- list(names(pops), as.numeric(pops))

#Filter and prepare data
flt_genos <- filter_snps(raw_genos, 3, 0.05, 0.55, floor((ncol(raw_genos)-3)/2), 0.5, pop = pops) #just filtering snps.
pa_genos <- format_snps(flt_genos, 3, 7)
pa_genos$samp <- as.character(pa_genos$samp)
pops <- substr(pa_genos$samp, 1, 3)
pa_genos <- cbind(pa_genos[,1], Population = pops, pa_genos[,2:ncol(pa_genos)])

saveRDS(pa_genos, "~/GitHub/tSNE_data/monarch/pa_genotypes.RDS")


#===================================================Steelhead microsats=============
rm(list = ls())

#get data and metadata
setwd("~/GitHub/tSNE_data/steelhead/msats/")
genos <- read.table("steelhead_genotypes.txt", sep='', header=FALSE, colClasses = "character")
genos <- t(genos)
pops <- genos[1,]
genos <- genos[-1,]
colnames(genos) <- pops
lnam <- c("Oki23", "Ssa407", "mSsa408", "Ots209", "OtsG249b", "OtsG85", "Omy27", "Omy1001", "Ots243", "Ots409", "OtsG3", "Ots212", "Omm1087")
genos <- cbind(locus = lnam, lnum = 1:13, genos)

#any to filter?
mloc <- rowSums(ifelse(genos[,3:ncol(genos)] == "0000", 1, 0))
any(mloc > 0.5*(ncol(genos) - 2)) #no loci to remove
mind <- colSums(ifelse(genos[,3:ncol(genos)] == "0000", 1, 0))
any(mind > 0.5*(nrow(genos))) #no individuals to remove

#reformat and save
pa_genos <- format_snps(genos, 2, 7, "msat_2", "00", lnames = lnam)
pa_genos <- cbind(samp = 1:nrow(pa_genos), pop = pa_genos$samp, pa_genos[,-1])

saveRDS(pa_genos, "pa_genos.RDS")

#===================================================Steelhead SNP data===============
rm(list = ls())

#metadata:
bamlist <- read.table("../bamlist.txt")
bamlist <- as.vector(t(bamlist))

#get ind IDs
samps <- gsub(".+bam/", "", bamlist)
samps <- gsub("_.+$", "", samps)

#get pops
pops <- gsub(".+[0-9]{4}_", "", bamlist)
pops <- gsub("_.+", "", pops)

#get families
fams <- gsub(".+_", "", bamlist)
fams <- gsub("\\.bam", "", fams)




#get data
setwd("~/GitHub/tSNE_data/steelhead/snps/RAPTURE/new/")
genos <- read.table("rapGenos_clean.geno")
colnames(genos) <- c("group", "position", paste0(samps, "_", pops, "_", fams))
flt_genos <- filter_snps(genos, 2, 0.05, 0.55, floor((ncol(genos)-2)/2), .5)
pa_genos <- format_snps(flt_genos, 2, 7)
pa_genos$samp <- as.character(pa_genos$samp)
meta <- matrix(unlist(strsplit(pa_genos$samp, "_")), nrow(pa_genos), 3, T)
pa_genos <- cbind(data.frame(pop = meta[,2], fam = meta[,3], stringsAsFactors = F), pa_genos, stringsAsFactors = F)

saveRDS(pa_genos, "pa_genos.RDS")

#===================================================Chinook microsats================
rm(list = ls())

setwd("~/GitHub/tSNE_data/chinook/")
#prepare data

sample_table <- read.table("14JUL2015_cougar_juvenile_master.txt",sep='', header=TRUE, colClasses = "character", stringsAsFactors = F)
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

# Throwing loci and sample names in
colnames(format_input) <- samps
format_input <- cbind(loci = 1:nrow(format_input), marker = lnam, format_input)

#any to filter?
mloc <- rowSums(ifelse(format_input[,3:ncol(format_input)] == "000000", 1, 0))
any(mloc > 0.5*(ncol(format_input) - 2)) #no loci to remove
mind <- colSums(ifelse(format_input[,3:ncol(format_input)] == "000000", 1, 0))
any(mind > 0.5*(nrow(format_input))) #no individuals to remove


# Will's reformatting function
pa_genos <- format_snps(format_input, 2, output = 7,
                               input_form = "msat_3", miss = "000", lnames = lnam)

# Creating a new table of those 2011 juvs from the master table
single_year_class <- pa_genos[grep("MR11",pa_genos[,1]),]


#merge with sib data.
sibdata <- read.table("aj.pedigree.10.mendel.checked.txt", sep = "\t", header=T)
colnames(single_year_class)[1] <- "off"
comb <- merge(x = sibdata, single_year_class, by = "off")

#save
saveRDS(comb, "pa11_genos.RDS")



#===================================================Stickleback SNPs====================
#import data and make pop info list
genos <- read.table("~/Stickleback/Full Data/2017_reruns/Rewrites/github_repo/Hemstrom_et_al_2018/data/snps_numeric_filt.txt", header = T, colClasses = "character")
pops <- substr(colnames(genos)[-c(1:3)], 1, 3)
pops <- table(pops)
pops <- list(names(pops), as.numeric(pops))

#filter data
# genos <- filter_snps(genos, 3, 0.05, 0.55, floor((ncol(genos)-3)/2), 0.5, pop = pops, mDat = "0000")
#note, this data has already been filtered, at similar parameters, so no change)

#format data, multiple forms.
#full
pa_genos <- format_snps(genos, 3, 7, "0000", "00")
pa_genos$pop <- gsub("\\.\\w+", "", pa_genos$samp)
pa_genos <- pa_genos[,c(1, ncol(pa_genos), 2:ncol(pa_genos))]
pa_genos <- pa_genos[,-ncol(pa_genos)]

#noXIX
pa_nXIX <- format_snps(genos[genos$group != "groupXIX",], 3, 7, "0000", "00")
pa_nXIX$pop <- gsub("\\.\\w+", "", pa_nXIX$samp)
pa_nXIX <- pa_nXIX[,c(1, ncol(pa_nXIX), 2:ncol(pa_nXIX))]
pa_nXIX <- pa_nXIX[,-ncol(pa_nXIX)]

#XIX only
pa_XIX <- format_snps(genos[genos$group == "groupXIX",], 3, 7, "0000", "00")
pa_XIX$pop <- gsub("\\.\\w+", "", pa_XIX$samp)
pa_XIX <- pa_XIX[,c(1, ncol(pa_XIX), 2:ncol(pa_XIX))]
pa_XIX <- pa_XIX[,-ncol(pa_XIX)]

#inversion only
rFST <- read.table("~/Stickleback/Full Data/2017_reruns/PlotData/FST/rawfst.txt", header = T)
isnps <- rFST[rFST$group == "groupIX" & rFST$position >= 13*1000000& rFST$position <= 18.5*1000000 & rFST$comp == "ASP_OPL",]
isnps <- isnps[isnps$Fst >= 0.15, 1]
pa_IX <- format_snps(genos[genos$snp %in% isnps,], 3, 7, "0000", "00")
pa_IX$pop <- gsub("\\.\\w+", "", pa_IX$samp)
pa_IX <- pa_IX[,c(1, ncol(pa_IX), 2:ncol(pa_IX))]
pa_IX <- pa_IX[,-ncol(pa_IX)]


##get inversion state from input data.
igenos <- genos[genos$snp %in% isnps,]
ac_igenos <- format_snps(igenos, 3, 1, "0000", miss = "00", pop = pops)

###find common allele in ASP (will be the non-inversion allele)
pstate <- character(nrow(igenos))
agenos <- igenos[,grepl("ASP", colnames(igenos))]
for(i in 1:length(pstate)){
  v <- c(substr(agenos[i,], 1, 2), substr(agenos[i,], 3, 4))
  v <- table(v)
  v <- v[names(v) != "00"]
  pstate[i] <- names(v)[which.max(v)]
}

###determine inversion haplotype for each individual.
haps <- character(ncol(igenos) - 3)
as <- numeric(nrow(igenos))
midbounds <- c(0.4, 0.6) #set lower and upper boudary common allele frequencies for heterozygotes
ubounds <- c(0.1, 0.9) #set common allele frequencies for homozygous q and p respectively
for(i in 1:length(haps)){
  #get the number of common alleles at each locus in this ind.
  as <- substr(igenos[,i+3], 1, 2) == pstate #is allele one p?
  as <- as + (substr(igenos[,i+3], 3, 4) == pstate) #is allele two p?
  if(any(igenos[,i+3] == "0000")){
    as <- as[-which(igenos[,i+3] == "0000")] #remove missing genotypes
  }
  
  #call genotype
  as <- mean(as)
  as <- as/2
  
  
  if(length(as) == 0){
    haps[i] <- "uncalled"
  }
  if(as <= ubounds[1]){
    haps[i] <- "homozygous q"
  }
  else if (as >= ubounds[2]){
    haps[i] <- "homozygous p"
  }
  else if (as >= midbounds[1] & as <= midbounds[2]){
    haps[i] <- "heterozygous"
  }
  else{
    haps[i] <- "uncalled"
  }
}

pa_IX <- cbind(pa_IX[,1:2], haplotype = haps, pa_IX[,3:ncol(pa_IX)])

saveRDS(pa_genos, "../tSNE_data/stickleback/pa_genos.RDS")
saveRDS(pa_nXIX, "../tSNE_data/stickleback/pa_nXIX.RDS")
saveRDS(pa_XIX, "../tSNE_data/stickleback/pa_XIX.RDS")
saveRDS(pa_IX, "../tSNE_data/stickleback/pa_IX.RDS")

#===================================================Sturgeon pop structure============
Fr.genos <- read_delim("../tSNE_data/sturgeon/Fraser River data for t-SNE test.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#format into pa

##need to make a nested list containing locus information
lnames <- unlist(strsplit(colnames(Fr.genos)[2:ncol(Fr.genos)], "_"))
lnames <- lnames[nchar(lnames) > 3]
lnames <- unique(lnames)
lopts <- list()
for(i in 1:length(lnames)){
  tcols <- grep(paste0(lnames[i], "_"), colnames(Fr.genos))
  tcols <- sort(c(which(colnames(Fr.genos) == lnames[i]),tcols))
  topts <- as.numeric(Fr.genos[1,tcols])
  lopts[[i]] <- list(locus = lnames[i], alleles = topts, cols = tcols)
}

#function to get this data into presence absence. Returns this plus a count of missing data per ind
#Arguments:
# x: input data. make sure that the FIRST row is the first row of DATA (not header info)!
#    Note that the program expects absences to be noted by NA. Anything at all in the cell will be noted as a presence.
# ecs: number of header COLUMNS before the first column of data.
# linfo: locus info. A nested list. At the outer level, each element is a locus. Each inner level list
#        contains locus information, in a list as given by:
#        list(locus = "locus name, alleles = c("a1", "a2", "a3", ...), cols = numeric vector of the columns containing the allele presence absences for the locus))
format_mp_pa <- function(x, ecs, linfo){
  #grab header info and convert the NAs to zeros and the others to 1
  x[,(ecs+1):ncol(x)] <- ifelse(is.na(x[,(ecs+1):ncol(x)]), 0, 1)

  #initialize missing count
  miss <- numeric(length(nrow(x)))

  #get observed allele frequencies
  for(i in 1:length(linfo)){
    cat("Formating locus:", linfo[[i]]$locus, "\n")
    xi <- x[,linfo[[i]]$cols]
    xm <- colMeans(xi)
    miss <- miss + (rowSums(xi) == 0)
    xi[rowSums(xi) == 0,] <- xm
    x[,linfo[[i]]$cols] <- xi
    colnames(x)[linfo[[i]]$cols] <- paste0(linfo[[i]]$locus, "_", linfo[[i]]$alleles)
  }
  return(list(pam = x, miss = miss))
}

##reformat and filter poorly sequenced individuals. No crappy loci in this data set, they were already removed.
pa.geno <- format_mp_pa(Fr.genos[2:nrow(Fr.genos),], 1, lopts)
pa.geno <- pa.geno$pam[-which(pa.geno$miss >= floor((length(lopts))/2)),]
pa.geno <- as.data.frame(pa.geno)

##add pop info
pa.geno <- cbind(samp = pa.geno[,1], Population = gsub("[0-9]+", "", pa.geno$X1), pa.geno[,2:ncol(pa.geno)])

#save
saveRDS(pa.geno, "../tSNE_data/sturgeon/formatted_frasier_pa.RDS")



#===================================================Chinook IBS and Posteriors============
setwd("../tSNE_data/steelhead/snps/RAPTURE/new/")

#import metadata
bamlist <- read.table("../bamlist.txt")
bamlist <- as.vector(t(bamlist))

##get ind IDs
samps <- gsub(".+bam/", "", bamlist)
samps <- gsub("_.+$", "", samps)

##get pops
pops <- gsub(".+[0-9]{4}_", "", bamlist)
pops <- gsub("_.+", "", pops)

##get families
fams <- gsub(".+_", "", bamlist)
fams <- gsub("\\.bam", "", fams)



#import IBS and process
IBS <- read.table("rapIBS_clean.ibs", header = T)

##correctly fromat and filter IBS

##function to do this. Arguments:
##  IBS:input IBS
##  lcut1: missing data locus filter
##  icut: missing data individual fitlter
##  meta: metadata for samples
##  missing: identifier for missing data
prepIBS <- function(IBS, lcut1, icut, meta, missing = -1){
  #filter poor loci
  IBS[IBS == missing] <- NA
  IBS <- IBS[-which(rowSums(is.na(IBS)) >= (ncol(IBS))*(1-lcut1)),]
  
  #remove poorly sequenced individuals
  IBS <- t(IBS)
  IBS <- cbind(meta, IBS)
  IBS <- IBS[-which(rowSums(is.na(IBS)) >= (ncol(IBS) - ncol(meta))*(1-icut)),]
  
  miss <- rowSums(is.na(IBS))
  
  #interpolate missing--see snpR::format_snps for the source of this code
  dat <- as.matrix(IBS[,-c(1:ncol(meta))])
  afs <- colMeans(dat, TRUE)
  temp <- which(is.na(dat))/nrow(dat)
  fill_element <- floor(temp) + 1 #get the column for each missing data point
  fill_element[which(temp %% 1 == 0)] <- fill_element[which(temp %% 1 == 0)] - 1 #correct for anything in the last row
  dat[which(is.na(dat))] <- afs[fill_element] #fill with the appropriate allele frequency.
  
  IBS <- cbind(as.data.frame(cbind(IBS[,1:3], miss = miss)), dat)
  return(IBS)
}

##reformat and filter
IBS <- prepIBS(IBS[,-c(1:4)], .5, .5, data.frame(samp = samps, pop = pops, fam = fams))

saveRDS(IBS, "pa_IBS.RDS")


#import and process Posteriors

##import
Posts <- read.table("rapPosts_clean.geno")
Posts <- Posts[,-c(1:2)]
meta <- data.frame(samp = samps, pop = pops, fam = fams)

##remove crappy samples

###get the expected number of each allele
PG1 <- t(Posts[,seq(1,ncol(Posts), 3)])
PG2 <- t(Posts[,seq(2,ncol(Posts), 3)])
PG3 <- t(Posts[,seq(3,ncol(Posts), 3)])

###figure out likely missing data
missing <- ifelse(round(PG1, 3) == 0.333 & round(PG2, 3) == 0.333 & round(PG3, 3) == 0.333, 1, 0)

###get the expected number of each allele
PMaj <- PG1*2 + PG2
PMin <- PG3*2 + PG2
remove(PG2)

###remove crappy alleles
cuttoff <- 0.5
cl <- which(colSums(missing)/nrow(missing) > (1-cuttoff))
if(length(cl) > 0){
  PMaj <- PMaj[,-cl]
  PMin <- PMin[,-cl]
  missing <- missing[,-cl]
}

###remove crappy samples
cuttoff <- 0.5
cl <- which(rowSums(missing)/ncol(missing) > (1-cuttoff))
if(length(cl) > 0){
  PMaj <- PMaj[-cl,]
  PMin <- PMin[-cl,]
  meta <- meta[-cl,]
  missing <- missing[-cl,]
}

###rebind
PGs <- cbind(meta, miss = rowSums(missing), as.data.frame(cbind(PMaj, PMin)))

##save
saveRDS(PGs, "pa_Posts.RDS")


