#!/bin/R
### Map QTLs 1 of 3
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
## For plotting
marker_dens <- list()

# Table of Chroms with sig QTLs
test.QTLs <- read.table(file.path(basedir, "rQTL/metadata/QTLs.txt"), sep = "\t", 
  header = T)

## Get chrom number vector
test.QTLs$chrm.n <- gsub("chr", "", test.QTLs$chrom)

print(paste(pop, X, sep = " "))
############ 

## read in the QTL cross
cross.18 <- read.cross.jm(file = file.path(indpops, paste(pop, ".unphased.f2.csvr", 
  sep = "")), format = "csvr", geno = c(1:3), estimate.map = FALSE)

### Pull names from plinkfile
path <- file.path(indpops, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross.18$pheno$ID <- paste(popname, indname, sep = "_")

## Subset and drop parents
cross.pars <- subset(cross.18, ind = is.na(cross.18$pheno$Phen))

## Remove problematic individuals (found by kinship analysis)
con <- file(file.path(popdir, "kinship.keep.ind.txt"), open = "r")
keepers <- readLines(con)
close(con)

print("Dropping kinship outliers")
cross.18 <- subset(cross.18, ind = cross.18$pheno$ID %in% keepers)
cross.18 <- subset(cross.18, ind = !is.na(cross.18$pheno$Phen))

marker.warning()

print("Dropping all chromosomes except the one to map")
## Map each QTL chro independently
if (mapped.only == T) {
  allbut <- c(1:24)[-X]
  subset.qtl <- chrnames(cross.18)[!chrnames(cross.18) %in% allbut]
  cross.18 <- subset(cross.18, chr = subset.qtl)
  cross.pars <- subset(cross.pars, chr = subset.qtl)
  marker.warning()
}

print("marker dump")

if (pop == "ELR") {
  con <- file(file.path(dirso, "elr.dump.txt"))
  dump <- readLines(con)
  close(con)
  cross.18 <- drop.markers(cross.18, dump)
  
  fileConn <- file(file.path(dirso, "elr.down.txt"))
  elr.down <- readLines(fileConn)
  close(fileConn)
  cross.18 <- drop.markers(cross.18, markernames(cross.18)[markernames(cross.18) %in% 
    elr.down])
  marker.warning()
  
  print("Removing duplicates")
  dups <- findDupMarkers(cross.18, exact.only = F, adjacent.only = F)
  cross.18 <- drop.markers(cross.18, unlist(dups))
  
} else {
  con <- file(file.path(dirso, "nor.dump.txt"))
  dump <- readLines(con)
  close(con)
  cross.18 <- drop.markers(cross.18, dump)
}

marker.warning()

dirso <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/data/"
con <- file(file.path(dirso, "toSwitch.txt"))
toSwitch <- readLines(con)
close(con)

cross.18 <- switchAlleles(cross.18, toSwitch)

gt.missing <- geno.table(cross.18)

gt.cross.pars <- geno.table(cross.pars)

ind <- rownames(gt.cross.pars)[which(gt.cross.pars$AA == 2 | gt.cross.pars$BB == 
  2)]
ind <- ind[ind %in% rownames(gt.missing)]
cross.18 <- drop.markers(cross.18, ind)

marker.warning()

print("pval filter")
gt.missing <- geno.table(cross.18)
cross.18 <- drop.markers(cross.18, rownames(gt.missing[gt.missing$P.value < cutoff, 
  ]))

marker.warning()

print(paste("dropping marks more than", miss, "missing"))
cross.18 <- drop.missing(cross.18, miss)

marker.warning()

gt.cross.pars <- geno.table(cross.pars)[markernames(cross.18), ]

swit <- checkAlleles(cross.18, threshold = 6, verbose = F)
swit.v <- swit[!swit[, 1] %in% toSwitch, 1]
cross.18 <- switchAlleles(cross.18, swit.v)

swit <- checkAlleles(cross.18, threshold = 7, verbose = F)
swit <- swit[!swit[, 1] %in% toSwitch, ]
swit.m <- as.matrix(dist(swit$index))
swit.v <- which(swit.m == 1, arr.ind = TRUE)
swit.v <- apply(swit.v, 1, max)
cross.18 <- switchAlleles(cross.18, swit[unique(swit.v[1]), 1])
swit <- checkAlleles(cross.18, threshold = 10, verbose = F)[, 1]
cross.18 <- drop.markers(cross.18, swit)

cross.18 <- formLinkageGroups(cross.18, max.rf = 0.25, min.lod = 12, reorgMarkers = TRUE)
cross.18 <- switchAlleles(cross.18, markernames(cross.18, chr = 1))
cross.18 <- formLinkageGroups(cross.18, max.rf = 0.25, min.lod = 12, reorgMarkers = TRUE)
cross.18 <- switchAlleles(cross.18, markernames(cross.18, chr = 1))
cross.18 <- subset(cross.18, chr = which.max(nmar(cross.18)))
names(cross.18$geno) <- X

print("estimating map with markers at physical positions")
ord <- order(as.numeric(gsub(paste(X, ":", sep = ""), "", markernames(cross.18, chr = X))))

cross.18 <- switch.order(cross.18, chr = X, ord, error.prob = 0.01, map.function = "kosambi", 
  maxit = 1000, tol = 0.001, sex.sp = F)

marker.warning()

cross.18 <- removeDoubleXO(cross.18, verbose = T)
print("Done removing dxo..")

cross.18 <- drop.missing(cross.18, miss)

marker.warning()

POS.map.18 <- est.map(cross.18, error.prob = 0.1, map.function = "kosambi", chr = X, 
  maxit = 1000)

cross.18 <- replace.map(cross.18, POS.map.18)

print(summary(pull.map(cross.18))[as.character(X), ])
cross.7 <- cross.18

print("Writing the markers to rQTL format")
write.cross(cross.18, filestem = paste(popdir, "/chr", X, "_", outname, ".manymarks.QTLmap", 
  sep = ""), format = "csv", chr = X)
