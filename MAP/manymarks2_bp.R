#!/bin/R
### marker filtering manymarks2 debug.cross <- T
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
if (!pop == "BRP") {
  con <- file(file.path(popdir, "kinship.keep.ind.txt"), open = "r")
  keepers <- readLines(con)
  close(con)
  print("Dropping kinship outliers")
  cross.18 <- subset(cross.18, ind = cross.18$pheno$ID %in% keepers)
  cross.18 <- subset(cross.18, ind = !is.na(cross.18$pheno$Phen))
}

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

BP <- subset(cross.18, ind = cross.18$pheno$ID %in% c("BRP_BRP1F", "BRP_BRP8F"))
BI <- subset(cross.18, ind = cross.18$pheno$ID %in% c("BRP_BRP1M", "BRP_BRP8M"))
gt.BI <- geno.table(BI)
gt.BP <- geno.table(BP)
swt <- c(rownames(gt.BP[which(gt.BP$BB == 2), ]), rownames(gt.BI[which(gt.BP$AA == 
  2), ]))
cross.18 <- switchAlleles(cross.18, swt)
BI <- switchAlleles(BI, swt)
BP <- switchAlleles(BP, swt)
gt.BI <- geno.table(BI)
gt.BP <- geno.table(BP)

dirso <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/data/"
con <- file(file.path(dirso, "toSwitch.N.txt"))
toSwitch <- readLines(con)
close(con)

cross.18 <- switchAlleles(cross.18, toSwitch)

cross.18 <- drop.missing(cross.18, miss1)

marker.warning()

### To confirm markers later
g <- pull.geno(cross.18, chr = X)[1, swt]

gt.missing <- geno.table(cross.18)

gt.cross.pars <- geno.table(cross.pars)

marker.warning()

print("pval filter")
gt.missing <- geno.table(cross.18)
cross.18 <- drop.markers(cross.18, rownames(gt.missing[gt.missing$P.value < cutoff, 
  ]))

marker.warning()

print(paste("dropping marks more than", miss, "missing"))
cross.18 <- drop.missing(cross.18, miss1)

marker.warning()

cross.18 <- subset(cross.18, ind = !is.na(cross.18$pheno$Phen))

cross.18 <- formLinkageGroups(cross.18, max.rf = 0.2, min.lod = 12, reorgMarkers = TRUE)
cross.18 <- switchAlleles(cross.18, markernames(cross.18, chr = 1))
cross.18 <- switchAlleles(cross.18, markernames(cross.18, chr = 2))
cross.18 <- formLinkageGroups(cross.18, max.rf = 0.2, min.lod = 12, reorgMarkers = TRUE)
cross.18 <- switchAlleles(cross.18, markernames(cross.18, chr = 1))
cross.18 <- switchAlleles(cross.18, markernames(cross.18, chr = 2))
cross.18 <- formLinkageGroups(cross.18, max.rf = 0.2, min.lod = 12, reorgMarkers = TRUE)
cross.18.1 <- subset(cross.18, chr = 1)
cross.18.2 <- subset(cross.18, chr = 2)
names(cross.18.1$geno) <- X
names(cross.18.2$geno) <- X

if (getmax == T) {
  us <- which.max(c(nmar(cross.18.1), nmar(cross.18.2)))
  cross.18 <- get(c("cross.18.1", "cross.18.2")[us])
} else {
  us <- which.min(c(nmar(cross.18.1), nmar(cross.18.2)))
  cross.18 <- get(c("cross.18.1", "cross.18.2")[us])
}

print("Removing duplicates")
dups <- findDupMarkers(cross.18, exact.only = T, adjacent.only = F)
cross.18 <- drop.markers(cross.18, unlist(dups))

print("estimating map with markers at physical positions")
ord <- order(as.numeric(gsub(paste(X, ":", sep = ""), "", markernames(cross.18, chr = X))))

cross.18 <- switch.order(cross.18, chr = X, ord, error.prob = 0.01, map.function = "kosambi", 
  maxit = 1000, tol = 0.001, sex.sp = F)

marker.warning()

print("Removing duplicates")
dups <- findDupMarkers(cross.18, exact.only = F, adjacent.only = T)
cross.18 <- drop.markers(cross.18, unlist(dups))

marker.warning()

cross.18 <- removeDoubleXO(cross.18, verbose = T)
print("Done removing dxo..")

marker.warning()

cross.18 <- drop.missing(cross.18, miss2)

marker.warning()

print("Re-estimating the map")
POS.map.18 <- est.map(cross.18, error.prob = 0.08, map.function = "kosambi", chr = X, 
  maxit = 3000)

cross.18 <- replace.map(cross.18, POS.map.18)

print(summary(pull.map(cross.18))[as.character(X), ])

print("Writing the markers to rQTL format")
write.cross(cross.18, filestem = paste(popdir, "/chr", X, "_", outname, ".manymarks.QTLmap", 
  sep = ""), format = "csv", chr = X)
