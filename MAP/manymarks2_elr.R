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
cross.pars.BI <- subset(cross.18, ind = 'BLI_BI1124M')

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

print("Removing duplicates")
##dups <- findDupMarkers(cross.18, exact.only = F, adjacent.only = F)
##cross.18 <- drop.markers(cross.18, unlist(dups))
##confirm ahr2a 343745   343931 AHR2a
##mid is 343835


marker.warning()

gt.missing <- geno.table(cross.18)


##############################
##############################
pos <- as.numeric(gsub(".*:","",rownames(gt.missing)))
names(pos) <- rownames(gt.missing)
head(sort(abs(pos -  343835)))

1:317181

1:363497

chr1gts <- pull.geno(cross.18, 1)

chr1phn <- pull.pheno(cross.18, 1)

AHR <- cbind(chr1phn,chr1gts[,'1:317181'],chr1gts[,'1:363497'])
AHR <- AHR[order(AHR[,1]),]

table(AHR[AHR[,1]<2,3])
table(AHR[AHR[,1]>2,3])

table(AHR[AHR[,1]==0,3])
table(AHR[AHR[,1]==1,3])
table(AHR[AHR[,1]==4,3])
table(AHR[AHR[,1]==5,3])

chr1.pars <- pull.geno(cross.pars, 1)
rbind(chr1.pars[,'1:317181'], chr1.pars[,'1:363497'])

TAKE THE HOMZYGOUS GENOTYPES FOR THE ONE PARENT AND SEE IF THEY TEND TOWARD 1:2:1 compared to
het in parent.

AHR[,1] <- as.factor(AHR[,1])
AHR[,2] <- as.factor(AHR[,2])
AHR[,3] <- as.factor(AHR[,3])

png('~/public_html/ER_AHR.png')
plot(table(AHR[AHR[,1]<2,2]))
dev.off()

for(
table(AHR[,1])

##############################
##############################

gt.cross.BI <- geno.table(cross.pars.BI)
keep.BI <- rownames(gt.cross.BI)[which(gt.cross.BI$AA==1 | gt.cross.BI$BB==1)]









cutoff <- 1e-04

print("pval filter")
gt.missing <- geno.table(cross.18)
cross.18 <- drop.markers(cross.18, rownames(gt.missing[gt.missing$P.value < cutoff,
  ]))

marker.warning()

print(paste("dropping marks more than", miss, "missing"))
cross.18 <- drop.missing(cross.18, miss)

marker.warning()

gt.cross.pars <- geno.table(cross.pars)[markernames(cross.18), ]
swit <- rownames(gt.cross.pars[gt.cross.pars$AA == 1, ])
cross.pars <- switchAlleles(cross.pars, markers = swit)
cross.18 <- switchAlleles(cross.18, markers = swit)

fileConn <- file(file.path(popdir, paste(X, "BB_markes.txt", sep = "_")))
writeLines(swit, fileConn)
close(fileConn)

print("estimating map with markers at physical positions")
ord <- order(as.numeric(gsub(paste(X, ":", sep = ""), "", markernames(cross.18, chr = X))))

cross.18 <- switch.order(cross.18, chr = X, ord, error.prob = 0.01, map.function = "kosambi",
  maxit = 1000, tol = 0.001, sex.sp = F)

swit <- checkAlleles(cross.18, threshold = 6, verbose = F)
cross.18 <- switchAlleles(cross.18, swit)

cross.18 <- formLinkageGroups(cross.18, max.rf = 0.25, min.lod = 12, reorgMarkers = TRUE)
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

print("Writing the markers to rQTL format")
write.cross(cross.18, filestem = paste(popdir, "/chr", X, "_", outname, ".manymarks.QTLmap",
  sep = ""), format = "csv", chr = X)
