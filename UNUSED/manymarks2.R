### manymarks2

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

dirso <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/data/"

if (!pop == "ELR") {
  con <- file(file.path(dirso, "toSwitch.N.txt"))
  toSwitch <- readLines(con)
  close(con)
} else {
  con <- file(file.path(dirso, "toSwitch.txt"))
  toSwitch <- readLines(con)
  close(con)
}

cross.18 <- switchAlleles(cross.18, toSwitch)

marker.warning()

gt.missing <- geno.table(cross.18)

gt.cross.pars <- geno.table(cross.pars)
