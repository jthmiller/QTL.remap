############ TRY TO MAP ADDnl AHRs #####3

cross.18 <- read.cross(file = file.path(indpops, paste(pop, ".um.unmapped.f2.csvr", 
  sep = "")), format = "csvr", geno = c(1:3), estimate.map = FALSE)

kep <- c(markernames(cross.NBH), markernames(cross.18, chr = "NW_012234311.1"))
drops <- markernames(cross.18)[!markernames(cross.18) %in% kep]
cross.18 <- drop.markers(cross.18, drops)
cross.18 <- formLinkageGroups(cross.18, max.rf = 0.1, min.lod = 12, reorgMarkers = TRUE)

### Pull names from plinkfile
path <- file.path(indpops, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross.18$pheno$ID <- paste(popname, indname, sep = "_")

sex <- read.table(file = file.path(dirso, "sex.txt"))
rownames(sex) <- sex$ID
cross.18$pheno$sex <- sex[as.character(cross.18$pheno$ID), 2]
names(cross.18$pheno) <- c("pheno", "sex", "ID")
cross.18$pheno$binary <- as.numeric(cross.18$pheno$pheno >= 3)
