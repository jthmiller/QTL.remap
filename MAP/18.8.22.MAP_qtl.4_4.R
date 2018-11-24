#!/bin/R

### Map QTLs 4 of 3
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
load(paste(popdir, "/chr", X, "_", outname, ".QTLmap.Rsave", sep = ""))

cross.18 <- read.cross(format = "csv", dir = popdir, file = paste("/chr", X, "_", 
  outname, "_3.QTLmap.csv", sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", 
  "B"))

cross.18 <- formLinkageGroups(cross.18, max.rf = 0.25, min.lod = 12, reorgMarkers = TRUE)
cross.18 <- subset(cross.18, chr = which.max(nmar(cross.18)))
names(cross.18$geno) <- X

print("estimating map with markers at physical positions")
ord <- order(as.numeric(gsub(paste(X, ":", sep = ""), "", markernames(cross.18, chr = X))))

if (!reorder.marks == T) {
  cross.18 <- switch.order(cross.18, chr = X, ord, error.prob = 0.025, map.function = "kosambi", 
    maxit = 1000, tol = 0.001, sex.sp = F)
}

print("Re-estimating the final map with many iterations...")
POS.map.18 <- est.map(cross.18, error.prob = 0.025, map.function = "kosambi", chr = X, 
  maxit = 10000)
cross.18 <- replace.map(cross.18, POS.map.18)
print("Done mapping..")

vec <- as.numeric(gsub(paste(X, ":", sep = ""), "", markernames(cross.18)))
print(paste("physical positions from", min(vec), "to", max(vec)))

print("Adding un-genotyped individuals for stratified analysis")
pheno.all <- phen <- read.csv("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/metadata/ALL_phenotype_Dist.csv", 
  header = T)
phen$Pheno_05 <- phen$pheno_all
phen$pheno_all <- as.numeric(trsl.bin[as.character(phen$pheno_all)])
index <- which(phen$pop_all == pop)
rownames(phen) <- paste(phen$pop_all, phen$IND, sep = "_")
set <- rownames(phen)[which(phen$pop_all == pop)]
ngos <- set[which(!set %in% cross.18$pheno$ID)]

## SET genotyped ind label
cross.18$pheno$gt <- "GT"

no_genos <- data.frame(Pheno = phen[ngos, 4], sex = 0, ID = ngos, gt <- "NG", markers = matrix("-", 
  nrow = length(ngos), ncol = as.numeric(nmar(cross.18))))

write.table(no_genos, file = file.path(popdir, paste(X, "no_genos.csv", sep = "_")), 
  col.names = F, row.names = F, quote = F, sep = ",")

write.cross(cross.18, filestem = paste(popdir, "/chr", X, "_", outname, "_4.QTLmap", 
  sep = ""), format = "csv", chr = X)

system(paste("cat ", popdir, "/chr", X, "_", outname, "_4.QTLmap.csv ", popdir, "/", 
  paste(X, "no_genos.csv", sep = "_"), " > ", popdir, "/temp.", X, sep = ""))

print("saving... done with mapping ind chromosomes")

save.image(paste(popdir, "/chr", X, "_", outname, ".QTLmap.Rsave", sep = ""))
