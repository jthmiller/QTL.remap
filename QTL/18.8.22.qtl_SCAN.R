#!/bin/R
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")

cross.18 <- reconst(X = chrms, pop = popq, temp.dir = popdir, a = 2)

cross.18 <- qtl::clean(cross.18)

print("Writing the merged chromosome markers to rQTL format")
write.cross(cross.18, filestem = paste(popdir, "/", outname, ".BACKUP.QTL_chr.QTLmap", 
  sep = ""), format = "csv")

pheno.all <- phen <- read.table("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/metadata/ALL_phenotype_Dist.txt", 
  header = T)
phen$Pheno_05 <- phen$pheno_all
index <- which(phen$pop_all == popq)
count.pheno <- sapply(0:5, function(pt) {
  pt <- as.character(pt)
  total <- sum(phen$Pheno_05[index] == pt)
  incross <- sum(cross.18$pheno$Pheno_05 == pt)
  return(as.numeric(total - incross))
})
names(count.pheno) <- as.character(0:5)
count.pheno <- count.pheno[!is.na(count.pheno)]
count.pheno <- rep.int(names(count.pheno), times = as.numeric(count.pheno))
trsl.bin <- c(0, 0, 0, 1, 1, 1)
names(trsl.bin) <- as.character(0:5)
phenos <- data.frame(count.pheno, Pheno = as.numeric(trsl.bin[as.character(count.pheno)]))
rownames(phenos) <- paste("NG", 1:length(count.pheno), sep = "")

ind.inx <- grep("NG", cross.18$pheno$ID)
repl <- phenos[as.character(cross.18$pheno$ID[grep("NG", cross.18$pheno$ID)]), 1]
cross.18$pheno$pheno_05[ind.inx] <- as.character(repl)
cross.18$pheno$pheno_05 <- as.numeric(cross.18$pheno$pheno_05)
cross.18$pheno$pheno <- as.numeric(cross.18$pheno$pheno_05 >= 3)  ### Split pheotypes

### Ids for strata permutations
cross.18$pheno$stata <- gsub("[0-9]", "", cross.18$pheno$ID)
#### Transform the phenotype for model flexibility. Imputations are better for
#### selective genotyping than HK
cross.18$pheno$pheno_norm <- nqrank(cross.18$pheno$pheno_05)
#### Missing individuals have NA for phenotype
indy <- cross.18$pheno$ID[grep("NG", cross.18$pheno$ID)]
cross.18$pheno$pheno_miss01 <- cross.18$pheno$pheno
cross.18$pheno$pheno_miss05 <- cross.18$pheno$pheno_05
cross.18$pheno$pheno_miss01[indy] <- NA
cross.18$pheno$pheno_miss05[indy] <- NA
### Error rate and genoprobs
ers <- 0.002
cross.18 <- calc.genoprob(cross.18, step = 1, error.prob = ers, map.function = "kosambi")
cross.18 <- sim.geno(cross.18, error.prob = ers, step = 0.5, n.draws = 250)
write.cross(cross.18, filestem = paste(qtldir, "NO_DUP_MARKERS.QTLmap", sep = ""), 
  format = "csv")

#### Nonparametric scan
scan.np.em <- scanone(cross.18, method = "em", model = "np", pheno.col = 2, maxit = 5000)
### 2 part model
scan.2p.em <- scanone(cross.18, method = "em", model = "2part", pheno.col = 1, maxit = 5000)
### Imputation on score
scan.norm_05.imp <- scanone(cross.18, method = "imp", model = "normal", maxit = 5000, 
  pheno.col = 2)

### Binary model with EM
scan.bin.em <- scanone(cross.18, method = "em", model = "binary", pheno.col = 1, 
  maxit = 5000)
### Binary model with marker regression
scan.bin.mr <- scanone(cross.18, method = "mr", model = "binary", pheno.col = 1)
### Binary model with haley knott (heavy inflation of LOD)
scan.bin.hk <- scanone(cross.18, method = "hk", model = "binary", pheno.col = 1)
### Binary model without ungenotyped
scan.bin.em.gt <- scanone(cross.18, method = "em", model = "binary", pheno.col = 7, 
  maxit = 5000)

### Normal scan on transformed phenotype
scan.norm.em <- scanone(cross.18, method = "em", model = "normal", maxit = 5000, 
  pheno.col = 6)
### Normal scan on transformed phenotype w/Extended haley knott (better for
### selective/missing genos at non-gt'ed ind)
scan.norm.ehk <- scanone(cross.18, method = "ehk", model = "normal", maxit = 5000, 
  pheno.col = 6)
### Normal scan on transformed phenotype fast haley knott (not robust to missing
### data. LOD inflation)
scan.norm.hk <- scanone(cross.18, method = "hk", model = "normal", maxit = 5000, 
  pheno.col = 6)
### Imputation on transformed phenotype
scan.norm.imp <- scanone(cross.18, method = "imp", model = "normal", pheno.col = 6)
### Imputation on transformed phenotype with missing ind removed
scan.norm.imp.gt <- scanone(cross.18, method = "imp", model = "normal", pheno.col = 8)
### EM scan without ungenotyped
scan.norm.em.gt <- scanone(cross.18, method = "em", model = "normal", maxit = 5000, 
  pheno.col = 8)

## Find LOD thresholds to get a null dist. of remove markers that are exactly the
## same to speed up (same results)

perms.np.em <- scanone(cross.18, model = "np", pheno.col = 2, n.perm = 2000, method = "em", 
  perm.strata = cross.18$pheno$stata, n.cluster = slurmcore)
perms.2p.em <- scanone(cross.18, model = "2part", n.perm = 2000, perm.strata = cross.18$pheno$stata, 
  pheno.col = 2, n.cluster = slurmcore)

perms.bin.em <- scanone(cross.18, method = "em", model = "binary", maxit = 2000, 
  n.perm = 500, perm.strata = cross.18$pheno$stata, pheno.col = 1, n.cluster = slurmcore)
perms.bin.mr <- scanone(cross.18, method = "mr", model = "binary", n.perm = 2000, 
  perm.strata = cross.18$pheno$stata, n.cluster = slurmcore)
# perms.bin.em <- scanone(cross.18,
# method='hk',model='binary',maxit=2000,n.perm=500,perm.strata=cross.18$pheno$stata,pheno.col=1,
# n.cluster=slurmcore)

perms.norm.em <- scanone(cross.18, method = "em", model = "normal", maxit = 500, 
  n.perm = 500, perm.strata = cross.18$pheno$stata, pheno.col = 6, n.cluster = slurmcore)
perms.norm.ehk <- scanone(cross.18, method = "ehk", model = "normal", maxit = 500, 
  n.perm = 500, perm.strata = cross.18$pheno$stata, pheno.col = 6, n.cluster = slurmcore)
perms.norm.hk <- scanone(cross.18, method = "hk", model = "normal", maxit = 500, 
  n.perm = 500, perm.strata = cross.18$pheno$stata, pheno.col = 6, n.cluster = slurmcore)
perms.norm.imp <- scanone(cross.18, method = "imp", model = "normal", n.perm = 500, 
  perm.strata = cross.18$pheno$stata, pheno.col = 6, n.cluster = slurmcore)
perms.norm.imp.2 <- scanone(cross.18, method = "imp", model = "normal", chr = -2, 
  n.perm = 500, perm.strata = cross.18$pheno$stata, pheno.col = 6, n.cluster = slurmcore)

save.image(file.path(popdir, "QTLmap.Rsave"))
## load(file.path(popdir,'QTLmap.Rsave'))

### Multi-QTL models
th <- summary(perms.norm.imp)[1, ]
norm.qtl <- summary(scan.norm.imp, perms = perms.norm.imp, alpha = 0.05)
qtl.uns <- makeqtl(cross.18, chr = norm.qtl$chr, pos = norm.qtl$pos)
full <- stepwiseqtl(cross.18, additive.only = T, method = "imp", pheno.col = 6, scan.pairs = T)

### CIM
out.cim.40 <- cim(cross.18, n.marcovar = 3, window = 40, pheno.col = 6, method = "imp", 
  error.prob = 0.002)
out.cim.inf <- cim(cross.18, n.marcovar = 3, window = Inf, pheno.col = 6, method = "imp", 
  error.prob = 0.002)
# add.cim.covar(out.cim.inf, chr=c(1,4,6,15))

save.image(file.path(popdir, "QTLmap.Rsave"))

cross.18 <- subset(cross.18, ind = cross.18$pheno$stata == "ind")
gt <- geno.table(cross.18)
barks <- rownames(gt[which(gt$missing > 2), ])
cross.18 <- drop.markers(cross.18, markers = barks)
print("Re-estimating the final map with many iterations...")
POS.map.18 <- est.map(cross.18, error.prob = ers, map.function = "kosambi", chr = X, 
  maxit = 1000)
cross.18 <- replace.map(cross.18, POS.map.18)
print("Done mapping..")
cross.18 <- sim.geno(cross.18, error.prob = ers, step = 1, n.draws = 25)
full.2 <- stepwiseqtl(cross.18, qtl = qtl.uns, additive.only = F, method = "imp", 
  pheno.col = 6, scan.pairs = T, keeptrace = FALSE)

Impute <- mqmaugment(cross.18, minprob = 0.01, strategy = "impute", verbose = TRUE)
print("mqmscan additive")
mqm <- mqmscan(Impute, pheno.col = 6, model = "additive", forceML = F, outputmarkers = F)
autocofactors <- mqmautocofactors(Impute, 15)
mqm_auto <- mqmscan(Impute, autocofactors)
mqm.d <- mqmscan(Impute, pheno.col = 6, model = "dominance", forceML = F, outputmarkers = F)
results <- mqmpermutation(Impute, scanfunction = mqmscan, pheno.col = 6, cofactors = autocofactors, 
  n.cluster = slurmcore, n.perm = 24, batchsize = 24)
resultsrqtl <- mqmprocesspermutation(results)

save.image(file.path(popdir, "QTLmap.Rsave"))
