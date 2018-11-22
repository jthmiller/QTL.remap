#!/bin/R
debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")

# NBH <- new.env()
# load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NBH/REMAPS/QTLmap.Rsave', envir =
# NBH) ELR <- new.env()
# load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/ELR/REMAPS/QTLmap.Rsave', envir =
# ELR) NEW <- new.env()
# load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NEW/REMAPS/QTLmap.Rsave', envir =
# NEW)

popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NBH/REMAPS"
cross.NBH <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv", 
  sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))

popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NEW/REMAPS"
cross.NEW <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv", 
  sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))

popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/ELR/REMAPS"
cross.ELR <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv", 
  sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))

cross.ELR$pheno$ID_pop <- "ELR"
cross.NEW$pheno$ID_pop <- "NEW"
cross.NBH$pheno$ID_pop <- "NBH"

#### Sim and reduce to grid ####
NBH.x <- qtl::clean(cross.NBH)
NEW.x <- qtl::clean(cross.NEW)
ELR.x <- qtl::clean(cross.ELR)

cross <- c(NBH.x, NEW.x)
cross$pheno$pheno_norm <- nqrank(cross$pheno$pheno_05)

# POS.map.18 <- est.map(cross, error.prob = 0.025, map.function = 'kosambi',
# maxit = 1000, n.cluster=slurmcore) cross <- replace.map(cross,POS.map.18)

## Use rqtl2 to determine map
cross2 <- convert2cross2(cross)
cross2.map <- est_map(cross2, error_prob = 0.01, map_function = "kosambi", lowmem = FALSE, 
  maxit = 10000, tol = 1e-04, quiet = TRUE, save_rf = T, cores = slurmcore)


### ln model between position and bp for each chrm
models <- lapply(cross2.map, function(x) {
  pos <- as.numeric(x)
  bp <- as.numeric(gsub("*.:", "", names(x)))
  chr <- data.frame(pos, bp)
  mod <- lm(pos ~ bp, data = chr)
})

## positions to predict CM position for
positions <- lapply(cross2.map, function(x) {
  data.frame(bp = as.numeric(unlist(gsub("*.:", "", names(x)))))
})

## Test; predict CM position at map markers
predicted <- mapply("predict", models, positions)
## 91% correlation
elr.positions <- lapply(ELR.x$geno, function(x) {
  data.frame(bp = as.numeric(unlist(gsub("*.:", "", names(x$map)))))
})
elr.predicted <- mapply("predict", models, elr.positions)

map <- pull.map(cross.ELR)
marks <- markernames(cross.ELR)

for (i in 1:24) {
  names(elr.predicted[[i]]) <- names(cross.ELR$geno[[i]]$map)
  # class(elr.predicted[[i]]) <- 'A' marks <- names(cross.ELR$geno[[i]]$map)
  # map[[i]][marks] <- elr.predicted[[i]]
}

# Try to replace ELR map
library(plyr)
df <- ldply(elr.predicted, data.frame)
rownames(df) <- markernames(cross.ELR)
elr.map <- table2map(df)
zero.map <- shiftmap(elr.map)
cross.ELR <- replacemap(cross.ELR, zero.map)

######################### 

## Cross object of all to place psuedomarkers
cross.all <- c(NBH.x, NEW.x, ELR.x)
cross.all$pheno$pheno_norm <- nqrank(cross.all$pheno$pheno_05)

## Simulate genos with high error rate for ELR (should smooth out genotyping
## error)
cross.high <- sim.geno(cross.all, n.draws = 250, step = 1, off.end = 0, error.prob = 0.08, 
  map.function = "kosambi", stepwidth = "fixed")

## Simulate genos with low error rate for NBH (should smooth out genotyping error)

cross <- replace.map(cross, cross2.map)
cross.low <- sim.geno(cross, n.draws = 250, step = 1, off.end = 0, error.prob = 0.01, 
  map.function = "kosambi", stepwidth = "fixed")

save.image("~/NEW.NBH.MAP.QTLmap.Rsave")

cross.nbh <- subset(cross.low, ind = cross.low$pheno$cross1 == 1)
cross.new <- subset(cross.low, ind = cross.low$pheno$cross2 == 1)
cross.elr <- subset(cross.high, ind = cross.high$pheno$cross3 == 1)

nbh.grid <- reduce2grid(cross.nbh)
new.grid <- reduce2grid(cross.new)
elr.grid <- reduce2grid(cross.elr)

scan.norm.imp.NBH.gd <- scanone(nbh.grid, method = "imp", model = "normal", pheno.col = 9)
scan.norm.imp.NEW.gd <- scanone(new.grid, method = "imp", model = "normal", pheno.col = 9)
scan.norm.imp.ELR.gd <- scanone(elr.grid, method = "imp", model = "normal", pheno.col = 10)

scan.norm.imp.NBH <- scanone(cross.nbh, method = "imp", model = "normal", pheno.col = 9)
scan.norm.imp.NEW <- scanone(cross.new, method = "imp", model = "normal", pheno.col = 9)
scan.norm.imp.ELR <- scanone(cross.elr, method = "imp", model = "normal", pheno.col = 10)

save.image("~/NEW.NBH.MAP.QTLmap.Rsave")

### Melt the data to plot in ggplot2
melted.nbh <- data.frame(pop = "NBH", chr = scan.norm.imp.NBH$chr, pos = scan.norm.imp.NBH$pos, 
  lod = scan.norm.imp.NBH$lod)
melted.new <- data.frame(pop = "NEW", chr = scan.norm.imp.NEW$chr, pos = scan.norm.imp.NEW$pos, 
  lod = scan.norm.imp.NEW$lod)
melted.elr <- data.frame(pop = "ELR", chr = scan.norm.imp.ELR$chr, pos = scan.norm.imp.ELR$pos, 
  lod = scan.norm.imp.ELR$lod)
melted <- rbind(melted.nbh, melted.new, melted.elr)

melted$pop <- factor(melted$pop, levels = rev(c("NBH", "NEW", "ELR")))

##### Make QTL Objects to Estimate Effect Size under 2 and full QTL models

nbh.q1.pos <- summary(scan.norm.imp.NBH)[c(2, 18), 2]
nbh.q1.chr <- summary(scan.norm.imp.NBH)[c(2, 18), 1]
nbh.2.qtl <- makeqtl(nbh.grid, chr = nbh.q1.chr, pos = nbh.q1.pos)
nbh.2.fit <- fitqtl(nbh.grid, qtl = nbh.2.qtl, formula = y ~ Q1 + Q2)

new.q1.pos <- summary(scan.norm.imp.NEW.gd)[c(2, 18), 2]
new.q1.chr <- summary(scan.norm.imp.NEW.gd)[c(2, 18), 1]
new.2.qtl <- makeqtl(cross.new, chr = new.q1.chr, pos = new.q1.pos)
new.2.fit <- fitqtl(cross.new, qtl = new.2.qtl, formula = y ~ Q1 + Q2)

elr.q1.pos <- summary(scan.norm.imp.ELR)[c(2, 18), 2]
elr.q1.chr <- summary(scan.norm.imp.ELR)[c(2, 18), 1]
elr.2.qtl <- makeqtl(elr.grid, chr = elr.q1.chr, pos = elr.q1.pos)
elr.2.fit <- fitqtl(elr.grid, qtl = elr.2.qtl, formula = y ~ Q1 + Q2)



########## Expression ######

PAH <- read.table("/home/jmiller1/blast_probes/PAH.txt", stringsAsFactors = F)
PAH$pos <- as.numeric(PAH$pos)
PAH <- PAH[PAH$chr %in% 1:24, ]
PAH$chr <- as.numeric(PAH$chr)
newpos.pah.a <- numeric(length = nrow(PAH))
newpos.pah.b <- numeric(length = nrow(PCB))

for (i in seq_along(PAH$pos)) {
  chr <- map[which(map$chr == PAH$chr[i]), ]
  chr$pos <- as.numeric(chr$pos)
  chr$bp <- as.numeric(chr$bp)
  mod <- lm(pos ~ bp, data = chr)
  new <- data.frame(bp = as.numeric(PAH$pos[i]))
  newpos.pah.a[i] <- predict(mod, newdata = new)
  new <- data.frame(bp = as.numeric(PAH$posb[i]))
  newpos.pah.b[i] <- predict(mod, newdata = new)
}
PAH$pos1 <- newpos.pah.a
PAH$pos2 <- newpos.pah.b
PAH$chem <- "PAH"

PCB <- read.table("/home/jmiller1/blast_probes/PCB.txt", stringsAsFactors = F)
PCB$pos <- as.numeric(PCB$pos)
PCB <- PCB[PCB$chr %in% 1:24, ]
PCB$chr <- as.numeric(PCB$chr)
newpos.pcb.a <- numeric(length = nrow(PCB))
newpos.pcb.b <- numeric(length = nrow(PCB))

for (i in seq_along(PCB$pos)) {
  chr <- map[which(map$chr == PCB$chr[i]), ]
  chr$pos <- as.numeric(chr$pos)
  chr$bp <- as.numeric(chr$bp)
  mod <- lm(pos ~ bp, data = chr)
  new <- data.frame(bp = as.numeric(PCB$pos[i]))
  newpos.pcb.a[i] <- predict(mod, newdata = new)
  new <- data.frame(bp = as.numeric(PCB$posb[i]))
  newpos.pcb.b[i] <- predict(mod, newdata = new)
}
PCB$pos1 <- newpos.pcb.a
PCB$pos2 <- newpos.pcb.b
PCB$chem <- "PCB"

melting <- rbind(PCB, PAH)

##### AHR gene location ##### use linear interpolation to find location ###

AHRs <- read.table("~/QTL_Map_Raw/popgen/tables/AHRs.lifted", stringsAsFactors = F, 
  header = T)
AHRs$chrom <- as.numeric(gsub("chr", "", AHRs$chrom))
AHRs$str <- as.numeric(AHRs$str)
AHRs$stp <- as.numeric(AHRs$stp)
AHRs <- AHRs[!is.na(AHRs$chrom), ]
for (i in seq_along(AHRs$str)) {
  chr <- map[which(map$chr == AHRs$chrom[i]), ]
  chr$pos <- as.numeric(chr$pos)
  chr$bp <- as.numeric(chr$bp)
  mod <- lm(pos ~ bp, data = chr)
  new <- data.frame(bp = as.numeric(AHRs$str[i]))
  AHRs$str[i] <- predict(mod, newdata = new)
  new <- data.frame(bp = as.numeric(AHRs$stp[i]))
  AHRs$stp[i] <- predict(mod, newdata = new)
}

####### 














# save.image('~/NEW.NBH.ELR.QTLmap.Rsave')

load("~/NEW.NBH.ELR.QTLmap.Rsave")
