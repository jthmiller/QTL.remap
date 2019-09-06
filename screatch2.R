fix genotypes ER
cross.18 <- read.cross.jm(file = file.path(indpops, paste(pop, ".unphased.f2.csvr",
  sep = "")), format = "csvr", geno = c(1:3), estimate.map = FALSE)

### Pull names from plinkfile
path <- file.path(indpops, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross.18$pheno$ID <- paste(popname, indname, sep = "_")

## Subset and drop parents
cross.pars <- subset(cross.18, ind = is.na(cross.18$pheno$Phen))
cross.pars <- subset(cross.pars, ind = c(T,F))
gtp <- geno.table(cross.pars)
otp <- geno.table(subset(cross.18, ind = 'NEW_11296'))

A.match <- rownames(gtp[gtp$AA == otp$AA,])
A.no.match <- rownames(gtp[!gtp$AA == otp$AA,])

before.flip <- sapply(1:24,function(X){

    z <- gtp[which(gtp$chr==X),]
    y <- otp[which(otp$chr==X),]

    if((sum(z$BB == y$BB) + sum(z$AA == y$AA)) / (sum(z$BB == y$AA) + sum(z$AA == y$BB)) > 1){return("B")}
    else {return("E")}
    }
)



otp <- geno.table(subset(cross.18, ind = cross.18$pheno$ID == 'NEW_11296'))

after.flip <- sapply(1:24,function(X){

    z <- gtp[which(gtp$chr==X),]

    y <- otp[which(otp$chr==X),]

    z <- z[rownames(z) %in% rownames(y),]


    if((sum(z$BB == y$BB) + sum(z$AA == y$AA)) / (sum(z$BB == y$AA) + sum(z$AA == y$BB)) > 1){return("B")}
    else {return("E")}
    }
)


flips <- which(!before.flip == after.flip)






debug.cross <- T
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')
library('qtlbim')
clean <- qtl::clean
genotyped.only <- T


load("/home/jmiller1/public_html/ELR.bim.Rsave")

cross.18 <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv",
  sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))
sex <- read.table(file = file.path(dirso, "sex.txt"))
rownames(sex) <- sex$ID
cross.18$pheno$sex <- sex[as.character(cross.18$pheno$ID), 2]
cross.18$pheno$binary <- as.numeric(cross.18$pheno$pheno >= 3)

pheno <- read.csv('~/QTL_Map_Raw/popgen/rQTL/data/PhenoDist.csv')
rownames(pheno) <- paste(pheno$pop_all,pheno$Sample,sep='_')
cross.18$pheno$gt <- pheno[as.character(cross.18$pheno$ID),6]

#### IF GENOTYPED IND ONLY
if (genotyped.only == T) cross.18 <- subset(cross.18, ind = cross.18$pheno$gt ==
  'GT')
#### Pheno to GT/NG
## Remove problematic individuals (found by kinship analysis)
con <- file(file.path(popdir, "kinship.keep.ind.txt"), open = "r")
keepers <- readLines(con)
close(con)

print("Dropping kinship outliers")
cross.18 <- subset(cross.18, ind = cross.18$pheno$ID %in% keepers)
cross.18 <- subset(cross.18, ind = !is.na(cross.18$pheno$phen))

ints <- c(1,2,6,8,9,13,14,18,24)
qtls <- subset(cross.18,chr=ints)
qtls <- formLinkageGroups(qtls, max.rf = 0.15, min.lod = 12, reorgMarkers = TRUE)
names(qtls$geno) <- c(14,6,24,1,2,13,9,8,18)

qtls <- orderMarkers(qtls, chr = 18, window = 5, use.ripple = F, error.prob = 0.025,
  map.function = "kosambi", sex.sp = F, maxit = 15000, tol = 0.001)

POS.map.18 <- est.map(qtls, error.prob = 0.025, map.function = "kosambi", maxit = 30)
qtls <- replace.map(qtls, POS.map.18)


crOb <- qtls
qbData.b <- qb.data(crOb, pheno.col = 5, trait = "binary")
qbModel.b <- qb.model(crOb, epistasis = T, main.nqtl = 2, mean.nqtl = 4, depen = FALSE,
  max.qtl = 0)
mc.b <- qb.mcmc(crOb, qbData.b, qbModel.b, pheno.col = 5, n.iter = 30000)
so.b <- qb.scanone(mc.b, epistasis = T, type.scan = "heritability")
so.b.bf <- qb.scanone(mc.b, epistasis = T, type.scan = "2logBF")
best.b <- qb.BayesFactor.jm(mc.b,items = c("pattern", "nqtl"))
two.b <- qb.scantwo(mc.b)
close.b <- qb.close(mc.b)

qtls <- flip.order(qtls,chr=18)


png("/home/jmiller1/public_html/elr_scan_diagnost.so.png", width = 1000)
plot(qb.hpdone(mc.b,smooth = 5, level = 0.95,profile= "2logBF"))
dev.off()



### Slice on 14
temp <- qb.sliceone(mc.b, slice = 1,smooth=20,weight ="count",type.scan = "cellmean")
png("/home/jmiller1/public_html/slice.14_6.18.png")
plot(temp,chr=c(2,9) ,ylim=c(-2,2))
dev.off()
### Slice on 6
temp <- qb.sliceone(mc.b, slice = 2,smooth=20,weight ="count",type.scan = "cellmean")
png("/home/jmiller1/public_html/slice.6_14.18.png")
plot(temp,chr=c(1,9) ,ylim=c(-2,2))
dev.off()

### slice on 13
temp <- qb.sliceone(mc.b, slice =6, smooth=20,weight ="count",type.scan = "cellmean")
png("/home/jmiller1/public_html/slice.13_24.1.8.18.png")
plot(temp,chr=c(3,4,8,9),ylim=c(-2,2),sub="slice on 13")
dev.off()

temp <- qb.sliceone(mc.b, slice =6, smooth=20,weight ="count",type.scan = "LPD")
png("/home/jmiller1/public_html/slice.13_24.1.8.18.LPDf.png")
plot(temp,chr=c(3,4,8,9),sub="slice on 13")
dev.off()

#slice on 18
temp <- qb.sliceone(mc.b, slice=6,smooth=20,weight ="count",type.scan = "cellmean")
png("/home/jmiller1/public_html/18_1_elr_slice.png")
plot(temp,chr=c(1,2,6,8,14,18,23),ylim=c(-2,2))
dev.off()





#slice on 24
temp <- qb.sliceone(mc.b, slice =6,smooth=20,weight ="count",type.scan = "cellmean")
png("/home/jmiller1/public_html/slice.png")
plot(temp,ylim=c(-2,2))
dev.off()


scan = c("sum","mean","epistasis")
chr=c(1,2,6,8,9,13,18),


png("/home/jmiller1/public_html/elr_coda.png")
plotRF(qtls)
dev.off()

png("/home/jmiller1/public_html/elr_coda.png")
plotRF(qtls)
dev.off()



a <- find.marker(cross.18,13,27 )
b <- find.marker(cross.18, 8, 76)
png("/home/jmiller1/public_html/elr_13by8.png")
effectplot(crOb,pheno.col=1, mname1 = b, mname2 = a,ylim=c(0,5))
dev.off()

a <- find.marker(cross.18,13,27 )
b <- find.marker(cross.18, 8, 76)
png("/home/jmiller1/public_html/elr8by13.png")
effectplot(crOb,pheno.col=1, mname1 = a, mname2 = b,ylim=c(0,5))
dev.off()





png("/home/jmiller1/public_html/NEW_2_pxg.png")
plotPXG(cross.18, a, pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = a)
dev.off()
png("/home/jmiller1/public_html/NEW_18_pxg.png")
plotPXG(cross.18, b, pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = b)
dev.off()








confirm ahr2a 343745   343931 AHR2a
mid is 343835




module load bcftools
bcftools view -r 1:340000-345000 SOMM.vcf.gz | less -S

vcftools --gzvcf SOMM.vcf.gz --chr 'chr1' --from-bp 340000 --to-bp 345000 --stdout | less -S
