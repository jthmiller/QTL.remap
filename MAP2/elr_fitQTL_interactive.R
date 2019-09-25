#!/bin/R
### Map QTLs 1 of 3
#debug.cross <- T
#source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library('qtl')
pop <- 'ELR'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

################################################################################
## put chromosomes together
################################################################################

file_list <- gsub('.csv','',list.files(mpath, '*downsmpl_map*'))
file_list <- list.files(mpath, '*downsmpl_map*',include.dirs = T)

chr <- gsub("ELR_gts_CHR",'',file_list)
chr <- as.numeric(gsub("_downsmpl_map",'',chr))

elr <- lapply(file_list,function(X){ read.cross(file=X,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)})

gnos <- lapply(elr,function(X){
 data.frame(X[[1]][[1]][['data']],stringsAsFactors=F)
})
gnos <- do.call(cbind,gnos)
gnos <- cbind(elr[[1]]$pheno,gnos)
gnos$ID <- as.character(gnos$ID)

m_names <- unlist(sapply(elr,function(X){
 markernames(X)
}))

colnames(gnos) <- c('Pheno','sex','ID','bin','pheno_norm',m_names)
rownames(gnos) <- elr[[1]]$pheno$ID

map <- c(colnames(elr[[1]]$pheno),unname(unlist(sapply(elr,pull.map))))
chr <- c(colnames(elr[[1]]$pheno),gsub(":.*","",m_names))
info <- c(colnames(elr[[1]]$pheno),m_names)
headers <- rbind(info,chr,map)
colnames(headers) <- headers[1,]
headers[2:3,1:5] <- ''

headers.u <- unname(data.frame(headers,row.names=NULL,stringsAsFactors=FALSE))
gnos.u <- unname(data.frame(lapply(gnos, as.character),row.names=NULL,stringsAsFactors=FALSE))
colnames(headers.u) <- colnames(gnos.u) <- headers.u[1,]
to_write <- rbind(headers.u,gnos.u)


write.table(to_write,file.path(mpath,'elr.mapped.1_24.csv'),sep=',',row.names=F,quote=F,col.names = F)


fl <- file.path(mpath,'elr.mapped.1_24.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

cross <- sim.geno(cross)


################################################################################
## SCAN
################################################################################

cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno))

################################################################################
perms.bin.em <- scanone(cross, method = "em", model = "binary", maxit = 100,
  n.perm = 100, pheno.col = 4, n.cluster = 6)
perms.norm.em <- scanone(cross, method = "em", model = "normal", maxit = 100,
  n.perm = 100, pheno.col = 5, n.cluster = 6)

## binary
scan.bin.em <- scanone(cross, method = "em", model = "binary", pheno.col = 4)
scan.bin.imp <- scanone(cross, method = "imp", model = "binary", pheno.col = 4)
scan.bin.mr <- scanone(cross, method = "mr", model = "binary", pheno.col = 4)

## normal
scan.norm.em <- scanone(cross, method = "em", model = "normal", pheno.col = 5)
scan.norm.mr <- scanone(cross, method = "mr", model = "normal", pheno.col = 5)
scan.norm.imp <- scanone(cross, method = "imp", model = "normal", pheno.col = 5)
scan.norm.ehk <- scanone(cross, method = "ehk", model = "normal", pheno.col = 5)
## non-parametric
scan.np.bin.em <- scanone(cross, method = "em", model = "np", pheno.col = 4, maxit = 5000)
scan.np.norm.em <- scanone(cross, method = "em", model = "np", pheno.col = 5, maxit = 5000)

bins <- data.frame(summary(scan.bin.em), imp=summary(scan.bin.imp)[,'lod'],mr=summary(scan.bin.mr)[,'lod'],np=summary(scan.np.bin.em)[,'lod'])
norms <- data.frame(summary(scan.norm.em), imp=summary(scan.norm.imp)[,'lod'],mr=summary(scan.norm.mr)[,'lod'],np=summary(scan.np.norm.em)[,'lod'],ehk=summary(scan.norm.ehk)[,'lod'])

## step-wise
full.norm <- stepwiseqtl(cross, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=6)
full.bin <- stepwiseqtl(cross, additive.only = T, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=6)

## manual
qtl <- makeqtl(cross, chr=c(18,13), pos=c(3.63,119.685701),  what="draws")
fitted <- fitqtl(cross,qtl=qtl, formula=y~Q1+Q2)

out.a <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2, method="imp",model="binary",pheno.col=4)
out.ia <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2+Q1:Q3, method="imp",model="binary",pheno.col=4)
out.ib <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2+Q2:Q3, method="imp",model="binary",pheno.col=4)

qtl <- makeqtl(cross, chr=c(13), pos=c(119.68569),  what="draws")
out.i.18 <- addqtl(cross, qtl=qtl, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=5)
out.a.18 <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=5)
####################################################################################

6 and 8 allude to interaction
6 and 18 as well

max(try, what=c( "full"))

 sapply(what=c("best", "full", "add", "int")

####################################################################################
z <- cbind(summary(scan.norm.em),a=summary(scan.norm.mr)[,3],b=summary(scan.norm.imp)[,3])
z <- cbind(z, (z$lod + z$a + z$b)/3)
cbind(summary(scan.bin.em),a=summary(scan.norm.em)[,3],b=summary(scan.bin.mr)[,3])

full.norm <- stepwiseqtl(cross, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=6)
full.bin <- stepwiseqtl(cross, additive.only = T, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=6)

full.norm <- stepwiseqtl(cross, additive.only = T, model='normal', method = "hk", pheno.col = 5, scan.pairs = T, max.qtl=6)
full.bin <- stepwiseqtl(cross, additive.only = T, model='binary', method = "hk", pheno.col = 4, scan.pairs = T, max.qtl=6)


cross <- sim.geno(cross)
full.norm <- stepwiseqtl(cross, additive.only = F, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=6)
full.bin <- stepwiseqtl(cross, additive.only = F, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=6)
save.image(file.path(mpath,'scans.elr.rsave'))

################################################################################
cross <- sim.geno(cross)
qtl <- makeqtl(cross, chr=c(18,13), pos=c(119.685701,119.685701),  what="draws")
out.i.18 <- addqtl(cross, qtl=qtl, formula=y~Q1*Q2, method="imp",model="binary",pheno.col=4)
out.a.18 <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2*Q3, method="imp",model="binary",pheno.col=4)

cross <- sim.geno(cross)
qtl <- makeqtl(cross, chr=c(13), pos=c(119.68569),  what="draws")
out.i.18 <- addqtl(cross, qtl=qtl, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=5)
out.a.18 <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=5)

cross <- sim.geno(cross)
qtl <- makeqtl(cross, chr=c(13,18,23), pos=c(90.5,8.3,66.4), what="draws")

## scan for more additive qtl
out.a <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2+Q3, method="imp",model="normal",pheno.col=5)
### suggests weak qtl on 2?
qtl2 <- makeqtl(cross, chr=c(13,18,23,2), pos=c(90.5,8.3,66.4,52.6229), what="draws")
out.add.4 <- addqtl(cross, qtl=qtl2, formula=y~Q1+Q2+Q3+Q4, method="imp",model="normal",pheno.col=5)


qtl <- makeqtl(cross, chr=1,pos=78,what="draws")


qtl <- makeqtl(cross, chr=c(13,18,23), pos=c(90.5,8.3,66.4), what="draws")
ahr_add <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2+Q3+Q2:Q3, method="imp",model="normal",pheno.col=5)
fitqtl(cross,qtl=qtl, formula=y~Q1+Q2+Q3+Q2:Q3)


qtl <- makeqtl(cross, chr=c(13,18,23), pos=c(90.5,8.3,66.4), what="draws")
ahr_add <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2+Q2:Q3, method="imp",model="normal",pheno.col=5)
fitqtl(cross,qtl=qtl, formula=y~Q1+Q2+Q2:Q3)

cross <- calc.genoprob(cross, error=0.01)



### simple additive
qtl <- makeqtl(cross, chr=c(13,18,23), pos=c(90.5,8.3,66.4), what="draws")
q3.add <- fitqtl(cross,qtl=qtl, formula=y~Q1+Q2+Q3)




qtl <- makeqtl(cross, chr=c(1,13,18,23), pos=c(78,90.5,8.3,66.4))
q4.add <- fitqtl(cross,qtl=qtl, formula=y~Q1+Q2+Q3+Q4)
q4.add.int <- addint(cross, qtl=qtl, formula=y~Q1+Q2+Q3+Q4)


qtl <- makeqtl(cross, chr=c(1,13,18,23,2), pos=c(78,90.5,8.3,66.4,52.6229), what="draws")
ahr_add <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2+Q3+Q4+Q5, method="imp",model="normal",pheno.col=5)
addint(cross, qtl=qtl, formula=y~Q1+Q2+Q3+Q4+Q5)







ahr_add <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2+Q3+Q4, method="imp",model="normal",pheno.col=5)



addint(cross,qtl=qtl2, formula=y~Q1+Q2+Q3+Q4)

out.aq <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2+Q3+Q1:Q2)

addint(cross, qtl=qtl, formula=y~Q1+Q2+Q3+Q1:Q2)

fitqtl(cross,qtl=qtl2)


int <- summary(out.i.18)
add <- summary(out.a.18)

cbind(int,add, int$lod-add$lod)
################################################################################
chr=c(1,13,18,23), pos=c(78,90.5,8.3,66.4)

    # take out several QTLs
    qc <- c(1, 8, 13)
    fake.f2 <- subset(fake.f2, chr=qc)
    # imputate genotypes
    fake.f2 <- calc.genoprob(cross, step=5, err=0.001)
    # 2-dimensional genome scan with additive 3-QTL model
    pos <- list(c(15,35), c(45,65), 28)

chr <- c(1,13,18,23)
pos <- list(c(25,90), c(45,100), c(0,28), c(40,60))

    result <- scanqtl(fake.f2, pheno.col=5, chr=chr, pos=pos,
                      formula=y~Q1+Q2+Q3+Q4, method="imp")


png(paste0('~/public_html/ELR_MR.png'),width=2000,height=2000)
plot(try,chr=c(1,2,6,8,13,18,23),col.scheme ="redblue")
dev.off()

col.scheme = c("viridis", "redblue"

try <- scantwo(cross, pheno.col=5, model="normal",
            method="mr",
            addcovar=NULL, intcovar=NULL, weights=NULL,
            use="complete.obs",
            incl.markers=TRUE, clean.output=FALSE,
            clean.nmar=1, clean.distance=0,
            maxit=1000, tol=1e-4,
            verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
            assumeCondIndep=FALSE, batchsize=250, n.cluster=6)


scntwo <- summary(try,what='best')

scntwo[order(scntwo$lod.int,decreasing=T),]

c1 :c17 5.35e+01  32.6456     8.66    5.92   5.111    52.360  53.70    3.55


    # image of the results
    # chr locations
    chr1 <- as.numeric(matrix(unlist(strsplit(colnames(result),"@")),
                       ncol=2,byrow=TRUE)[,2])
    chr8 <- as.numeric(matrix(unlist(strsplit(rownames(result),"@")),
                       ncol=2,byrow=TRUE)[,2])
    image(chr1, chr8, t(result), las=1, col=rev(rainbow(256,start=0,end=2/3)))


################################################################################


cross.bk <- cross
cross <- switchAlleles(cross, markers = markernames(cross,chr=8))


> find.marker(cross,c(1,17),pos=c( 5.35e+01 ,  32.6456))
[1] "1:16357536"  "17:15233568"


png(paste0('~/public_html/ELR_multi_imputation_interaction_lod.png'))
plotPXG(cross,c("17:15233568","1:16357536"),pch=18,jitter=2,infer=F,pheno.col=4)
dev.off()

png(paste0('~/public_html/ELR_MR.png'),width=2000,height=2000)
plot(try)
dev.off()

png(paste0('~/public_html/ELR_multi_imputation_interaction_lod.png'))
plotPXG(cross,c('8:16768182'),pch=18,jitter=2,infer=F)
dev.off()

2:33206497
23:29198861
'8:16768182'
'8:35102512'
c('18:20509118','13:13194338')
c('13:3456858','20:2544767')

geno.crosstab(cross.bk,c('13:3456858','23:30443925'))

geno.crosstab(cross.bk,c('23:29198861','2:33206497'))


scan.norm.em <- scanone(cross.18, method = "em", model = "normal", maxit = 5000,
  pheno.col = 5)

### Normal scan on transformed phenotype w/Extended haley knott (better for
### selective/missing genos at non-gted ind)
scan.norm.ehk <- scanone(cross, method = "ehk", model = "normal", maxit = 5000, pheno.col = 5)

### Normal scan on transformed phenotype fast haley knott (not robust to missing
### data. LOD inflation)
scan.norm.hk <- scanone(cross.18, method = "hk", model = "normal", maxit = 5000,
  pheno.col = 6)

################################################################################
################################################################################























##############################
##############################
pos <- as.numeric(gsub(".*:","",rownames(gt.missing)))
names(pos) <- rownames(gt.missing)
head(sort(abs(pos -  343835)))

1:317181

1:363497

crs.bk

chr1gts <- pull.geno(crs.bk, 1)

chr1phn <- pull.pheno(crs.bk, 1)



chr1gts <- pull.geno(cross.18, 1)

chr1phn <- pull.pheno(cross.18, 1)


AHR <- cbind(chr1phn,chr1gts[,'1:317181'],chr1gts[,'1:363497'])

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

print("Removing duplicates")
##dups <- findDupMarkers(cross.18, exact.only = F, adjacent.only = F)
##cross.18 <- drop.markers(cross.18, unlist(dups))
##confirm ahr2a 343745   343931 AHR2a
##mid is 343835
