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

file_list <- list.files(mpath, 'ELR_all_mark_?[0-9]?[0-9]_tsp.csv')

chr <- gsub("ELR_all_mark_",'',file_list)
chr <- as.numeric(gsub("_tsp.csv",'',chr))

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

write.table(to_write,file.path(mpath,'elr.mapped.tsp.csv'),sep=',',row.names=F,quote=F,col.names = F)

################################################################################

################################################################################
fl <- file.path(mpath,'elr.mapped.tsp.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

cross <- sim.geno(cross)
cross <- calc.genoprob(cross,step=1,error.prob=0.01,off.end=5)

## binary
scan.bin.em <- scanone(cross, method = "em", model = "binary", pheno.col = 4)
scan.bin.mr <- scanone(cross, method = "mr", model = "binary", pheno.col = 4)

## normal transformed
scan.normT.em <- scanone(cross, method = "em", model = "normal", pheno.col = 5)
scan.normT.mr <- scanone(cross, method = "mr", model = "normal", pheno.col = 5)
scan.normT.imp <- scanone(cross, method = "imp", model = "normal", pheno.col = 5)
scan.normT.ehk <- scanone(cross, method = "ehk", model = "normal", maxit = 5000, pheno.col = 5)

## normal
scan.norm.em <- scanone(cross, method = "em", model = "normal", pheno.col = 1)
scan.norm.mr <- scanone(cross, method = "mr", model = "normal", pheno.col = 1)
scan.norm.imp <- scanone(cross, method = "imp", model = "normal", pheno.col = 1)
scan.norm.ehk <- scanone(cross, method = "ehk", model = "normal", maxit = 5000, pheno.col = 1)

## non-parametric
scan.np.em.b <- scanone(cross, method = "em", model = "np", pheno.col = 4, maxit = 5000)
scan.np.em.n <- scanone(cross, method = "em", model = "np", pheno.col = 5, maxit = 5000)
#save.image(file.path(mpath,'scans.elr.rsave'))
################################################################################

################################################################################
##SEX
scan.bin.sex <- scanone(cross, method = "em", model = "binary", pheno.col = 2)
################################################################################

################################################################################
## step-wise
full.norm.add_only <- stepwiseqtl(cross, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=8)
################################################################################

################################################################################
qtl.2 <- makeqtl(cross, chr=c(2), pos=c(80),  what="draws")
out.a.2 <- addqtl(cross, qtl=qtl.2, formula=y~Q1+Q2, method="imp",model="binary",pheno.col=4)
out.i.2 <- addqtl(cross, qtl=qtl.2, formula=y~Q1*Q2, method="imp",model="binary",pheno.col=4)
## indicates possible interaction on 13 and 18

qtl.8 <- makeqtl(cross, chr=c(8), pos=c(88),  what="draws")
out.i.8 <- addqtl(cross, qtl=qtl.8, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=4)
out.a.8 <- addqtl(cross, qtl=qtl.8, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=4)

qtl.18 <- makeqtl(cross, chr=c(18), pos=c(26.7),  what="draws")
out.a.18 <- addqtl(cross, qtl=qtl.18, formula=y~Q1+Q2, method="imp",model="binary",pheno.col=4)
out.i.18 <- addqtl(cross, qtl=qtl.18, formula=y~Q1*Q2, method="imp",model="binary",pheno.col=4)
################################################################################

################################################################################
## manual at stepwise positions
qtl.18 <- makeqtl(cross, chr=c(18), pos=c(27),  what="draws")
out.i.18.b <- addqtl(cross, qtl=qtl.18, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=5)
out.a.18.b <- addqtl(cross, qtl=qtl.18, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=5)

qtl.2 <- makeqtl(cross, chr=c(2), pos=c(80),  what="draws")
out.i.2.b <- addqtl(cross, qtl=qtl.2, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=5)
out.a.2.b <- addqtl(cross, qtl=qtl.2, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=5)

qtl.8 <- makeqtl(cross, chr=c(8), pos=c(67),  what="draws")
out.i.8.b <- addqtl(cross, qtl=qtl.8, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=5)
out.a.8.b <- addqtl(cross, qtl=qtl.8, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=5)
################################################################################

####################################################################################
## Multi/interacting QTL

################################################################################
qtl <- makeqtl(cross, chr=c(2,8,18), pos=c(80,67,27),  what="draws")
fitted <- fitqtl(cross,qtl=qtl, formula=y~Q1+Q2+Q3, pheno.col=4)

qtl <- makeqtl(cross, chr=c(2,8,18), pos=c(80,67,27),  what="draws")
fitted <- fitqtl(cross, qtl=qtl,  formula=y~Q1+Q2+Q3, pheno.col=5, method="imp",model="normal")
fitted.noint.add <- addqtl(cross, qtl=qtl,  formula=y~Q1+Q2+Q3, pheno.col=5, method="imp",model="normal")

qtl <- makeqtl(cross, chr=c(2,18), pos=c(80,27))
fitted.int.only <- fitqtl(cross, qtl=qtl,  formula=y~Q1+Q2+Q1:Q2, pheno.col=5, method="imp",model="normal")
out.i <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2+Q3, pheno.col=5, method="imp",model="normal")

qtl <- makeqtl(cross, chr=c(2,18), pos=c(80,27))
out.ap13 <- addqtl(cross, qtl=qtl, formula = y~Q1 + Q1*Q3 + Q2,  model='normal', method = "imp", pheno.col = 5)
out.ap23 <- addqtl(cross, qtl=qtl, formula = y~Q1 + Q2 + Q2*Q3,  model='normal', method = "imp", pheno.col = 5)

qtl <- makeqtl(cross, chr=c(2,3,13,18,19), pos=c(80,11.18,38.09,27,41.3))
out.ap13 <- fitqtl(cross, qtl=qtl, formula = y~Q1 + Q1*Q2 + Q1*Q3 + Q4 + Q4*Q5,  model='normal', method = "imp", pheno.col = 5)

### 19 interacts w 2 (from picture)
### 18 and 23
####################################################################################

save.image(file.path(mpath,'scans.elr.rsave'))

####################################################################################

####################################################################################
## PERMS
perms.norm.imp <- scanone(cross, method = "imp", model = "normal", maxit = 1000,
  n.perm = 1000, pheno.col = 5, n.cluster = 10)

perms.bin.em <- scanone(cross, method = "em", model = "binary", maxit = 1000,
  n.perm = 1000, pheno.col = 4, n.cluster = 10)

perms.norm.em <- scanone(cross, method = "em", model = "normal", maxit = 1000,
  n.perm = 1000, pheno.col = 5, n.cluster = 10)
## LOD thresholds (1000 permutations)
##      lod
## 5%  4.48
## 10% 4.00
####################################################################################

save.image(file.path(mpath,'scans.elr.rsave'))

################################################################################

################################################################################
### LONG #######################################################################
full.norm <- stepwiseqtl(cross, additive.only = F, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=7)
save.image(file.path(mpath,'scans.elr.rsave'))

full.bin <- stepwiseqtl(cross, additive.only = F, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=7)
save.image(file.path(mpath,'scans.elr.rsave'))

full.norm.int <- stepwiseqtl(cross, additive.only = F, model='normal', method = "hk", pheno.col = 5, scan.pairs = T, max.qtl=7)
save.image(file.path(mpath,'scans.elr.rsave'))

full.bin.int <- stepwiseqtl(cross, additive.only = F, model='binary', method = "hk", pheno.col = 4, scan.pairs = T, max.qtl=7)
save.image(file.path(mpath,'scans.elr.rsave'))
################################################################################

################################################################################
gg <- sim.geno(cross,step=4)
gg_step2 <- reduce2grid(gg)

cgb <- calc.genoprob(cross, step=4)
cgb_step2 <- reduce2grid(cgb)

scan.norm.imp <- scanone(gg_step2, method = "imp", model = "normal", pheno.col = 5)

scan.norm.imp <- scanone(gg_step2, method = "imp", model = "normal", pheno.col = 5)

sc2_normal_imp <- scantwo(gg_step2, pheno.col=5, model="normal",
             method="imp",
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs",
             incl.markers=TRUE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=6)

save.image(file.path(mpath,'scans.elr.rsave'))
################################################################################

################################################################################
# MQM ##########################################################################
crossaug <- mqmaugment(cross,strategy='drop')  # Augmentation

result <- mqmscan(crossaug)    # Scan
    # show LOD interval of the QTL on chr 3

mq.add <- mqmscan(crossaug, cofactors=NULL, pheno.col = 5,
      model=c("additive"), forceML=FALSE,
      cofactor.significance=0.02, em.iter=1000,
      window.size=25.0, step.size=5.0,
      logtransform = FALSE, estimate.map = FALSE,
      plot=FALSE, verbose=FALSE, outputmarkers=TRUE,
      multicore=TRUE, batchsize=10, n.clusters=6, test.normality=FALSE,off.end=0
      )

mq.dom <- mqmscan(crossaug, cofactors=NULL, pheno.col = 5,
      model="dominance", forceML=FALSE,
      cofactor.significance=0.02, em.iter=1000,
      window.size=25.0, step.size=5.0,
      logtransform = FALSE, estimate.map = FALSE,
      plot=FALSE, verbose=FALSE, outputmarkers=TRUE,
      multicore=TRUE, batchsize=10, n.clusters=6, test.normality=FALSE,off.end=0
      )

save.image(file.path(mpath,'scans.elr.rsave'))
################################################################################
