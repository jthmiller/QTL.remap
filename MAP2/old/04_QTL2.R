#!/bin/R
### Map QTLs 1 of 3
#debug.cross <- T
#source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library('qtl')
pop <- 'NBH'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

################################################################################
## put chromosomes together
################################################################################

file_list <- list.files(mpath, 'NBH_gts_CHR.*_downsmpl_reordered.*')

chr <- gsub("NBH_gts_CHR",'',file_list)
chr <- as.numeric(gsub("_downsmpl_refined.csv",'',chr))

nbh <- lapply(file_list,function(X){ read.cross(file=X,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)})

gnos <- lapply(nbh,function(X){
 data.frame(X[[1]][[1]][['data']],stringsAsFactors=F)
})
gnos <- do.call(cbind,gnos)
gnos <- cbind(nbh[[1]]$pheno,gnos)
gnos$ID <- as.character(gnos$ID)

m_names <- unlist(sapply(nbh,function(X){
 markernames(X)
}))

colnames(gnos) <- c('Pheno','sex','ID','bin','pheno_norm',m_names)
rownames(gnos) <- nbh[[1]]$pheno$ID

map <- c(colnames(nbh[[1]]$pheno),unname(unlist(sapply(nbh,pull.map))))
chr <- c(colnames(nbh[[1]]$pheno),gsub(":.*","",m_names))
info <- c(colnames(nbh[[1]]$pheno),m_names)
headers <- rbind(info,chr,map)
colnames(headers) <- headers[1,]
headers[2:3,1:5] <- ''

headers.u <- unname(data.frame(headers,row.names=NULL,stringsAsFactors=FALSE))
gnos.u <- unname(data.frame(lapply(gnos, as.character),row.names=NULL,stringsAsFactors=FALSE))
colnames(headers.u) <- colnames(gnos.u) <- headers.u[1,]
to_write <- rbind(headers.u,gnos.u)

write.table(to_write,file.path(mpath,'nbh.mapped.1_24.csv'),sep=',',row.names=F,quote=F,col.names = F)
################################################################################

################################################################################
fl <- file.path(mpath,'nbh.mapped.1_24.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

cross <- sim.geno(cross)

png(paste0('~/public_html/NBH_RF_remap.png'))
plotRF(cross)
dev.off()

## binary
scan.bin.em <- scanone(cross, method = "em", model = "binary", pheno.col = 4)
scan.bin.imp <- scanone(cross, method = "imp", model = "binary", pheno.col = 4)
scan.bin.mr <- scanone(cross, method = "mr", model = "binary", pheno.col = 4)

## normal
scan.norm.em <- scanone(cross, method = "em", model = "normal", pheno.col = 5)
scan.norm.mr <- scanone(cross, method = "mr", model = "normal", pheno.col = 5)
scan.norm.imp <- scanone(cross, method = "imp", model = "normal", pheno.col = 5)
scan.norm.ehk <- scanone(cross, method = "ehk", model = "normal", maxit = 5000, pheno.col = 5)

## non-parametric
scan.np.em.b <- scanone(cross, method = "em", model = "np", pheno.col = 4, maxit = 5000)
scan.np.em.n <- scanone(cross, method = "em", model = "np", pheno.col = 5, maxit = 5000)
save.image(file.path(mpath,'scans.nbh.rsave'))
################################################################################

################################################################################
## step-wise
full.norm.add_only <- stepwiseqtl(cross, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=8)
full.bin.add_only <- stepwiseqtl(cross, additive.only = T, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=8)
save.image(file.path(mpath,'scans.nbh.rsave'))
### reduce2grid(cross)
################################################################################

qtl.2 <- makeqtl(cross, chr=c(2), pos=c(12.436),  what="draws")
out.a.2 <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2, method="imp",model="binary",pheno.col=4)
out.i.2 <- addqtl(cross, qtl=qtl.2, formula=y~Q1*Q2, method="imp",model="binary",pheno.col=4)
## indicates possible interaction on 13 and 18

qtl18 <- makeqtl(cross, chr=c(18), pos=c(29.396),  what="draws")
out.a.18 <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2, method="imp",model="binary",pheno.col=4)
out.i.18 <- addqtl(cross, qtl=qtl, formula=y~Q1*Q2, method="imp",model="binary",pheno.col=4)

## manual
qtl <- makeqtl(cross, chr=c(2,8,18), pos=c(12.436,19.031,29.396),  what="draws")
fitted <- fitqtl(cross,qtl=qtl, formula=y~Q1+Q2+Q3)

out.a <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2, method="imp",model="binary",pheno.col=4)
out.ia <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2+Q1:Q3, method="imp",model="binary",pheno.col=4)

out.ia <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2+Q1:Q3, method="imp",model="binary",pheno.col=4)
out.ib <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2+Q2:Q3, method="imp",model="binary",pheno.col=4)

qtl <- makeqtl(cross, chr=c(13), pos=c(119.68569),  what="draws")
out.i.18 <- addqtl(cross, qtl=qtl, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=5)
out.a.18 <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=5)
####################################################################################

################################################################################
## PERMS
perms.bin.em <- scanone(cross, method = "em", model = "binary", maxit = 1000,
  n.perm = 1000, pheno.col = 4, n.cluster = 6)
perms.norm.em <- scanone(cross, method = "em", model = "normal", maxit = 1000,
  n.perm = 1000, pheno.col = 5, n.cluster = 6)
save.image(file.path(mpath,'scans.nbh.rsave'))
################################################################################

################################################################################
### LONG #######################################################################
full.norm <- stepwiseqtl(cross, additive.only = F, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=6)
save.image(file.path(mpath,'scans.nbh.rsave'))
full.bin <- stepwiseqtl(cross, additive.only = F, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=6)
save.image(file.path(mpath,'scans.nbh.rsave'))
full.norm.int <- stepwiseqtl(cross, additive.only = F, model='normal', method = "hk", pheno.col = 5, scan.pairs = T, max.qtl=6)
save.image(file.path(mpath,'scans.nbh.rsave'))
full.bin.int <- stepwiseqtl(cross, additive.only = F, model='binary', method = "hk", pheno.col = 4, scan.pairs = T, max.qtl=6)
save.image(file.path(mpath,'scans.nbh.rsave'))
################################################################################

################################################################################
sc2_normal_mr <- scantwo(cross, pheno.col=5, model="normal",
             method="mr",
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs",
             incl.markers=TRUE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=6)

sc2_binary_mr <- scantwo(cross, pheno.col=4, model="binary",
             method="mr",
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs",
             incl.markers=TRUE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=6)

save.image(file.path(mpath,'scans.nbh.rsave'))
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

save.image(file.path(mpath,'scans.nbh.rsave'))
################################################################################
