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

file_list <- list.files(mpath, 'ELR.*downsmpl_map.*')

file_list <- file_list[!file_list %in% c("ELR_gts_CHR11_downsmpl_map.csv","ELR_gts_CHR14_downsmpl_map.csv")]

chr <- gsub("ELR_gts_CHR",'',file_list)
chr <- as.numeric(gsub("_downsmpl_map.csv",'',chr))

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

png(paste0('~/public_html/ELR_RF_remap.png'))
plotRF(cross)
dev.off()


################################################################################
#### ADD IN AHR GENOS

#AHR2a_del <- 343745
#AIP_261 <- 29370504
#AIP_252 <- 29370500

#pull.map(cross)[[1]]["1:363497"]
#pull.map(cross)[[2]]["2:29667555"]

fla <-file.path(mpath, 'ER_ahr_aip_whoi_gt.csv')
cross.ahr <- read.cross(file = fla,format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

ahr2 <- pull.geno(cross.ahr)[,"AHR2a_del"]
aip252 <- pull.geno(cross.ahr)[,"AIP_252"]
aip261 <- pull.geno(cross.ahr)[,"AIP_261"]

add_gts <- data.frame(ahr2,aip252,aip261,stringsAsFactors=F)
rownames(add_gts) <- cross.ahr$pheno$ID
add_gts <- add_gts[as.character(cross$pheno$ID),]

cross <- addmarker(cross,add_gts[,'ahr2'],'ahr2a',chr=1,pos=2.643695)
cross <- addmarker(cross,add_gts[,'aip261'],'aip261',chr=2,pos=6.5)
cross <- addmarker(cross,add_gts[,'aip252'],'aip252',chr=2,pos=6.6)

##i <- 1
##cr1 <- subset(cross,chr=1)
##cr1s <- as.character(cr1$pheno$ID)[order(as.numeric(cr1$pheno$Pheno))]
##
##png(paste0('~/public_html/ELR_gts_pheno_order',i,'.png'),height=2500,width=4500)
##plotGeno(cr1, idv=cr1s)
##dev.off()
##
##n_ahr <- pull.geno(cross)[,"1:363497"]
##cbind(add_gts,n_ahr,as.character(cross$pheno$ID))
##cbind(cr1s,cr1$pheno)
################################################################################

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

## step-wise
full.norm.add_only <- stepwiseqtl(cross, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=5)
full.bin.add_only <- stepwiseqtl(cross, additive.only = T, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=5)

save.image(file.path(mpath,'scans.elr.rsave'))

## PERMS
perms.bin.em <- scanone(cross, method = "em", model = "binary", maxit = 1000,
  n.perm = 1000, pheno.col = 4, n.cluster = 6)
perms.norm.em <- scanone(cross, method = "em", model = "normal", maxit = 1000,
  n.perm = 1000, pheno.col = 5, n.cluster = 6)

save.image(file.path(mpath,'scans.elr.rsave'))

full.norm <- stepwiseqtl(cross, additive.only = F, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=6)

save.image(file.path(mpath,'scans.elr.rsave'))

full.bin <- stepwiseqtl(cross, additive.only = F, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=6)

save.image(file.path(mpath,'scans.elr.rsave'))

full.norm.int <- stepwiseqtl(cross, additive.only = F, model='normal', method = "hk", pheno.col = 5, scan.pairs = T, max.qtl=6)

save.image(file.path(mpath,'scans.elr.rsave'))

full.bin.int <- stepwiseqtl(cross, additive.only = F, model='binary', method = "hk", pheno.col = 4, scan.pairs = T, max.qtl=6)

save.image(file.path(mpath,'scans.elr.rsave'))


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

save.image(file.path(mpath,'scans.elr.rsave'))

# MQM
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
