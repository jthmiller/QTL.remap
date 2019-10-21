#!/bin/R
### Map QTLs 1 of 3
#debug.cross <- T
#source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library('qtl')
pop <- 'ELR'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

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
zd <- as.numeric(gsub(":.*","",m_names))

zd[is.na(zd)] <- c(1,2,2)
chr <- c(colnames(elr[[1]]$pheno),zd)
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
## scan
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

## normal
scan.norm.em <- scanone(cross, method = "em", model = "normal", pheno.col = 1)
scan.norm.mr <- scanone(cross, method = "mr", model = "normal", pheno.col = 1)
scan.norm.imp <- scanone(cross, method = "imp", model = "normal", pheno.col = 1)
scan.norm.ehk <- scanone(cross, method = "ehk", model = "normal", maxit = 5000, pheno.col = 1)

## normal transform
scan.normT.em <- scanone(cross, method = "em", model = "normal", pheno.col = 5)
scan.normT.mr <- scanone(cross, method = "mr", model = "normal", pheno.col = 5)
scan.normT.imp <- scanone(cross, method = "imp", model = "normal", pheno.col = 5)
scan.normT.ehk <- scanone(cross, method = "ehk", model = "normal", maxit = 5000, pheno.col = 5)

## non-parametric
scan.np.em.b <- scanone(cross, method = "em", model = "np", pheno.col = 4, maxit = 5000)
scan.np.em.n <- scanone(cross, method = "em", model = "np", pheno.col = 5, maxit = 5000)

##SEX
scan.bin.sex <- scanone(cross, method = "em", model = "binary", pheno.col = 2)
################################################################################

save.image(file.path(mpath,'single_scans.elr.rsave'))

################################################################################
## step-wise
full.norm.add_only <- stepwiseqtl(cross, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=4)
##full.bin.add_only <- stepwiseqtl(cross, additive.only = T, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=3)
################################################################################

save.image(file.path(mpath,'single_scans.elr.rsave'))

################################################################################
## PERMS WITH ALL LOCI
perms.norm.imp <- scanone(cross, method = "imp", model = "normal", maxit = 10000,
  n.perm = 1000, pheno.col = 5, n.cluster = 10)

perms.bin.em <- scanone(cross, method = "em", model = "binary", maxit = 10000,
  n.perm = 1000, pheno.col = 4, n.cluster = 10)
################################################################################

save.image(file.path(mpath,'single_scans.elr.rsave'))
