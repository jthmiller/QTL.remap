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
###############################################################################

file_list <- list.files(mpath, 'NBH_all_mark_?[0-9]?[0-9]_tsp.csv')

chr <- gsub("NBH_all_mark_",'',file_list)
chr <- as.numeric(gsub("_tsp.csv",'',chr))

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

write.table(to_write,file.path(mpath,'nbh.mapped.tsp.csv'),sep=',',row.names=F,quote=F,col.names = F)

################################################################################
################################################################################
fl <- file.path(mpath,'nbh.mapped.tsp.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

cross <- sim.geno(cross)
cross <- calc.genoprob(cross,step=1,error.prob=0.01,off.end=5)

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

##SEX
scan.bin.sex <- scanone(cross, method = "em", model = "binary", pheno.col = 2)
################################################################################

################################################################################
save.image(file.path(mpath,'single_scans.nbh.rsave'))
################################################################################

################################################################################
## step-wise
full.norm.add_only <- stepwiseqtl(cross, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=5)
##full.bin.add_only <- stepwiseqtl(cross, additive.only = T, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=5)
################################################################################

################################################################################
save.image(file.path(mpath,'single_scans.nbh.rsave'))
################################################################################

####################################################################################
## PERMS WITH ALL LOCI

perms.norm.imp <- scanone(cross, method = "imp", model = "normal", maxit = 1000,
  n.perm = 1000, pheno.col = 5, n.cluster = 10)

perms.bin.em <- scanone(cross, method = "em", model = "binary", maxit = 1000,
  n.perm = 1000, pheno.col = 4, n.cluster = 10)
################################################################################


####################################################################################

save.image(file.path(mpath,'single_scans.nbh.rsave'))

################################################################################
