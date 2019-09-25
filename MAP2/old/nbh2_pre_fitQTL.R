#!/bin/R

pop <- 'NBH'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")

library('qtl')

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]

################################################################################
## put chromosomes together
################################################################################

file_list <- gsub('.csv','',list.files(mpath, 'NBH*downsmpl_map*'))
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


#####################################################################
