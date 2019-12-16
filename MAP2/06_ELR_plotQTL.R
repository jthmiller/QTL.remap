#!/bin/R
### Map QTLs 1 of 3
#debug.cross <- T
#source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library('qtl')
pop <- 'ELR'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

##load(file.path(mpath,'scans.elr.rsave'))
save.image(file.path(mpath,'single_scans.elr.rsave'))

bins <- data.frame(
 em=summary(scan.bin.em),
 imp=summary(scan.bin.imp)[,'lod'],
 mr=summary(scan.bin.mr)[,'lod'],
 np=summary(scan.np.em.b)[,'lod'])

binpo <- data.frame(
 empo=rownames(summary(scan.bin.em)),
 impo=rownames(summary(scan.bin.imp)),
 mrpo=rownames(summary(scan.bin.mr)),
 nppo=rownames(summary(scan.np.em.b)))

norms <- data.frame(
 em=summary(scan.norm.em),
 imp=summary(scan.norm.imp)[,'lod'],
 mr=summary(scan.norm.mr)[,'lod'],
 np=summary(scan.np.em.n)[,'lod'],
 ehk=summary(scan.norm.ehk)[,'lod'])

normpo <- data.frame(
 empo=rownames(summary(scan.norm.em)),
 impo=rownames(summary(scan.norm.imp)),
 mrpo=rownames(summary(scan.norm.mr)),
 nppo=rownames(summary(scan.np.em.n)))

png(paste0('~/public_html/ELR_np.png'), width=1000)
 plot(scan.np.em.b)
dev.off()

png(paste0('~/public_html/ELR_rf_9.png'))
 plotRF(cross,chr=9)
dev.off()

png(paste0('~/public_html/ELR_scan.norm.mr.png'))
 plot(scan.norm.mr)
dev.off()


cbind(summary(out.ap13),summary(out.ap23))

save.image(file.path(mpath,'scans.elr.rsave'))
