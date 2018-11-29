#!/bin/R
### first run combine pops for multi-pop cross objects

debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library("ggridges")
library('plyr')
library('scales')
library('ggrepel')
#### AHRs #####
AHR.bed <- read.table("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/data/lift_AHR_genes.bed",
  stringsAsFactors = F, header = F)
colnames(AHR.bed) <- c("chrom", "str", "stp", "gene")
AHR.bed$chrom <- as.numeric(gsub("chr", "", AHR.bed$chrom))
AHR.bed$str <- as.numeric(AHR.bed$str)
AHR.bed$stp <- as.numeric(AHR.bed$stp)
AHR.notmap <- AHR.bed[is.na(AHR.bed$chrom), ]
AHR.bed <- AHR.bed[!is.na(AHR.bed$chrom), ]
# add arnts (forgot to scan for them)

#######
popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NBH/REMAPS"
cross.NBH <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv",
  sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))

popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/ELR/REMAPS"
cross.ELR <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv",
  sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))

popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NEW/REMAPS"
cross.NEW <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv",
  sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))
#################
mak <- markernames(cross.NEW)
cross.NEW <- switchAlleles(cross.NEW, markers=mak)
#################
#########
cross.nbh <- sim.geno(cross.NBH, n.draws = 500, step = 5, off.end = 10, error.prob = 0.01,
  map.function = "kosambi", stepwidth = "fixed")
cross.new <- sim.geno(cross.NEW, n.draws = 500, step = 5, off.end = 10, error.prob = 0.01,
  map.function = "kosambi", stepwidth = "fixed")
cross.elr <- sim.geno(cross.ELR, n.draws = 500, step = 5, off.end = 10, error.prob = 0.05,
  map.function = "kosambi", stepwidth = "fixed")
#################
plotInfo(hyper)


sex <- read.table(file = file.path(dirso, "sex.txt"))
rownames(sex) <- sex$ID
cross.nbh$pheno$sex <- sex[as.character(cross.nbh$pheno$ID), 2]
cross.nbh$pheno$binary <- as.numeric(cross.nbh$pheno$pheno >= 3)
cross.new$pheno$sex <- sex[as.character(cross.new$pheno$ID), 2]
cross.new$pheno$binary <- as.numeric(cross.new$pheno$pheno >= 3)
cross.elr$pheno$sex <- sex[as.character(cross.elr$pheno$ID), 2]
cross.elr$pheno$binary <- as.numeric(cross.elr$pheno$pheno >= 3)

cross.nbh <- reduce2grid(cross.nbh)
cross.new <- reduce2grid(cross.new)
cross.elr <- reduce2grid(cross.elr)

#scan.norm.imp.NBH <- scanone(cross.nbh, model = "normal", pheno.col = 1, method = "imp", addcovar = cross.nbh$pheno$sex)
scan.norm.imp.NBH <- scanone(cross.nbh, model = "normal", pheno.col = 1, method = "imp")
#scan.norm.imp.NEW <- scanone(cross.new, model = "normal", pheno.col = 1, method = "imp", addcovar = cross.new$pheno$sex)
scan.norm.imp.NEW <- scanone(cross.new, model = "normal", pheno.col = 1, method = "imp")
#scan.norm.imp.ELR <- scanone(cross.elr, model = "normal", pheno.col = 1, method = "imp", addcovar = cross.elr$pheno$sex)
scan.norm.imp.ELR <- scanone(cross.elr, model = "normal", pheno.col = 1, method = "imp")
### use scanone for plots
themelt.nbh <- scan.norm.imp.NBH
themelt.new <- scan.norm.imp.NEW
themelt.elr <- scan.norm.imp.ELR
themelt.nbh$pop <- 'NBH'
themelt.new$pop <- 'NEW'
themelt.elr$pop <- 'ELR'

##save.image('/home/jmiller1/public_html/QTL_plot.Rsave')

### ggplot format AHR genes
cnv.ahrs <- function(cross2,AHRdf,EXP){
  mapo <- convert2cross2(cross2)
  mapo$pmap <- mapo$gmap

  mapo$pmap <- lapply(mapo$pmap,function(X){
    return(as.numeric(gsub('[0-9]+:','',names(X))))
    }
  )
  for (i in 1:24){
    names(mapo$pmap[[i]]) <- names(mapo$gmap[[i]])
  }
  AHR.list <- split(AHRdf,AHRdf$chrom)
  AHR.pos1 <- lapply(AHR.list,'[[',2)
  AHR.pos2 <- lapply(AHR.list,'[[',3)
  AHR.pos1 <- lapply(AHR.pos1,as.numeric)
  AHR.pos2 <- lapply(AHR.pos2,as.numeric)
  pos1 <- interp_map(AHR.pos1,mapo$pmap,mapo$gmap)
  pos2 <- interp_map(AHR.pos2,mapo$pmap,mapo$gmap)
  AHR.list <- mapply(cbind,AHR.list,pos1=pos1,SIMPLIFY=FALSE)
  AHR.list <- mapply(cbind,AHR.list,pos=pos2,SIMPLIFY=FALSE)
  ah.gens <- ldply(AHR.list, data.frame, .id = NULL)
  ah.gens$chrom <- factor(ah.gens$chrom,levels=as.character(1:24))
  colnames(ah.gens)[1] <- 'chr'
  #ah.gens$lod <- rep_len(c(-1:-10), length(ah.gens[,1]))
  ah.gens$lod <- 0
  if(EXP==F){ ah.gens <- ah.gens[which(!ah.gens$gene=='EXPRESSED'),]}
  ah.gens <- ah.gens[!grepl('*many*',ah.gens$gene),]
  ah.gens <- ah.gens[!grepl('*many*',ah.gens$gene),]
  return(ah.gens)
}
nbh.gens <- cnv.ahrs(cross2=cross.nbh,AHRdf=AHR.bed,EXP=F)
new.gens <- cnv.ahrs(cross.new,AHRdf=AHR.bed,EXP=F)
elr.gens <- cnv.ahrs(cross.elr,AHRdf=AHR.bed,EXP=F)
qtl.gens <-nbh.gens[which(nbh.gens$chr %in% c(1,2,5,8,10,12,13,18,24)),]
minor.gens  <- nbh.gens[which(nbh.gens$chr %in% c(8,13,23,24)),]

### ggplot popgen locations
dir <- '/home/jmiller1/QTL_Map_Raw/popgen/tables'
nbh.popgen <- read.table(file.path(dir,'outliersNBH.txt.ncbi.lifted'),sep='\t',header=T)
new.popgen <- read.table(file.path(dir,'outliersNYC.txt.ncbi.lifted'),sep='\t',header=T)
elr.popgen <- read.table(file.path(dir,'outliersER.txt.ncbi.lifted'),sep='\t',header=T)

### Use nbh coords but elr and new popgen
cnv.popgen <- function(cross2,popgen,top){
  nbhmap <- convert2cross2(cross2)
  nbhmap$pmap <- nbhmap$gmap
  nbhmap$pmap <- lapply(nbhmap$pmap,function(X){
    return(as.numeric(gsub('[0-9]+:','',names(X))))
    }
  )
  for (i in 1:24){
    names(nbhmap$pmap[[i]]) <- names(nbhmap$gmap[[i]])
  }

  popgen$chrom <- gsub('chr','',popgen$chrom)
  colnames(popgen)[1] <- 'chr'
  pop.list <- split(popgen,popgen$chr)
  pop.list <- pop.list[as.character(c(1:24))]
  popgen.pos1 <- lapply(pop.list,'[[',2)
  popgen.pos2 <- lapply(pop.list,'[[',3)
  popgen.pos1 <- lapply(popgen.pos1,as.numeric)
  popgen.pos2 <- lapply(popgen.pos2,as.numeric)
  pos1 <- interp_map(popgen.pos1,nbhmap$pmap,nbhmap$gmap)
  pos2 <- interp_map(popgen.pos2,nbhmap$pmap,nbhmap$gmap)
  pop.list <- mapply(cbind,pop.list,pos1=pos1,SIMPLIFY=FALSE)
  pop.list <- mapply(cbind,pop.list,pos=pos2,SIMPLIFY=FALSE)
  fraps <- ldply(pop.list, data.frame, .id = NULL)
  fraps$chr <- factor(fraps$chr,levels=as.character(1:24))
  fraps <- fraps[which(fraps$rank<top),]
  fraps$lod <- rescale(fraps$rank, to = c(20,0))
  fraps$sz <- rescale(fraps$rank, to = c(4,1))
  return(fraps)
}
new.rank <- cnv.popgen(cross.nbh,new.popgen,top=50)
nbh.rank <- cnv.popgen(cross.nbh,nbh.popgen,top=50)
elr.rank <- cnv.popgen(cross.nbh,elr.popgen,top=50)
nbh.rank$pop <- 'NBH'
new.rank$pop <- 'NEW'
elr.rank$pop <- 'ELR'
all.rank <- rbind(new.rank,nbh.rank,elr.rank)
all.rank$pop <- factor(all.rank$pop, levels = c("NBH", "NEW", "ELR"))
qtl.rank <- all.rank[which(all.rank$chr %in% c(1,2,5,8,10,12,13,18,24)),]
minor.rank <- all.rank[which(all.rank$chr %in% c(8,13,23,24)),]

### GGriges plot
melted.nbh <- data.frame(pop = "NBH", chr = scan.norm.imp.NBH$chr, pos = scan.norm.imp.NBH$pos,
  lod = scan.norm.imp.NBH$lod)
melted.new <- data.frame(pop = "NEW", chr = scan.norm.imp.NEW$chr, pos = scan.norm.imp.NEW$pos,
  lod = scan.norm.imp.NEW$lod)
melted.elr <- data.frame(pop = "ELR", chr = scan.norm.imp.ELR$chr, pos = scan.norm.imp.ELR$pos,
  lod = scan.norm.imp.ELR$lod)

melted <- rbind(melted.nbh, melted.new, melted.elr)
melted$pop <- factor(melted$pop, levels = rev(c("NBH", "NEW", "ELR")))

##Total CM length of NBH. Rescale to NBH
mxes <- sapply(1:24,function(X){max(themelt.nbh$pos[which(themelt.nbh$chr==X)])})

ts <- themelt.new[which(themelt.new$chr==1),]
ts$pos <- rescale(ts$pos, to=c(-10,mxes[1]))
new.rescale <- ts
for(i in 2:24){
 ts <- themelt.new[which(themelt.new$chr==i),]
 ts$pos <- rescale(ts$pos, to=c(-10,mxes[i]))
 new.rescale <- rbind(new.rescale,ts)
}

ts <- themelt.elr[which(themelt.elr$chr==1),]
ts$pos <- rescale(ts$pos, to=c(-10,mxes[1]))
elr.rescale <- ts
for(i in 2:24){
 ts <- themelt.elr[which(themelt.elr$chr==i),]
 ts$pos <- rescale(ts$pos, to=c(-10,mxes[i]))
 elr.rescale <- rbind(elr.rescale,ts)
}

allmelt <- rbind(themelt.nbh,new.rescale,elr.rescale)
allmelt$pop <- factor(allmelt$pop, levels = c("NBH", "NEW", "ELR"))
qtlmelt <- allmelt[which(allmelt$chr %in% c(1,2,5,8,10,12,13,18,24)),]
qtlminor  <- allmelt[which(allmelt$chr %in% c(8,13,23,24)),]
#### Subest for only qtl plots
##qtl.rank <- nbh.rank[which(nbh.rank$chr %in% c(1,2,5,8,10,12,13,18,24)),]

#save.image('/home/jmiller1')

# All chr Not scaled either
png("/home/jmiller1/public_html/ridges_accurate map_length.qtl.png", width = 3000)
ggplot(melted, aes(y = as.factor(pop), x = pos, height = lod, fill = pop)) +
  geom_ridgeline(stat = "identity",alpha = 0.5, scale = 0.25) +
  facet_wrap(~as.factor(chr), scales = "free_x", nrow = 1,ncol = 24, strip.position = NULL) +
  # scale_x_continuous('pos') + scale_y_continuous('lod') +
  scale_y_discrete("lod") + theme(axis.text = element_text(size = 10))
dev.off()

p <- ggplot(themelt.nbh, aes(x = pos, y = lod))
png("/home/jmiller1/public_html/all_popgen_rank_scaled.qtl.png", width = 3000)
  p + facet_wrap(~chr, scales = "free_x",nrow = 1,ncol = 24) +
  scale_y_continuous(limits = c(-12, 23)) +
  geom_line(size = 2, alpha = 0.6) +
  theme_minimal() +
  theme(axis.text = element_text(size = 10)) +
  labs(x = "Chromosome",y = "LOD",linetype = "") +
  geom_label_repel(aes(x = pos, y=-0.1,label = gene),
    box.padding = unit(0.25, "lines") , parse = T,point.padding =unit(0.2, "lines"),force = 10,
    label.padding= unit(0.2, "lines"),ylim = c(0,-12),  segment.size =1,max.iter = 6000,
    data = nbh.gens, direction = 'y', size = 4, seed=666,
    nudge_y = -0.01, vjust=3) +
  geom_label_repel(aes(x = pos, y=-0.1,label = gene),
    box.padding = unit(0.25, "lines") , parse = T,point.padding =unit(0.2, "lines"),force = 10,
    label.padding= unit(0.2, "lines"),ylim = c(0,-12),  segment.size =0,max.iter = 6000,
    data = nbh.gens, direction = 'y', size = 4, seed=666,
    nudge_y = -0.01, vjust=3) +
  geom_point(aes(size = rank), data = nbh.rank) +
  geom_label_repel(aes(x = pos, y = lod,label = rank),
    data = nbh.rank, size = 4, box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"), vjust=1)
dev.off()

p <- ggplot(themelt.new, aes(x = pos, y = lod))
png("/home/jmiller1/public_html/new_genes_below_rank_above.qtl.png", width = 3000)
  p + facet_wrap(~chr, scales = "free_x",nrow = 1,ncol = 24) +
  scale_y_continuous(limits = c(-12, 23)) +
  geom_line(size = 2, alpha = 0.6) +
  theme_minimal() +
  theme(axis.text = element_text(size = 10)) +
  labs(x = "Chromosome",y = "LOD",linetype = "") +
  geom_label_repel(aes(x = pos, y=-0.1,label = gene),
    box.padding = unit(0.25, "lines") , parse = T,point.padding =unit(0.2, "lines"),force = 10,
    label.padding= unit(0.2, "lines"),ylim = c(0,-12),  segment.size =1,max.iter = 6000,
    data = new.gens, direction = 'y', size = 4, seed=666,
    nudge_y = -0.01, vjust=3) +
  geom_label_repel(aes(x = pos, y=-0.1,label = gene),
    box.padding = unit(0.25, "lines") , parse = T,point.padding =unit(0.2, "lines"),force = 10,
    label.padding= unit(0.2, "lines"),ylim = c(0,-12),  segment.size =0,max.iter = 6000,
    data = new.gens, direction = 'y', size = 4, seed=666,
    nudge_y = -0.01, vjust=3) +
  geom_label_repel(aes(x = pos, y = lod,label = rank),
    data = nbh.rank, size = 4, box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"), vjust=1)
dev.off()

p <- ggplot(themelt.elr, aes(x = pos, y = lod,color=pop))
png("/home/jmiller1/public_html/rank_top_scaled_random.qtl.png", width = 3000)
  p + facet_wrap(~chr, scales = "free_x",nrow = 1,ncol = 24) +
  scale_y_continuous(limits = c(-12, 23)) +
  scale_color_manual(values=popcol)+
  geom_line(size = 2, alpha = 0.6) +
  theme_minimal() +
  theme(axis.text = element_text(size = 10)) +
  labs(x = "Chromosome",y = "LOD",linetype = "") +
  geom_label_repel(aes(x = pos, y=-0.1,label = gene),
    box.padding = unit(0.25, "lines") , parse = T,point.padding =unit(0.2, "lines"),force = 10,
    label.padding= unit(0.2, "lines"),ylim = c(0,-12),  segment.size =1,max.iter = 6000,
    data = elr.gens, direction = 'y', size = 4, seed=666,
    nudge_y = -0.01, vjust=3) +
  geom_label_repel(aes(x = pos, y=-0.1,label = gene),
    box.padding = unit(0.25, "lines") , parse = T,point.padding =unit(0.2, "lines"),force = 10,
    label.padding= unit(0.2, "lines"),ylim = c(0,-12),  segment.size =0,max.iter = 6000,
    data = elr.gens, direction = 'y', size = 4, seed=666,
    nudge_y = -0.01, vjust=3) +
  geom_label_repel(aes(x = pos, y = lod,label = rank),
    data = elr.rank, size = 4, box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"), vjust=1)
dev.off()

p <- ggplot(allmelt, aes(x = pos, y = lod, colour=pop))
png("/home/jmiller1/public_html/no_annot_.qtl.png", width = 2000)
  p + facet_wrap(~chr,nrow = 1,scales = "free_x",ncol = 24) +
  scale_y_continuous(limits = c(0, 22)) +
  scale_color_manual(values=popcol)+
  geom_line(size = 2, alpha = 0.6) +
  theme_minimal() +
  theme(axis.text = element_text(size = 10)) +
  labs(x = "Chromosome",y = "LOD",linetype = "")
dev.off()

p <- ggplot(allmelt, aes(x = pos, y = lod, colour=pop))
png("/home/jmiller1/public_html/ahr_genes_only.qtl.png", width = 2000)
  p + facet_wrap(~chr,nrow = 1,scales = "free_x",ncol = 24) +
  scale_y_continuous(limits = c(0, 22)) +
  scale_color_manual(values=popcol)+
  geom_line(size = 2, alpha = 0.6) +
  theme_minimal() +
  theme(axis.text = element_text(size = 10)) +
  geom_label_repel(aes(x = pos, y=0.1,label = gene),
    box.padding = unit(0.5, "lines") , parse = T,point.padding =unit(0.8, "lines"),force = 1,
    label.padding= unit(0.2, "lines"),ylim = c(5,20),  segment.size =1,max.iter = 6000,
    data = nbh.gens, direction = 'y', size = 3, seed=666,
    nudge_y = 0.01, vjust=.1) +
  geom_label_repel(aes(x = pos, y=0.1,label = gene),
    box.padding = unit(0.5, "lines") , parse = T,point.padding =unit(0.8, "lines"),force = 1,
    label.padding= unit(0.2, "lines"),ylim = c(5,20),  segment.size =0,max.iter = 6000,
    data = nbh.gens, direction = 'y', size = 3, seed=666,
    nudge_y = 0.01, vjust=.1) +
  labs(x = "Chromosome",y = "LOD",linetype = "")
dev.off()

p <- ggplot(allmelt, aes(x = pos, y = lod, colour =pop))
png("/home/jmiller1/public_html/All_chr.qtl.png", width = 2000)
  p + facet_wrap(~chr,nrow = 1,scales = "free_x",ncol = 24) +
  scale_color_manual(values=popcol)+
  scale_y_continuous(limits = c(-5, 22)) +
  theme_minimal() +
    geom_label_repel(aes(x = pos, y=0.1,label = gene),color='black',segment.size =1,label.padding= unit(0.2, "lines"),
    data = nbh.gens, direction = 'y', size = 3, seed=666,ylim = c(5,20),segment.alpha=.5) +
  geom_line(size = 2, alpha = 0.5) +
  geom_label_repel(aes(x = pos, y=0.1,label = gene),color='black',label.padding= unit(0.2, "lines"),
    ylim = c(5,20),  segment.size =0, max.iter = 6000,nudge_y=2,
    data = nbh.gens, direction = 'y', size = 3.5, seed=666) +
  geom_label_repel(aes(x = pos, y = 0,fill=pop,label = rank),
    data = all.rank, size = 3.5,segment.size =1,force=4,min.segment.length=0.1,
    point.padding = unit(0.4, "lines"),direction = 'both',ylim = c(0,-8),seed=666,box.padding=0.1)+
  geom_label_repel(aes(x = pos, y = 0,label = rank,), fontface = "bold",
    data = all.rank, size = 3.5, segment.size =0,force=4,min.segment.length=0.1,
    point.padding = unit(0.4, "lines"),direction = 'both',ylim = c(0,-8),seed=666,box.padding=0.1)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")+
  labs(x = "Chromosome",y = "LOD",linetype = "")
dev.off()

p <- ggplot(qtlmelt, aes(x = pos, y = lod, colour =pop))
png("/home/jmiller1/public_html/all_pop_qtl_only.png", width = 2000)
  p + facet_wrap(~chr,nrow = 1,scales = "free_x",ncol = 9) +
  scale_color_manual(values=popcol)+
  scale_y_continuous(limits = c(-5, 22)) +
  theme_minimal() +
    geom_label_repel(aes(x = pos, y=0.1,label = gene),color='black',segment.size =1,label.padding= unit(0.2, "lines"),
    data = qtl.gens , direction = 'y', size = 5, seed=666,nudge_y=2,ylim = c(5,20),segment.alpha=.5) +
  geom_line(size = 2, alpha = 0.5) +
  geom_label_repel(aes(x = pos, y=0.1,label = gene),color='black',label.padding= unit(0.2, "lines"),
    ylim = c(5,20),  segment.size =0, max.iter = 6000,nudge_y=2,
    data = qtl.gens , direction = 'y', size = 5, seed=666) +
  geom_label_repel(aes(x = pos, y = 0,color=pop,label = rank),
    data = qtl.rank, size = 5,segment.size =1,force=4,min.segment.length=0.1,
    point.padding = unit(0.4, "lines"),direction = 'both',ylim = c(0,-8),seed=666,box.padding=0.1)+
  geom_label_repel(aes(x = pos, y = 0,label = rank), fontface = "bold",
    data = qtl.rank, size = 5, segment.size =0,force=4,min.segment.length=0.1,
    point.padding = unit(0.4, "lines"),direction = 'both',ylim = c(0,-8),seed=666,box.padding=0.1)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")+
  labs(x = "Chromosome",y = "LOD",linetype = "")
dev.off()
## Why would chromosomes with top ranked loci have fewest signatures?
##

p <- ggplot(qtlminor, aes(x = pos, y = lod, colour =pop))
png("/home/jmiller1/public_html/minor_qtl_only.png", width = 2000)
  p + facet_wrap(~chr,nrow = 1,scales = "free_x",ncol = 9) +
  scale_color_manual(values=popcol)+
  scale_y_continuous(limits = c(-5, 22)) +
  theme_minimal() +
    geom_label_repel(aes(x = pos, y=0.1,label = gene),color='black',segment.size =1,label.padding= unit(0.2, "lines"),
    data = minor.gens , direction = 'y', size = 5, seed=666,nudge_y=2,ylim = c(5,20),segment.alpha=.5) +
  geom_line(size = 2, alpha = 0.5) +
  geom_label_repel(aes(x = pos, y=0.1,label = gene),color='black',label.padding= unit(0.2, "lines"),
    ylim = c(5,20),  segment.size =0, max.iter = 6000,nudge_y=2,
    data = minor.gens , direction = 'y', size = 5, seed=666) +
  geom_label_repel(aes(x = pos, y = 0,color=pop,label = rank),
    data = minor.rank, size = 5,segment.size =1,force=4,min.segment.length=0.1,
    point.padding = unit(0.4, "lines"),direction = 'both',ylim = c(0,-8),seed=666,box.padding=0.1)+
  geom_label_repel(aes(x = pos, y = 0,label = rank), fontface = "bold",
    data = minor.rank, size = 5, segment.size =0,force=4,min.segment.length=0.1,
    point.padding = unit(0.4, "lines"),direction = 'both',ylim = c(0,-8),seed=666,box.padding=0.1)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")+
  labs(x = "Chromosome",y = "LOD",linetype = "")
dev.off()










rabbit hole
i <- 1
disto <- sapply(1:24,function(X){
length(which(all.rank$chr==X & all.rank$rank<(max(all.rank$rank)*i)))
})

for (i in seq.int(0.01,0.99,.025)){

  bindo <- sapply(1:24,function(X){
  length(which(all.rank$chr==X & all.rank$rank<(max(all.rank$rank)*i)))
})
  disto <- rbind(disto,bindo)
}

disto <- t(disto)
colnames(disto) <- c(1,seq.int(0.01,0.99,.025))
rownames(disto) <- 1:24

png('/home/jmiller1/public_html/rabithole.png')
plot(as.numeric(colnames(disto)),disto[1,])
sapply(2:24, function(X){ points(as.numeric(colnames(disto)),disto[X,])})
dev.off()
































plotfig <- function(gens,melts){

p <- ggplot(melts, aes(x = pos, y = lod)) +
  #geom_hline(data = lod, aes(yintercept = lod, linetype = group)) +
  facet_wrap(~chr, scales = "free_x",nrow = 1,ncol = 24) +
  scale_y_continuous(limits = c(-12, 23)) +
  geom_line(size = 2, alpha = 0.6) +
  theme_minimal() +
  theme(axis.text = element_text(size = 10)) +
  labs(x = "Chromosome",y = "LOD",linetype = "") +
  geom_label_repel(aes(x = pos, y=-0.1,label = gene),
    box.padding = unit(0.25, "lines") , parse = T,point.padding =unit(0.2, "lines"),force = 10,
    label.padding= unit(0.2, "lines"),ylim = c(0,-12),  segment.size =1,max.iter = 6000,
    data = gens, direction = 'y', size = 4, seed=666,
    nudge_y = -0.01, vjust=3) +
  geom_label_repel(aes(x = pos, y=-0.1,label = gene),
    box.padding = unit(0.25, "lines") , parse = T,point.padding =unit(0.2, "lines"),force = 10,
    label.padding= unit(0.2, "lines"),ylim = c(0,-12),  segment.size =0,max.iter = 6000,
    data = gens, direction = 'y', size = 4, seed=666,
    nudge_y = -0.01, vjust=3) +
  geom_label_repel(aes(x = pos, y = lod,label = rank),
    data = nbh.rank, size = 4, box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"), vjust=1)
}


plotfig(gens=nbh.gens,melts=themelt.nbh)

plotfig(new.gens, themelt.new)

a <- plotfig(elr.gens, themelt.elr)
png("/home/jmiller1/public_html/ggplot2.qtl.png", width = 3000)
a
dev.off()




size = 4, box.padding = 2,direction = 'y',
fill = 'grey',ylim = c(-8,-0),force=10,segment.size=0.2


,direction ='y' force=10,
point.padding = unit(0.3, "lines"),box.padding = unit(0.35, "lines"),vjust=1,direction = 'y',

png("/home/jmiller1/public_html/ggplot2.qtl.png", width = 3000)
  ggplot(themelt, aes(x = pos, y = lod)) +
  #geom_hline(data = lod, aes(yintercept = lod, linetype = group)) +
  facet_wrap(~ chr, scales = "free_x",nrow = 1,ncol = 24) +
  geom_line(size = 2, alpha = 0.6) +
  theme_minimal() +
  theme(axis.text = element_text(size = 10)) +
  labs(x = "Chromosome",y = "LOD",color = "",linetype = "") +
  geom_text(aes(x = pos, y = lod,label = gene,group=NULL),
    data = nbh.gens)
dev.off()

!!!! Check the factor levels in this!!

load(file.path('~/NEW.NBH.ELR.QTLmap.Rsave'))
# load(file.path('~/NEW.NBH.ELR.QTLmap.Rsave'))

## Lod on 2
png("/home/jmiller1/public_html/chr2.lod.qtl.png", width = 600)
plotRF(cross.NEW,chr=2)
dev.off()

# All chr not scaled
png("/home/jmiller1/public_html/ggplot2.qtl.png", width = 3000)
ggplot(melted, aes(y = as.factor(pop), x = pos, height = lod, fill = pop)) + ylab(NULL) +
  geom_ridgeline(stat = "identity", scale = 0.2) + facet_wrap(~as.factor(chr),
  scales = "free_x", nrow = 1, ncol = 24, strip.position = NULL) + scale_y_discrete("lod") +
  scale_fill_manual(name = "Populations", guide = F, values = popcol) + theme(axis.text = element_text(size = 10),
  axis.text.y = element_blank())
dev.off()

# All chr Not scaled either
png("/home/jmiller1/public_html/ggplot2.qtl.png", width = 3000)
ggplot(melted, aes(y = as.factor(pop), x = pos, height = lod, fill = pop)) + geom_ridgeline(stat = "identity",
  alpha = 0.5, scale = 0.25) + facet_wrap(~as.factor(chr), scales = "free_x", nrow = 1,
  ncol = 24, strip.position = NULL) + # scale_x_continuous('pos') + scale_y_continuous('lod') +
scale_y_discrete("lod") + theme(axis.text = element_text(size = 10))
dev.off()

# gene expression
png("/home/jmiller1/public_html/ggplot2.expr.qtl.png", width = 3000)
ggplot(melting, aes(x = pos1,y=-logp.DNT, height = logp.DNT)) +
geom_point( aes(color=as.factor(chem)), alpha=0.8, size=1.3) +
facet_wrap(~as.factor(chr), scales = "free_x", nrow = 1,
  ncol = 24, strip.position = NULL) + # scale_x_continuous('pos') + scale_y_continuous('lod') +

scale_x_continuous(label = chr) +
scale_y_continuous(expand = c(0, 0) ))    # remove space between plot area and x axis
dev.off()


# Phenotype figure
##### For Phenotype plots Phenotype plots Phenotype density
pheno.all <- phen <- read.csv("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/metadata/ALL_phenotype_Dist.csv",
  header = T)
phen <- phen[!is.na(phen$pheno_all), ]
phen$Pheno_05 <- phen$pheno_all
# index <- which(phen$pop_all == 'BRP')
phen$pop_all <- factor(phen$pop_all, levels = c("NBH", "BP", "NEW", "ELR"))

png("/home/jmiller1/public_html/phenotypes.png", width = 600)
ggplot(phen, aes(Pheno_05, y = pop_all, fill = pop_all, group = pop_all, height = ..density..,
  show.legend = F, order = as.numeric(pop_all))) +
  geom_density_ridges(scale = 1.33, bandwidth = 0.5) + xlab("Phenotype Score") + ylab("Density") + scale_alpha(guide = "none") +
  scale_y_discrete(limits = rev(levels(phen$pop_all)), expand = c(0.01, 0.01)) +
  scale_x_continuous(expand = c(0.02, 0.02), limits = c(0, 5)) + scale_fill_manual(name = "Populations",
  guide = "legend", values = popcol) + theme(text = element_text(size = 20), axis.text.y = element_blank(),
  panel.background = element_rect(fill = "white", colour = "black", size = 0.1,
    linetype = 1), )
dev.off()


#### Phenotype by genotype plots
outdir <-  '~/QTL_Map_Raw/popgen/rQTL/plots'
plot_pxg





### line for gene expression (each pcb/pah), pcb popgen outliers, AHR gene locations












#
# guess_phase - turn imputed genotypes into phased genotypes along chromosomes
# sim_geno - multiple imputations of underlying genotypes given marker data
# predict_snpgeno - predict SNP genotypes in a multiparent population from
# inferred genotypes plus founder strains’ SNP alleles.

perms.unstrat <- scan1perm(pr, cross2$pheno, model = "binary", cores = 12, n_perm = 200)

par(mar = c(4.1, 4.1, 0.6, 0.6))
plot_pxg(nbh2, iron$pheno[, "liver"], ylab = "Liver phenotype")


pxg1 <- find_marker(nbh2, chr = "2", 198)

png("/home/jmiller1/public_html/Q1.png", width = 300)
plotPXG(nbh2, pheno.col = 2)
