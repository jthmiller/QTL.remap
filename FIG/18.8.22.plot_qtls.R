#!/bin/R
### first run combine pops for multi-pop cross objects
debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library("ggridges")
library("plyr")
library("scales")
library("ggrepel")
#### AHRs #####
AHR.bed <- read.table("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/data/lift_AHR_genes.bed",
  stringsAsFactors = F, header = F)
colnames(AHR.bed) <- c("chrom", "str", "stp", "gene")
AHR.bed$chrom <- as.numeric(gsub("chr", "", AHR.bed$chrom))
AHR.bed$str <- as.numeric(AHR.bed$str)
AHR.bed$stp <- as.numeric(AHR.bed$stp)
AHR.notmap <- AHR.bed[is.na(AHR.bed$chrom), ]
AHR.bed <- AHR.bed[!is.na(AHR.bed$chrom), ]
AHR.bed$gene <- gsub(":158640", "", AHR.bed$gene)
# add arnts (forgot to scan for them)
## Phenotypes
##############
popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NBH/REMAPS"
cross.NBH <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv",
  sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))

popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/ELR/REMAPS"
cross.ELR <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv",
  sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))

popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NEW/REMAPS"
cross.NEW <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv",
  sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))

popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/BRP/REMAPS"
cross.BRP <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv",
  sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))

#################
mak <- markernames(cross.NEW)
cross.NEW <- switchAlleles(cross.NEW, markers = mak)
#################
cross.nbh <- sim.geno(cross.NBH, n.draws = 500, step = 5, off.end = 10, error.prob = 0.01,
  map.function = "kosambi", stepwidth = "fixed")
cross.new <- sim.geno(cross.NEW, n.draws = 500, step = 5, off.end = 10, error.prob = 0.01,
  map.function = "kosambi", stepwidth = "fixed")
cross.elr <- sim.geno(cross.ELR, n.draws = 500, step = 5, off.end = 10, error.prob = 0.05,
  map.function = "kosambi", stepwidth = "fixed")
cross.brp <- sim.geno(brp.remap, n.draws = 500, step = 5, off.end = 10, error.prob = 0.1,
  map.function = "kosambi", stepwidth = "fixed")
brp.remap
#################
sex <- read.table(file = file.path(dirso, "sex.txt"))
rownames(sex) <- sex$ID
cross.nbh$pheno$sex <- sex[as.character(cross.nbh$pheno$ID), 2]
cross.nbh$pheno$binary <- as.numeric(cross.nbh$pheno$pheno >= 3)
cross.new$pheno$sex <- sex[as.character(cross.new$pheno$ID), 2]
cross.new$pheno$binary <- as.numeric(cross.new$pheno$pheno >= 3)
cross.elr$pheno$sex <- sex[as.character(cross.elr$pheno$ID), 2]
cross.elr$pheno$binary <- as.numeric(cross.elr$pheno$pheno >= 3)
cross.brp$pheno$sex <- sex[as.character(cross.brp$pheno$ID), 2]
cross.brp$pheno$binary <- as.numeric(cross.brp$pheno$pheno >= 3)

cross.brp <- subset(cross.brp, ind = cross.brp$pheno$gt == 1)

cross.nbh <- reduce2grid(cross.nbh)
cross.new <- reduce2grid(cross.new)
cross.elr <- reduce2grid(cross.elr)
cross.brp <- reduce2grid(cross.brp)

# scan.norm.imp.NBH <- scanone(cross.nbh, model = 'normal', pheno.col = 1, method
# = 'imp', addcovar = cross.nbh$pheno$sex)
scan.norm.imp.NBH <- scanone(cross.nbh, model = "normal", pheno.col = 1, method = "imp")
# scan.norm.imp.NEW <- scanone(cross.new, model = 'normal', pheno.col = 1, method
# = 'imp', addcovar = cross.new$pheno$sex)
scan.norm.imp.NEW <- scanone(cross.new, model = "normal", pheno.col = 1, method = "imp")
# scan.norm.imp.ELR <- scanone(cross.elr, model = 'normal', pheno.col = 1, method
# = 'imp', addcovar = cross.elr$pheno$sex)
scan.norm.imp.ELR <- scanone(cross.elr, model = "normal", pheno.col = 1, method = "imp")
scan.norm.imp.BRP <- scanone(cross.brp, model = "normal", pheno.col = 1, method = "imp")
### use scanone for plots
themelt.nbh <- scan.norm.imp.NBH
themelt.new <- scan.norm.imp.NEW
themelt.elr <- scan.norm.imp.ELR
themelt.brp <- scan.norm.imp.BRP

themelt.nbh$pop <- "NBH"
themelt.new$pop <- "NEW"
themelt.elr$pop <- "ELR"
themelt.brp$pop <- "BRP"

### ggplot format AHR genes
nbh.gens <- cnv.ahrs(cross2 = cross.nbh, AHRdf = AHR.bed, EXP = F)
new.gens <- cnv.ahrs(cross.new, AHRdf = AHR.bed, EXP = F)
elr.gens <- cnv.ahrs(cross.elr, AHRdf = AHR.bed, EXP = F)
brp.gens <- cnv.ahrs(cross.brp, AHRdf = AHR.bed, EXP = F)
qtl.gens <- nbh.gens[which(nbh.gens$chr %in% c(1, 2, 5, 8, 10, 12, 13, 18, 24)),]
minor.gens <- nbh.gens[which(nbh.gens$chr %in% c(8, 13, 23, 24)), ]
incompat.gens <- nbh.gens[which(nbh.gens$chr %in% c(8, 13)), ]
qtl_pg <- c(2,8, 13, 18, 24)
ol.gens <- nbh.gens[which(nbh.gens$chr %in% qtl_pg), ]


### ggplot popgen locations
dir <- "/home/jmiller1/QTL_Map_Raw/popgen/tables"
nbh.popgen <- read.table(file.path(dir, "outliersNBH.txt.ncbi.lifted"), sep = "\t",
  header = T)
new.popgen <- read.table(file.path(dir, "outliersNYC.txt.ncbi.lifted"), sep = "\t",
  header = T)
elr.popgen <- read.table(file.path(dir, "outliersER.txt.ncbi.lifted"), sep = "\t",
  header = T)
brp.popgen <- read.table(file.path(dir, "outliersBP.txt.ncbi.lifted"), sep = "\t",
  header = T)

### Use nbh coords but elr and new popgen
new.rank <- cnv.popgen(cross.nbh, new.popgen, top = 50)
nbh.rank <- cnv.popgen(cross.nbh, nbh.popgen, top = 50)
elr.rank <- cnv.popgen(cross.nbh, elr.popgen, top = 50)
brp.rank <- cnv.popgen(cross.nbh, brp.popgen, top = 50)
nbh.rank$pop <- "NBH"
new.rank$pop <- "NEW"
elr.rank$pop <- "ELR"
brp.rank$pop <- "BRP"

all.rank <- rbind(new.rank, nbh.rank, elr.rank, brp.rank)
all.rank$pop <- factor(all.rank$pop, levels = c("NBH", "BRP", "NEW", "ELR"))
qtl.rank <- all.rank[which(all.rank$chr %in% c(1, 2, 5, 8, 10, 12, 13, 18, 23, 24)),
  ]
minor.rank <- all.rank[which(all.rank$chr %in% c(8, 13, 23, 24)), ]
incompat.rank <- all.rank[which(all.rank$chr %in% c(8, 13)), ]

qtl_pg <- c(2,8, 13, 18, 24)
ol.rank <- all.rank[which(all.rank$chr %in% qtl_pg), ]

### GGriges plot
melted.nbh <- data.frame(pop = "NBH", chr = scan.norm.imp.NBH$chr, pos = scan.norm.imp.NBH$pos,
  lod = scan.norm.imp.NBH$lod)
melted.new <- data.frame(pop = "NEW", chr = scan.norm.imp.NEW$chr, pos = scan.norm.imp.NEW$pos,
  lod = scan.norm.imp.NEW$lod)
melted.elr <- data.frame(pop = "ELR", chr = scan.norm.imp.ELR$chr, pos = scan.norm.imp.ELR$pos,
  lod = scan.norm.imp.ELR$lod)
melted.brp <- data.frame(pop = "BRP", chr = scan.norm.imp.BRP$chr, pos = scan.norm.imp.BRP$pos,
  lod = scan.norm.imp.BRP$lod)

melted <- rbind(melted.nbh, melted.new, melted.elr, melted.brp)
melted$pop <- factor(melted$pop, levels = rev(c("NBH", "BRP", "NEW", "ELR")))

## Total CM length of NBH. Rescale to NBH
mxes <- sapply(1:24, function(X) {
  max(themelt.nbh$pos[which(themelt.nbh$chr == X)])
})

ts <- themelt.new[which(themelt.new$chr == 1), ]
ts$pos <- rescale(ts$pos, to = c(-10, mxes[1]))
new.rescale <- ts
for (i in 2:24) {
  ts <- themelt.new[which(themelt.new$chr == i), ]
  ts$pos <- rescale(ts$pos, to = c(-10, mxes[i]))
  new.rescale <- rbind(new.rescale, ts)
}

ts <- themelt.elr[which(themelt.elr$chr == 1), ]
ts$pos <- rescale(ts$pos, to = c(-10, mxes[1]))
elr.rescale <- ts
for (i in 2:24) {
  ts <- themelt.elr[which(themelt.elr$chr == i), ]
  ts$pos <- rescale(ts$pos, to = c(-10, mxes[i]))
  elr.rescale <- rbind(elr.rescale, ts)
}

ts <- themelt.brp[which(themelt.brp$chr == 1), ]
ts$pos <- rescale(ts$pos, to = c(-10, mxes[1]))
brp.rescale <- ts
for (i in 2:24) {
  ts <- themelt.brp[which(themelt.brp$chr == i), ]
  ts$pos <- rescale(ts$pos, to = c(-10, mxes[i]))
  brp.rescale <- rbind(brp.rescale, ts)
}

allmelt <- rbind(themelt.nbh, new.rescale, elr.rescale, brp.rescale)
allmelt$pop <- factor(allmelt$pop, levels = c("NBH", "BRP", "NEW", "ELR"))
qtlmelt <- allmelt[which(allmelt$chr %in% c(1, 2, 5, 8, 10, 12, 13, 18,19, 23, 24)),
  ]
qtlminor <- allmelt[which(allmelt$chr %in% c(8, 13,19, 23, 24)), ]
incompat <- allmelt[which(allmelt$chr %in% c(8, 13)), ]

qtl_pg <- c(2,8, 13, 18, 24)
ol.melt <- allmelt[which(allmelt$chr %in% qtl_pg), ]

#### Subest for only qtl plots qtl.rank <- nbh.rank[which(nbh.rank$chr %in%
#### c(1,2,5,8,10,12,13,18,24)),]


#### MAP BRP to NBH before simcross
brp.remap <- cnv.premap(cross.nbh, cross.BRP)
##save.image("/home/jmiller1/public_html/BRP_remap.Rsave")


save.image('/home/jmiller1/public_html/QTL_plot.Rsave')


p <- ggplot(themelt.nbh, aes(x = pos, y = lod))
png("/home/jmiller1/public_html/all_popgen_rank_scaled.qtl.png", width = 3000)
p + facet_wrap(~chr, scales = "free_x", nrow = 1, ncol = 24) + scale_y_continuous(limits = c(-12,
  23)) + geom_line(size = 2, alpha = 0.6) + theme_minimal() + theme(axis.text = element_text(size = 10)) +
  labs(x = "Chromosome", y = "LOD", linetype = "") + geom_label_repel(aes(x = pos,
  y = -0.1, label = gene), box.padding = unit(0.25, "lines"), parse = T, point.padding = unit(0.2,
  "lines"), force = 10, label.padding = unit(0.2, "lines"), ylim = c(0, -12), segment.size = 1,
  max.iter = 6000, data = nbh.gens, direction = "y", size = 4, seed = 666, nudge_y = -0.01,
  vjust = 3) + geom_label_repel(aes(x = pos, y = -0.1, label = gene), box.padding = unit(0.25,
  "lines"), parse = T, point.padding = unit(0.2, "lines"), force = 10, label.padding = unit(0.2,
  "lines"), ylim = c(0, -12), segment.size = 0, max.iter = 6000, data = nbh.gens,
  direction = "y", size = 4, seed = 666, nudge_y = -0.01, vjust = 3) + geom_point(aes(size = rank),
  data = nbh.rank) + geom_label_repel(aes(x = pos, y = lod, label = rank), data = nbh.rank,
  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
  vjust = 1)
dev.off()

p <- ggplot(themelt.new, aes(x = pos, y = lod))
png("/home/jmiller1/public_html/new_genes_below_rank_above.qtl.png", width = 3000)
p + facet_wrap(~chr, scales = "free_x", nrow = 1, ncol = 24) + scale_y_continuous(limits = c(-12,
  23)) + geom_line(size = 2, alpha = 0.6) + theme_minimal() + theme(axis.text = element_text(size = 10)) +
  labs(x = "Chromosome", y = "LOD", linetype = "") + geom_label_repel(aes(x = pos,
  y = -0.1, label = gene), box.padding = unit(0.25, "lines"), parse = T, point.padding = unit(0.2,
  "lines"), force = 10, label.padding = unit(0.2, "lines"), ylim = c(0, -12), segment.size = 1,
  max.iter = 6000, data = new.gens, direction = "y", size = 4, seed = 666, nudge_y = -0.01,
  vjust = 3) + geom_label_repel(aes(x = pos, y = -0.1, label = gene), box.padding = unit(0.25,
  "lines"), parse = T, point.padding = unit(0.2, "lines"), force = 10, label.padding = unit(0.2,
  "lines"), ylim = c(0, -12), segment.size = 0, max.iter = 6000, data = new.gens,
  direction = "y", size = 4, seed = 666, nudge_y = -0.01, vjust = 3) + geom_label_repel(aes(x = pos,
  y = lod, label = rank), data = nbh.rank, size = 4, box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"), vjust = 1)
dev.off()

p <- ggplot(themelt.elr, aes(x = pos, y = lod, color = pop))
png("/home/jmiller1/public_html/rank_top_scaled_random.qtl.png", width = 3000)
p + facet_wrap(~chr, scales = "free_x", nrow = 1, ncol = 24) + scale_y_continuous(limits = c(-12,
  23)) + scale_color_manual(values = popcol) + geom_line(size = 2, alpha = 0.6) +
  theme_minimal() + theme(axis.text = element_text(size = 10)) + labs(x = "Chromosome",
  y = "LOD", linetype = "") + geom_label_repel(aes(x = pos, y = -0.1, label = gene),
  box.padding = unit(0.25, "lines"), parse = T, point.padding = unit(0.2, "lines"),
  force = 10, label.padding = unit(0.2, "lines"), ylim = c(0, -12), segment.size = 1,
  max.iter = 6000, data = elr.gens, direction = "y", size = 4, seed = 666, nudge_y = -0.01,
  vjust = 3) + geom_label_repel(aes(x = pos, y = -0.1, label = gene), box.padding = unit(0.25,
  "lines"), parse = T, point.padding = unit(0.2, "lines"), force = 10, label.padding = unit(0.2,
  "lines"), ylim = c(0, -12), segment.size = 0, max.iter = 6000, data = elr.gens,
  direction = "y", size = 4, seed = 666, nudge_y = -0.01, vjust = 3) + geom_label_repel(aes(x = pos,
  y = lod, label = rank), data = elr.rank, size = 4, box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"), vjust = 1)
dev.off()

p <- ggplot(allmelt, aes(x = pos, y = lod, colour = pop))
png("/home/jmiller1/public_html/no_annot_.qtl.png", width = 2000)
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 24) + scale_y_continuous(limits = c(0,
  22)) + scale_color_manual(values = popcol) + geom_line(size = 2, alpha = 0.85) +
  theme_minimal() + theme(axis.text = element_text(size = 10)) + labs(x = "Chromosome",
  y = "LOD", linetype = "")
dev.off()

p <- ggplot(allmelt, aes(x = pos, y = lod, colour = pop))
png("/home/jmiller1/public_html/ahr_genes_only.qtl.png", width = 2000)
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 24) + scale_y_continuous(limits = c(0,
  22)) + scale_color_manual(values = popcol) + geom_line(size = 2, alpha = 0.6) +
  theme_minimal() + theme(axis.text = element_text(size = 10)) + geom_label_repel(aes(x = pos,
  y = 0.1, label = gene), box.padding = unit(0.5, "lines"), parse = T, point.padding = unit(0.8,
  "lines"), force = 1, label.padding = unit(0.2, "lines"), ylim = c(5, 20), segment.size = 1,
  max.iter = 6000, data = nbh.gens, direction = "y", size = 3, seed = 666, nudge_y = 0.01,
  vjust = 0.1) + geom_label_repel(aes(x = pos, y = 0.1, label = gene), box.padding = unit(0.5,
  "lines"), parse = T, point.padding = unit(0.8, "lines"), force = 1, label.padding = unit(0.2,
  "lines"), ylim = c(5, 20), segment.size = 0, max.iter = 6000, data = nbh.gens,
  direction = "y", size = 3, seed = 666, nudge_y = 0.01, vjust = 0.1) + labs(x = "Chromosome",
  y = "LOD", linetype = "")
dev.off()

p <- ggplot(allmelt, aes(x = pos, y = lod, colour = pop))
png("/home/jmiller1/public_html/All_chr.qtl.png", width = 3000)
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 24) + scale_color_manual(values = popcol) +
  scale_y_continuous(limits = c(-5, 22)) + geom_label_repel(aes(x = pos, y = 0.1,
  label = gene), color = "black", segment.size = 1, label.padding = unit(0.2, "lines"),
  data = nbh.gens, direction = "y", size = 3, seed = 666, ylim = c(5, 20), segment.alpha = 0.5) +
  geom_line(size = 2, alpha = 0.5) + geom_label_repel(aes(x = pos, y = 0.1, label = gene),
  color = "black", label.padding = unit(0.2, "lines"), ylim = c(5, 20), segment.size = 0,
  max.iter = 6000, nudge_y = 2, data = nbh.gens, direction = "y", size = 3.5, seed = 666) +
  geom_label_repel(aes(x = pos, y = 0, fill = pop, label = rank), data = all.rank,
    size = 3.5, segment.size = 1, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
      "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  geom_label_repel(aes(x = pos, y = 0, label = rank, ), fontface = "bold", data = all.rank,
    size = 3.5, segment.size = 0, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
      "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = "none") + labs(x = "Chromosome", y = "LOD", linetype = "")
dev.off()


qtlmelt <- allmelt[which(allmelt$chr %in% c( 2, 13, 18)),]
qtl.rank <- qtl.rank[which(qtl.rank$rank <21),]

qtl.rank <- all.rank[which(all.rank$chr %in% c(2, 13, 18)),]
qtl.gens <- nbh.gens[which(nbh.gens$chr %in% c( 2, 13, 18)),]
qtl.gens <- qtl.gens[as.character(c(1,3,4,5,6,7,8,24,26,45,46,56,70,96,94,57)),]
qtl.gens <- qtl.gens[c(1,3,4,5,6,7),]
qtl.gens <- qtl.gens[-5,]

p <- ggplot(qtlmelt, aes(x = pos, y = lod, colour = pop))
pdf("/home/jmiller1/public_html/all_pop_qtl_only.pdf",width=15,height=11)
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 10) +
  scale_color_manual(values = popcol) +
  scale_y_continuous(limits = c(-5, 22)) +
  ### gene line
  geom_label_repel(aes(x = pos, y = 0.1, label = gene),fill='white', color = "black", segment.size = 1,
    label.padding = unit(0.2, "lines"), data = qtl.gens, direction = "y", size = 5,
    seed = 666, nudge_y = 2, ylim = c(5, 20), segment.alpha = 0.5) + ## plot lines

  geom_line(size = 2, alpha = 0.75) + ### gene label

  geom_label_repel(aes(x = pos, y = 0.1, label = gene), color = "black", segment.size = 0,
    label.padding = unit(0.2, "lines"), data = qtl.gens, direction = "y", size = 5,
    seed = 666, nudge_y = 2, ylim = c(5, 20), max.iter = 6000, ) + ### Rank line

  geom_label_repel(aes(x = pos, y = 0, color = pop, label = rank), data = qtl.rank,
    size = 5, segment.size = 1, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
    "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +
    ### Rank label

  geom_label_repel(aes(x = pos, y = 0, color = pop,label = rank), fontface = "bold", data = qtl.rank,
    size = 5, segment.size = 0, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
    "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +

  theme(axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", color = "black", size = 16),
    axis.text.y=element_text(face = "bold", color = "black", size = 16),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(face = "bold", color = "black",size=16),
    legend.position = "none") +
  labs(x = "Chromosome", y = "LOD", linetype = "")
dev.off()
## Why would chromosomes with top ranked loci have fewest signatures?

p <- ggplot(qtlminor, aes(x = pos, y = lod, colour = pop))
png("/home/jmiller1/public_html/minor_qtl_only.png", width = 2000)
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 9) + scale_color_manual(values = popcol) +
  scale_y_continuous(limits = c(-1, 8)) + theme_minimal() + geom_label_repel(aes(x = pos,
  y = 0.1, label = gene), color = "black", segment.size = 1, label.padding = unit(0.2,
  "lines"), data = minor.gens, direction = "y", size = 5, seed = 666, nudge_y = 2,
  ylim = c(5, 20), segment.alpha = 0.5) + geom_line(size = 2, alpha = 0.5) + geom_label_repel(aes(x = pos,
  y = 0.1, label = gene), color = "black", label.padding = unit(0.2, "lines"),
  ylim = c(5, 20), segment.size = 0, max.iter = 6000, nudge_y = 2, data = minor.gens,
  direction = "y", size = 5, seed = 666) + geom_label_repel(aes(x = pos, y = 0,
  color = pop, label = rank), data = minor.rank, size = 5, segment.size = 1, force = 4,
  min.segment.length = 0.1, point.padding = unit(0.4, "lines"), direction = "both",
  ylim = c(0, -8), seed = 666, box.padding = 0.1) + geom_label_repel(aes(x = pos,
  y = 0, label = rank), fontface = "bold", data = minor.rank, size = 5, segment.size = 0,
  force = 4, min.segment.length = 0.1, point.padding = unit(0.4, "lines"), direction = "both",
  ylim = c(0, -8), seed = 666, box.padding = 0.1) + theme(axis.title.x = element_blank(),
  axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  labs(x = "Chromosome", y = "LOD", linetype = "")
dev.off()

p <- ggplot(incompat, aes(x = pos, y = lod, colour = pop))
png("/home/jmiller1/public_html/elr_incompat.png", width = 750)
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 9) +
  theme(strip.text.x = element_text(size = 12)) +
  scale_color_manual(values = popcol) +
  scale_y_continuous(limits = c(-1, 8)) +
  geom_label_repel(aes(x = pos, y = 0.1, label = gene), color = "black", segment.size = 1,
    label.padding = unit(0.2, "lines"), data = incompat.gens, direction = "y",
    size = 5, seed = 666, nudge_y = 2, ylim = c(5, 20), segment.alpha = 0.5) +
  geom_line(size = 2, alpha = 0.5) + geom_label_repel(aes(x = pos, y = 0.1, label = gene),
  color = "black", label.padding = unit(0.2, "lines"), ylim = c(5, 20), segment.size = 0,
  max.iter = 6000, nudge_y = 2, data = incompat.gens, direction = "y", size = 5,
  seed = 666) +
  geom_label_repel(aes(x = pos, y = 0, color = pop, label = rank),
  data = incompat.rank, size = 5, segment.size = 1, force = 4, min.segment.length = 0.1,
  point.padding = unit(0.4, "lines"), direction = "both", ylim = c(0, -8), seed = 666,
  box.padding = 0.1) +
  geom_label_repel(aes(x = pos, y = 0, label = rank), fontface = "bold",
  data = incompat.rank, size = 5, segment.size = 0, force = 4, min.segment.length = 0.1,
  point.padding = unit(0.4, "lines"), direction = "both", ylim = c(0, -8), seed = 666,
  box.padding = 0.1) + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
  axis.ticks.x = element_blank(), legend.position = "none") + labs(x = "Chromosome",
  y = "LOD", linetype = "")
dev.off()

### QTL only

p <- ggplot(ol.melt, aes(x = pos, y = lod, colour = pop))
png("/home/jmiller1/public_html/ol_qtl_only.png", width = 2000)
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 9) +
scale_color_manual(values = popcol) +
  scale_y_continuous(limits = c(-1, 8)) +
  theme_minimal() + geom_label_repel(aes(x = pos,
  y = 0.1, label = gene), color = "black", segment.size = 1, label.padding = unit(0.2,
  "lines"), data = minor.gens, direction = "y", size = 5, seed = 666, nudge_y = 2,
  ylim = c(5, 20), segment.alpha = 0.5) +
  geom_line(size = 2, alpha = 0.5) +
  geom_label_repel(aes(x = pos,y = 0.1, label = gene),
  color = "black", label.padding = unit(0.2, "lines"),
  ylim = c(5, 20), segment.size = 0, max.iter = 6000, nudge_y = 2, data = minor.gens,
  direction = "y", size = 5, seed = 666) +
  geom_label_repel(aes(x = pos, y = 0,color = pop, label = rank),
  data = minor.rank, size = 5, segment.size = 1, force = 4,
  min.segment.length = 0.1, point.padding = unit(0.4, "lines"), direction = "both",
  ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  geom_label_repel(aes(x = pos,y = 0, label = rank),
  fontface = "bold", data = minor.rank, size = 5, segment.size = 0,
  force = 4, min.segment.length = 0.1, point.padding = unit(0.4, "lines"), direction = "both",
  ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  theme(axis.title.x = element_blank(),
  axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  labs(x = "Chromosome", y = "LOD", linetype = "")
dev.off()


names(popcol)[2] <- 'BRP'

p <- ggplot(ol.melt, aes(x = pos, y = lod, colour = pop))
pdf("/home/jmiller1/public_html/OL_chr.qtl.pdf")
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 24) +
  scale_color_manual(values = popcol) +
  scale_y_continuous(limits = c(-5, 22)) +
  geom_label_repel(aes(x = pos, y = 0.1,label = gene),data = ol.gens,
    color = "black", segment.size = 1, label.padding = unit(0.2, "lines"),
    direction = "y", size = 3, seed = 666, ylim = c(5, 20), segment.alpha = 0.5) +
  geom_line(size = 2, alpha = 0.5) +
  geom_label_repel(aes(x = pos, y = 0.1, label = gene), data = ol.gens,
    color = "black", label.padding = unit(0.2, "lines"), ylim = c(5, 20), segment.size = 0,
    max.iter = 6000, nudge_y = 2, direction = "y", size = 3.5, seed = 666) +
  geom_label_repel(aes(x = pos, y = 0, fill = pop, label = rank), data = ol.rank,
    size = 3.5, segment.size = 1, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
      "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  geom_label_repel(aes(x = pos, y = 0, label = rank, ), fontface = "bold", data = ol.rank,
    size = 3.5, segment.size = 0, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
      "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  theme(axis.title.y = element_text(face='bold',size = 12),
    strip.text.x = element_text(face='bold',size = 12),
    axis.text.x = element_text(face='bold', size = 12),
    axis.text.y = element_text(face='bold', size = 12),
    legend.position = "none") +
  labs(x = "Chromosome", y = "LOD", linetype = "")
dev.off()


qtl_pg <- c(1,5, 10, 12, 23)
ol.rank <- all.rank[which(all.rank$chr %in% qtl_pg), ]
ol.melt <- allmelt[which(allmelt$chr %in% qtl_pg), ]
ol.gens <- nbh.gens[which(nbh.gens$chr %in% qtl_pg), ]


p <- ggplot(ol.melt, aes(x = pos, y = lod, colour = pop))
pdf("/home/jmiller1/public_html/OL_noQTL_chr.qtl.pdf")
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 24) +
  scale_color_manual(values = popcol) +
  scale_y_continuous(limits = c(-5, 22)) +
  geom_label_repel(aes(x = pos, y = 0.1,label = gene),data = ol.gens,
    color = "black", segment.size = 1, label.padding = unit(0.2, "lines"),
    direction = "y", size = 3, seed = 666, ylim = c(5, 20), segment.alpha = 0.5) +
  geom_line(size = 2, alpha = 0.5) +
  geom_label_repel(aes(x = pos, y = 0.1, label = gene), data = ol.gens,
    color = "black", label.padding = unit(0.2, "lines"), ylim = c(5, 20), segment.size = 0,
    max.iter = 6000, nudge_y = 2, direction = "y", size = 3.5, seed = 666) +
  geom_label_repel(aes(x = pos, y = 0, fill = pop, label = rank), data = ol.rank,
    size = 3.5, segment.size = 1, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
      "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  geom_label_repel(aes(x = pos, y = 0, label = rank, ), fontface = "bold", data = ol.rank,
    size = 3.5, segment.size = 0, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
      "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  theme(axis.title.y = element_text(face='bold',size = 12),
    strip.text.x = element_text(face='bold',size = 12),
    axis.text.x = element_text(face='bold', size = 12),
    axis.text.y = element_text(face='bold', size = 12),
    legend.position = "none") +
  labs(x = "Chromosome", y = "LOD", linetype = "")
dev.off()


### Entropy

NBH <- subset(cross.NBH, ind = cross.NBH$pheno$gt == 1)
png("/home/jmiller1/public_html/nbh_entropy.png", width = 3000)
plotInfo(NBH, chr = c(1:24), main = "NBH", method = "both", include.genofreq = T)
dev.off()

NEW <- subset(cross.NEW, ind = cross.NEW$pheno$gt == 1)
png("/home/jmiller1/public_html/new_entropy.png", width = 3000)
plotInfo(NEW, chr = c(1:24), main = "NEW", method = "both", include.genofreq = T)
dev.off()

ELR <- subset(cross.ELR, ind = cross.ELR$pheno$gt == 1)
png("/home/jmiller1/public_html/elr_entropy.png", width = 3000)
plotInfo(ELR, chr = c(1:24), main = "ELR", method = "both", include.genofreq = T)
dev.off()

BRP <- subset(brp.remap, ind = brp.remap$pheno$gt == 1)
png("/home/jmiller1/public_html/brp_entropy.png", width = 3000)
plotInfo(BRP, chr = c(1:24), main = "BRP", method = "both", include.genofreq = T)
dev.off()

#### Density

pheno <- read.table("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/data/pheno.csv",
  stringsAsFactors = T, header = T,sep=',')
pheno$pop_all <- factor(pheno$pop_all, levels = rev(c('NBH','BRP','NEW','ELR')))
pheno$pheno_all <- factor(pheno$pheno_all, levels = c(NA,0:5))
pheno$gtd <- pheno$pheno_all
pheno[which(pheno$GT_NG_ALT=='NG'),6] <- NA
