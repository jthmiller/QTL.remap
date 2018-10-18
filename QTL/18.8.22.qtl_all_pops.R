#!/bin/bash
#source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')


NBH <- new.env()
load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NBH/REMAPS/QTLmap.Rsave', envir=NBH)
ELR <- new.env()
#ELR <- load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/ELR/REMAPS/QTLmap.Rsave', envir=ELR)
NEW <- new.env()
load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NEW/REMAPS/QTLmap.Rsave', envir=NEW)
### These markers worked in NBH,NEW. Maybe ELR?
marks <- unique(c(markernames(NBH$cross.18),markernames(NBH$cross.18)))

cross.18 <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
                format='csvr', geno=c(1:3),estimate.map=FALSE)

cross.18 <- drop.markers(cross.18,markernames(cross.18)[!markernames(cross.18) %in% marks])

gt <- geno.table(cross.18)
cutoff <- 1.0e-4
cross.18 <- drop.markers(cross.18,rownames(gt[gt$P.value<cutoff,]))
