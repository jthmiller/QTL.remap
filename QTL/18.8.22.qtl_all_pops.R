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
