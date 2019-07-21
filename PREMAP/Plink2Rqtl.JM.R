#!/bin/R

library(qtl)

source("~/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/PREMAP/PLINK2RQTL.f2.R")
setwd("/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops")
dir <- "/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops/"

for (name in c("NBH", "ELR", "BRP", "NEW")) {
  ped <- paste(dir, name, ".ped", sep = "")
  map <- paste(dir, name, ".map", sep = "")
  PLINKtoCSVR(ped = ped, map = map, out = paste(dir, name, ".unphased.f2.csvr", 
    sep = ""))
}

dir <- "/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops/"
for (name in c("NBH", "NEW", "BRP1", "BRP8", "ELR")) {
  ped <- paste(dir, name, ".filt.pk.recode.ped", sep = "")
  map <- paste(dir, name, ".filt.pk.recode.map", sep = "")
  PLINKtoCSVR(ped = ped, map = map, out = paste(dir, name, ".parents.csvr", sep = ""))
}

name <- "NBH.um"
ped <- paste(dir, name, ".ped", sep = "")
map <- paste(dir, name, ".map", sep = "")
PLINKtoCSVR(ped = ped, map = map, out = paste(dir, name, ".unmapped.f2.csvr", sep = ""))
