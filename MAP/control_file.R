#!/bin/R

## Each chrom/pop info

if (exists("debug.cross")) {
  ####### DEBUG ONLY ####
  slurmcore <- 12
  setwd("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/")
  ### Prompt
  pop <- function() c("NBH", "NEW", "ELR", "BRP")[menu(c("NBH", "NEW", "ELR", "BRP"), 
    title = "Which pop")]
  pop <- popq <- pop()
  X <- function() c(1:24)[menu(1:24, title = "Which Chromosome")]
  X <- X()
  ## Only use previously mapped markers?
  mapped.only <- function() c(TRUE, FALSE)[menu(c(TRUE, FALSE), title = "Mapped markers only?")]
  mapped.only <- mapped.only()
  ## Only use granparent confirmed markers?
  confirmed <- function() c(TRUE, FALSE)[menu(c(TRUE, FALSE), title = "Only use granparent confirmed markers?")]
  confirmed <- confirmed()
  
  ### Libraries for plots
  flib <- "/share/apps/rmodules"
  fpacks <- c("ggplot2", "reshape", "pheatmap", "lattice")
  lapply(fpacks, require, character.only = TRUE, lib.loc = flib)
  mylib <- "/home/jmiller1/R/x86_64-pc-linux-gnu-library/3.5"
  mpacks <- c("gplots", "qgraph")
  lapply(mpacks, require, character.only = TRUE, lib.loc = mylib)
} else {
  pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c("NBH", "ELR", "NEW", "BP")]
  if (pop == "BP") {
    pop <- "BRP"
  }
  X <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))  ## X is equal to chrom number
  slurmcore <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
  popq <- commandArgs(TRUE)[commandArgs(TRUE) %in% c("NBH", "ELR", "NEW", "BP")]
  plotting <- F
}
## Directories
basedir <- "/home/jmiller1/QTL_Map_Raw/popgen"
plotdir <- file.path(basedir, "rQTL/plots")
indpops <- file.path(basedir, "plinkfiles/ind.pops")
popdir <- file.path(basedir, "rQTL", pop, "REMAPS")
qtldir <- file.path(basedir, "rQTL/remap_out")
errfile <- file.path(qtldir, "genotyping_error_rate.txt")
dirso <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/data/"

## Funtions for processing rQTL map data
source(file.path(basedir, "rQTL/scripts/QTL_remap/MAP/source_file.R"))
# source(file.path(basedir, 'rQTL/scripts/QTL_remap/QTL/model_source_file.R'))

## Libraries
flib <- "/share/apps/rmodules"
fpacks <- c("devtools", "httr", "RColorBrewer", "qtl")
lapply(fpacks, require, character.only = TRUE, lib.loc = flib)

mylib <- "/home/jmiller1/R/x86_64-pc-linux-gnu-library/3.5"
mpacks <- c("qtl", "foreach", "qtl2", "qtlTools", "doParallel", "plyr")
lapply(mpacks, require, character.only = TRUE, lib.loc = mylib)

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if (trace) 
      cat(nm, ":")
    source(file.path(path, nm), ...)
    if (trace) 
      cat("\n")
  }
}

sourceDir("doParallel/R")
## Pop vars
chrms <- c(1:24)
pops <- c("NBH", "NEW", "ELR")
all.pops <- c("NBH", "BRP", "ELR", "NEW")
popcol <- brewer.pal(8, "Paired")[c(2, 4, 6, 8)]
names(popcol) <- all.pops

### Phenotype translation
trsl.bin <- c(0, 0, 0, 1, 1, 1)
names(trsl.bin) <- as.character(0:5)

hoods <- F
reprip <- F
## Parameters for rQTL for population specific datasets (NBH markers require at
## least 70% genotypes )
if (pop == "NBH") {
  ns <- "North"
  confirmed = T
  reorder.marks <- F
  mapped.only = TRUE
  grpLod <- 10  ## Standard LG form LOD
  finLod <- 12  ## Higher final NBH LOD
  grpRf <- 0.3
  finRf <- 0.15
  cutoff <- 1e-05
  miss <- 12
  miss1 <- 10
  miss2 <- 8
  droppo <- 3
} else if (pop == "ELR") {
  ns <- "South"
  confirmed = T
  reorder.marks <- F
  mapped.only = TRUE
  missing <- 0.9
  grpLod <- 10  ## Standard LG form LOD
  finLod <- 12  ## Higher final ELR LOD
  grpRf <- 0.3
  finRf <- 0.15
  cutoff <- 1e-04  ## Higher, need more power to detect seg distortion
  miss <- 10
  miss1 <- 10
  miss2 <- 10
  droppo <- 15
} else if (pop == "NEW") {
  ns <- "North"
  confirmed = T
  reorder.marks <- F
  mapped.only = TRUE
  inds <- c(NA)  # determined to be dropped low cov
  missing <- 0.8
  grpLod <- 10  ## Standard LG form LOD
  finLod <- 12  ## Higher final ELR LOD
  grpRf <- 0.3
  finRf <- 0.15
  cutoff <- 1e-05
  miss <- 15
  miss1 <- 10
  miss2 <- 8
  droppo <- 2
} else if (pop == "BRP") {
  ns <- "North"
  confirmed = T
  reorder.marks <- F
  mapped.only = TRUE
  inds <- c(NA)  # determined to be dropped low cov
  missing <- 0.8
  grpLod <- 10  ## Standard LG form LOD
  finLod <- 12  ## Higher final ELR LOD
  grpRf <- 0.3
  finRf <- 0.15
  cutoff <- 1e-05
  miss <- 15
  miss1 <- 10
  miss2 <- 8
  droppo <- 2
}

if (mapped.only == TRUE & reorder.marks == F) {
  outname <- "NW_dropped_physical"
} else if (mapped.only == TRUE & reorder.marks == T) {
  outname <- "NW_Mapped"
} else {
  outname <- "NW_ReMapped"
}

if (X %in% c(8, 5)) {
  reorder.marks <- T
}

## Try to get error exported by map
expr <- paste("tac ", errfile, " | grep -m 1 '", pop, " ", X, "' | awk '{print $5}'", 
  sep = "")
try(ers <- as.numeric(system(expr, intern = T)))

if (length(ers) == 0 | is.null(ers)) {
  print("couldnt find the error. Using 0.03")
  ers <- 0.03
} else {
  print(paste(ers, "genotyping error"))
}
