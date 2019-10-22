


library('qtl')
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'




new <- file.path(mpath,'new.mapped.tsp.csv')
nbh <- file.path(mpath,'nbh.mapped.tsp.csv')
elr <- file.path(mpath,'elr.mapped.tsp.csv')
brp <- file.path(mpath,'brp.mapped.tsp.csv')

new <- read.cross(file = new, format = "csv", genotypes=c("1","2","3"),estimate.map = FALSE)
nbh <- read.cross(file = nbh, format = "csv", genotypes=c("1","2","3"),estimate.map = FALSE)
elr <- read.cross(file = elr, format = "csv", genotypes=c("1","2","3"),estimate.map = FALSE)
brp <- read.cross(file = brp, format = "csv", genotypes=c("1","2","3"),estimate.map = FALSE)

gts.new <- geno.table(new)
gts.nbh <- geno.table(nbh)
gts.elr <- geno.table(elr)
gts.brp <- geno.table(brp)
