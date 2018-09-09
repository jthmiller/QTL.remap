## rQTL2
require(qtl2,lib.loc='/share/apps/rmodules')

## Directories
basedir <- '/home/jmiller1/QTL_Map_Raw/popgen'
plotdir <- file.path(basedir,'rQTL/plots')
indpops <- file.path(basedir,'plinkfiles/ind.pops')
qtldir <- file.path(basedir,'rQTL/remap_out')
setwd(qtldir)

## Funtions for processing rQTL data
source(file.path(basedir,'rQTL/scripts/QTL_remap/removeDoubleXO.R'))
source(file.path(basedir,'rQTL/scripts/QTL_remap/QTL_map_sourcefile.R'))

## QTL LGs to consider
X <- c(1,2,8,13,18,24)
X <- c(1,2,8,13,18)
