### Switching order of markers manually
source('~/Dropbox/QTL_Paper/CODE/QTL/local_debug.R')
load('~/Dropbox/QTL_Paper/Rough Figures/LOD profiles/MAP/NEW.mapped.Rsave')
1X   2   3   4   5   6   7   8   9  10  11  12
231X 261 194 205 178 254 176 265 206 194 166 134
#mapthis <- cross.18


i <- 2

plotRF(mapthis,chr=i)

abline(v=145,col='black')
abline(v=90,col='red')
abline(v=37,col='blue')
abline(v=95,col='black')

mapthis <- drop.markers(mapthis, markernames(mapthis,chr=i)[c(1:2)])

mapthis <- switch.order(mapthis, chr=i, c(1:150,157:158,151:156,159:261), error.prob=0.02,
  map.function="kosambi",maxit=1000, tol=1e-6, sex.sp=F)
plotRF(mapthis,chr=i,what='both',mark.diagonal=T,col.scheme="redblue")
save.image('~/Dropbox/QTL_Paper/Rough Figures/LOD/profiles/MAP/NEW.1_24.Rsave')





i <- 8
ord <- order(as.numeric(gsub(paste(i,":",sep=''),'',markernames(mapthis,chr=i))))

pos <- switch.order(mapthis, chr=i,ord , error.prob=0.02,
  map.function="kosambi",maxit=1000, tol=1e-6, sex.sp=F)


plotRF(pos,chr=i)
