### Switching order of markers manually
load('~/Dropbox/QTL_Paper/Rough Figures/LOD profiles/MAP/NEW.mapped.Rsave')
mapthis <- cross.18
12  13  14  15  16  17  18  19  20
134 144 181 142 223 182 179 249 274
21  22  23  24
250 179 253 204

i <- 12

plotRF(cross.18,chr=i)
abline(v=81,col='black')
abline(v=99,col='red')
abline(v=215,col='blue')
abline(v=210,col='green')

mapthis <- switch.order(mapthis, chr=i, c(1:79,99:134,80:98), error.prob=0.02,
  map.function="kosambi",maxit=1000, tol=1e-6, sex.sp=F)

plotRF(mapthis,chr=i)

try <- ripple(cross, chr, window=10, method="likelihood",
    error.prob=0.01, map.function="kosambi",maxit=1,
    tol=1e-6, sex.sp=TRUE, verbose=TRUE, n.cluster=1)
