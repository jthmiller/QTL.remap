try(
cross.no2 <- subset(cross,chr=c(1,3:24))
registerDoParallel(slurmcore)
operm.hk <- foreach(i=50,
  .combine=c,.packages = "qtl") %dopar% {
    operm <- scantwo(cross.no2, method="hk", n.perm=1,perm.strata=strata)
}
pen.em <- calc.penalties(operm.hk)
save.image(paste('QTLmap.Rsave',sep=''))
)
try(
### qtl scan2
registerDoParallel(slurmcore)
operm.em <- foreach(i=50,
  .combine=c,.packages = "qtl") %dopar% {
    operm <- scantwo(cross.no2, method="em", n.perm=1,perm.strata=strata)
}
pen.em <- calc.penalties(operm.em)
)
try(
registerDoParallel(slurmcore)
operm.uns <- foreach(i=50,
  .combine=c,.packages = "qtl") %dopar% {
    operm <- scantwo(cross.no2, method="em")
}
pen.em <- calc.penalties(operm.uns)
)
print('Done scanning. Saving...')
save.image(paste('QTLmap.Rsave',sep=''))
## probs_map <- interp_genoprob(pr, map)
