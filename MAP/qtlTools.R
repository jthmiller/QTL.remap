packs <- c('qtl','foreach','doParallel')
lapply(packs, require, character.only = TRUE)
require(qtl2,lib.loc='/share/apps/rmodules')

## For plotting
marker_dens <- list()

### parents ###
cross.pars <- read.cross.jm(file=file.path(indpops,paste(pop,'.parents.csvr',sep='')),
                format='csvr', geno=c(1:3),estimate.map=FALSE)

gt.par <- geno.table(cross.pars)
par.confirm.marks <- rownames(gt.par[which(gt.par$AA==1 & gt.par$BB==1),])

## read in the QTL cross
cross.18 <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
                format='csvr', geno=c(1:3),estimate.map=FALSE)

cross.pars <- subset(cross.pars, chr=c(1:24))

cross.18 <- c(cross.pars,cross.18)

## Map each QTL chro independently
allbut <- c(1:24)[-c(1,2,3,18)]
subset.qtl <- chrnames(cross.18)[!chrnames(cross.18) %in% allbut]
cross.18 <- subset(cross.18, chr=subset.qtl)

## Specific to 'outname'
if(mapped.only==T){
cross.18 <- drop.markers(cross.18,markernames(cross.18)[grep('NW',markernames(cross.18))])
}

#### Filter Conservative
print('Dropping markers with more than 5 genotypes missing')
cross.18 <- drop.missing(cross.18,7)
gt.missing <- geno.table(cross.18)

### invariants
gt.cross.par <- geno.table(cross.18)
cross.18 <- drop.markers(cross.18,rownames(gt.cross.par[gt.cross.par$P.value<cutoff,]))
gt.pval <- geno.table(cross.18)

cross.18 <- formLinkageGroups(cross.18, max.rf=0.1, reorgMarkers=TRUE)
cross.18 <- subset(cross.18, chr=c(1:10))


cross2 <- convert2cross2(cross.18)
map <- insert_pseudomarkers(cross2$gmap, step=1)
pr <- calc_genoprob(cross2, map, err=ers, cores=slurmcore)
pr <- clean_genoprob(pr)
apr <- genoprob_to_alleleprob(pr)

rela <- calc_kinship(pr, type = "overall", omit_x = FALSE,
  use_allele_probs = F, quiet = TRUE, cores = 12)

save.image('kinship.rsave')

############
## Only use genotyped individuals to select markers
#cross <- subset(cross.18,ind=cross.18$pheno$stata=='ind')
#cross <- dropSimilarMarkers(cross,chr=X,rf.threshold = 0.002,byChr = TRUE,re.est.map = FALSE)
#cross <- repRipple(cross, error.prob=ers, map.function="kosambi",window = 6)
#cross <- removeDoubleXO(cross, verbose=T)
#cross.map <- est.map(cross,error.prob=ers,map.function="kosambi",maxit=1000)
#cross <- replace.map(cross, cross.map)
#cross <- est.rf(cross)
#cross.2 <- repRipple(cross, error.prob=ers, map.function="kosambi",window = 6)
#cross.18 <- subset(cross.18,ind=cross.18$pheno$stata=='ind')




cross.18 <- read.cross(format='csv',dir=popdir,
   file=paste('chr',X,'_',outname,'.QTLmap.csv',sep=''),
   geno=c('AA','AB','BB'),alleles=c("A","B"))

cross <- dropSimilarMarkers(cross.18,chr=chrnames(cross)[1],rf.threshold = 0.002,byChr = TRUE,re.est.map = FALSE)
Z <- repRipple(cross.18, error.prob=ers, map.function="kosambi",window = 4)


pdf('~/rep_rip_8.pdf')
plot.rf(cross0, chr = 1, main = "recombination fractions before ripple")
dev.off()








interp_genoprob(probs, map, cores = 1)

NBH.cr <- convert2cross2(NBH$cross.18)
NEW.cr <- convert2cross2(NBH$cross.18)

map.h <- insert_pseudomarkers(NBH.cr$gmap, step=1)
map.w <- insert_pseudomarkers(NEW.cr$gmap, step=1)

pr.h <- calc_genoprob(NBH.cr, map.h, err=0.002, cores=slurmcore)
pr.w <- calc_genoprob(NEW.cr, map.w, err=0.002, cores=slurmcore)

probs_map.h <- interp_genoprob(pr.h, map.h, cores = slurmcore)
probs_map.w <- interp_genoprob(pr.w, map.h, cores = slurmcore)




out_bin <- scan1(probs_map.h,NBH.cr$pheno, model="normal", cores=slurmcore)

find_peaks(out_bin, map.h, threshold=4, drop=1.5)



probs_map.w <- interp_genoprob(pr.w, map.w, cores = slurmcore)








pr <- clean_genoprob(pr)
apr <- genoprob_to_alleleprob(pr)



iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

probs <- calc_genoprob(iron, iron$gmap, error_prob=0.002)

# you generally wouldn't want to do this, but this is an illustration
map <- insert_pseudomarkers(iron$gmap, step=1)
probs_map <- interp_genoprob(probs, map)







### These markers worked in NBH,NEW. Maybe ELR?
marks <- unique(c(markernames(NBH$cross.18),markernames(NBH$cross.18)))

cross.18 <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
                format='csvr', geno=c(1:3),estimate.map=FALSE)

cross.18 <- drop.markers(cross.18,markernames(cross.18)[!markernames(cross.18) %in% marks])

gt <- geno.table(cross.18)
cutoff <- 1.0e-4
cross.18 <- drop.markers(cross.18,rownames(gt[gt$P.value<cutoff,]))
