source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')

### Rsave with mapped markers
load(paste('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/',pop,'/REMAPS/QTLmap.Rsave',sep=''))

##### Re-do kinship analysis with markers that mapped
##### Drop idividuals that are to closely/distantly related and remap everything

cross <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
        format='csvr', geno=c(1:3),estimate.map=FALSE)

cross.18 <- subset(cross.18,ind=cross.18$pheno$stata=='ind')
path <- file.path(indpops,paste(pop,'.ped',sep=''))
popname <- system(paste('cut -f1 -d\' \'',path),intern = TRUE)
indname <- system(paste('cut -f2 -d\' \'',path),intern = TRUE)
cross$pheno$ID <- paste(popname,indname,sep='_')

drop <- markernames(cross)[!markernames(cross) %in% markernames(cross.18)]
cross  <- drop.markers(cross,drop)

### Taken from Karl Broman rQTL discussion
### https://groups.google.com/forum/#!searchin/rqtl-disc/Marker$20names$20don$27t$20match%7Csort:date/rqtl-disc/wcZZh0lfTiI/PLcpTV42yjQJ
# locations of x2 markers in cross x1
pos <- find.markerpos(cross.18, markernames(cross))
# for those found, move to position in x1
for(i in which(!is.na(pos[,1]))) cross <- movemarker(cross, rownames(pos)[i], pos[i,1], pos[i,2])
# perhaps drop the markers that weren't in x1
cross <- drop.markers(cross, rownames(pos)[is.na(pos[,1])])
# probably need to sort the chromosomes
cross$geno <- cross$geno[names(cross.18$geno)]
# use map to replace parents map
map <- pull.map(cross.18)
cross <- replace.map(cross,map)

### Lowdata
cross.min <- subset(cross, ind=nmissing(cross)>median(nmissing(cross)) | is.na(cross$pheno$Pheno))
### Drop ind with greater than 50% missing data
cross.max <- subset(cross, ind=nmissing(cross)<nmissing(cross)[1:(round(sum(nind(cross))/2,digits=0))] | is.na(cross$pheno$Pheno))

## Calculate matrix
rela <- rels(cross.max)
name <- paste(pop,'_kinship_mapped_markers.pdf',sep='')
main <- paste(pop,'kinship of good markers (proportion of shared genotypes, 0-1)')
feet(rela,name,main)

rela <- rels(cross.min)
name <- paste(pop,'_kinship_lowcov_mapped_markers.pdf',sep='')
main <- paste(pop,'kinship low cov w/parents (proportion of shared genotypes, 0-1)')
feet(rela,name,main)

save.image(paste('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/',pop,'/REMAPS/kinship_QTLmap.Rsave',sep=''))

#### All pops ####

cross.elr <- mega.cross('ELR')
cross.nbh <- mega.cross('NBH')
cross.new <- mega.cross('NEW')

marks <- intersect(intersect(markernames(cross.nbh),markernames(cross.new)),markernames(cross.elr))

cross.new <- drop.markers(cross.new,markernames(cross.new)[markernames(cross.new)%in%marks])
cross.elr <- drop.markers(cross.elr,markernames(cross.elr)[markernames(cross.elr)%in%marks])
cross.nbh <- drop.markers(cross.nbh,markernames(cross.nbh)[markernames(cross.nbh)%in%marks])

mega.cross <- c(cross.new,cross.elr,cross.nbh)
pop <- NA
rela <- rels(mega.cross)
name <- 'multipop.pdf'
main <- 'ELR,NEW,NBH w/parents (proportion of shared genotypes, 0-1)'
newt(rela,name,main,'~/')
name <- 'multipop_mds.pdf'
mdees(rela,name,main,'~/')

mega.cross <- c(cross.new,cross.elr)
pop <- NA
rela <- rels(mega.cross)
name <- 'ELR_NEW_forcedir.pdf'
main <- 'ELR,NEW w/parents (proportion of shared genotypes, 0-1)'
newt(rela,name,main,'~/')
name <- 'elr_new_mds.pdf'
mdees(rela,name,main,'~/')

mega.cross <- c(cross.nbh,cross.elr)
pop <- NA
rela <- rels(mega.cross)
name <- 'ELR_NBH_forcedir.pdf'
main <- 'ELR,NBH w/parents (proportion of shared genotypes, 0-1)'
newt(rela,name,main,'~/')
name <- 'ELR_NBH_mds.pdf'
mdees(rela,name,main,'~/')

pop <- 'NEW'
rela <- rels(cross.new)
name <- 'NEW_force_directed.pdf'
main <- 'NEW,w/parents (proportion of shared genotypes, 0-1)'
newt(rela,name,main,'~/')
name <- 'NEW_mds.pdf'
mdees.single(rela,name,main,'~/')

pop <- 'ELR'
rela <- rels(cross.elr)
name <- 'ELR_force_directed.pdf'
main <- 'ELR,w/parents (proportion of shared genotypes, 0-1)'
newt(rela,name,main,'~/')
name <- 'ELR_mds.pdf'
mdees.single(rela,name,main,'~/')

pop <- 'NBH'
rela <- rels(cross.nbh)
name <- 'NBH_forcedir.pdf'
main <- 'NBH w/parents (proportion of shared genotypes, 0-1)'
newt(rela,name,main,'~/')
name <- 'NBH_mds.pdf'
mdees(rela,name,main,'~/')

### IBD calculations https://www.cog-genomics.org/plink/1.9/ibd
ids <- read.table('/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops/ELR.mdist.id')
n.len <- file.info('/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops/ELR.mdist.bin')$size
dis <- readBin('/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops/ELR.mdist.bin',what='double',n=n.len)
#con <- file(file.path('~/Dropbox/QTL_Paper/Rough Figures/Kinship Analysis/NBH.kinship.keep.ind.txt'), open='r')

### Sim
ids <- read.table('/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops/NEW.mibs.id')
n.len <- file.info('/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops/NEW.mibs.bin')$size
dis <- readBin('/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops/NEW.mibs.bin',what='double',n=n.len)

dis <- as.numeric(dis)
dis <- matrix(dis, nrow=length(ids$V1),ncol=length(ids$V1))
colnames(dis) <- paste(ids$V1,ids$V2,sep="_")
rownames(dis) <- paste(ids$V1,ids$V2,sep="_")

d <- as.dist(dis)
#fit <- cmdscale(dis,eig=TRUE, k=2)


name <- 'NBH_force_directed_plink.pdf'
main <- 'NBH,w/parents (proportion of shared genotypes, 0-1)'
newt(dis,name,main,'~/')
name <- 'NBH_mds_plink.pdf'
mdees.single.IBS(dis,name,main,'~/',dist=F)

ids <- read.table('/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops/ALL.mibs.id')
n.len <- file.info('/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops/ALL.mibs.bin')$size
dis <- readBin('/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops/ALL.mibs.bin',what='double',n=n.len)
dis <- as.numeric(dis)
dis <- matrix(dis, nrow=length(ids$V1),ncol=length(ids$V1))
colnames(dis) <- paste(ids$V1,ids$V2,sep="_")
rownames(dis) <- paste(ids$V1,ids$V2,sep="_")

dis <- dis[rowSums(is.na(dis)) != (ncol(dis)-1),colSums(is.na(dis)) != (nrow(dis)-1) ]

pop <- 'ELR'
name <- 'ALL_force_directed_plink2.pdf'
main <- 'All Offspring ,w/parents'
newt(dis,name,main,'~/')
name <- 'ALLPOP_mds_plink.pdf'
mdees.single.IBS(dis,name,main,'~/')






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
