debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library("qtlbim")


if(pop=='NBH'){

  popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NBH/REMAPS"
  cross.18 <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv",
    sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))
  sex <- read.table(file = file.path(dirso, "sex.txt"))
  rownames(sex) <- sex$ID
  cross.18$pheno$sex <- sex[as.character(cross.18$pheno$ID), 2]
  cross.18$pheno$binary <- as.numeric(cross.18$pheno$pheno>=3)
  crOb <- cross.18
  crOb <- qb.genoprob(crOb, step = 10)
  qbData <- qb.data(crOb, pheno.col = 1, trait = "ordinal", rancov = 2)
  qbModel <- qb.model(crOb, epistasis = T, main.nqtl = 9, mean.nqtl = 10, depen = FALSE)
  mc <- qb.mcmc(crOb, qbData, qbModel, pheno.col = 1, n.iter = 30000)
  so <- qb.scanone(mc, epistasis = T, type = "2logBF")
  so.LPD <- qb.scanone(mc, epistasis = T, type = "LPD")
  scan <- list(upper="main",lower='epistasis')
  st <- qb.scantwo(mc,scan, type.scan="nqtl",chr = c(1, 2,6, 8,13, 18,20, 23,24),epistasis = TRUE)

  png("/home/jmiller1/public_html/nbh.scanone.so.png",width = 1000)
  plot(so, chr = c(1:24),cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=2.5,xlab=NA)
  dev.off()

  png("/home/jmiller1/public_html/nbh.model.png", width = 2000)
  plot(qb.close(mc, target = so))
  #plot(so.LPD, c(1:24))
  dev.off()

  png("/home/jmiller1/public_html/nbh.scanone.st.png", width = 2000, height = 2000)
  plot(st, c(1:24))
  dev.off()

  png("/home/jmiller1/public_html/nbh_interation_plot.png")
  vioplot(as.factor(cross2$geno[[as.character(chr)]][,mark][index]),cross2$pheno[,1][index])
  dev.off()

  ## locus specific

  png("/home/jmiller1/public_html/interation_plot.png")
    effectplot(scon, mname1="2:37728552", mname2="18:18310573")
  dev.off()


  chrp <- c(2,18)
  pos <- summary(so)[as.character(chrp),2]
  qtl <- find.marker(crOb,chr=chrp,pos=pos)
  names(chrp) <- qtl

  fake.f2 <- argmax.geno(crOb, step=5, off.end=5, err=0.1)
  cross2 <- subset(fake.f2, ind=(!is.na(fake.f2$pheno$sex)))
  cross2 <- convert2cross2(cross2)

for (i in 1:length(chrp)){
  chr <- chrp[i]
  mark <- names(chrp)[i]
  index <- cross2$geno[[as.character(chr)]][,mark]>0
  png(paste("/home/jmiller1/public_html/",pop,chr,"_interation_plot.png"))
    plot_pxg(geno=cross2$geno[[as.character(chr)]][,mark][index],pheno=cross2$pheno[,1][index],sort=F,SEmult=2)
  dev.off()

  png(paste("/home/jmiller1/public_html/",pop,chr,"_pxg.png"))
    plotPXG(crOb , mark, pheno.col=1, jitter=1, infer=F)
  dev.off()

  print(mark)
  print(table(cross2$geno[[as.character(chr)]][,mark]))
}




cross2 <- convert2cross2(crOb)
ents <- calc_entropy(genoprob)
 e <- do.call("cbind", ents)
mean_ind <- rowMeans(e)
mean_marker <- colMeans(e)





interp_map(map, oldmap, newmap)

tointerp <- list("7" = c(pos7.1= 5, pos7.2=15, pos7.3=25),
                  "9" = c(pos9.1=20, pos9.2=40))

 interp_map(tointerp, iron$gmap, iron$pmap)


interpPositions(oldpositions, oldmap, newmap)



cross.ELR

oldpositions <- data.frame(chr=   ,pos=   )


genoprob <- calc_genoprob(cross2, map = cross2.map, error_prob = 0.1, map_function = "kosambi", lowmem = FALSE, quiet = TRUE, cores = slurmcore)


inf <- data.frame(chr=as.character(gsub('\\..*:.*','',names(unlist(cross2.map)))), pos= gsub('.*:','',unlist(cross2.map)),sdp=3,snp=gsub('*.\\.','',names(unlist(cross2.map))))


### map position
inf <- data.frame(chr=as.character(gsub('\\..*:.*','',names(unlist(cross2.map)))), pos= gsub('.*:','',unlist(cross2.map)),sdp=3,snp=gsub('*.\\.','',names(unlist(cross2.map))))

psu <- insert_pseudomarkers(cross2.map)

try <- index_snps(cross2.map,inf , tol = 0.01)


try2 <- scan1snps(genoprob, map=cross2.map, pheno=cross2$pheno, keep_all_snps=T, model = "normal",snpinfo =try)



inf <- data.frame(chr=as.character(gsub(':.*','',marker_names(cross2))), pos= unlist(cross2$gmap),sdp=1,snp=marker_names(cross2))


unlist(map1)





try <- scan1snps(genoprob, map=map1,pheno=cross2$pheno, keep_all_snps=T, model = "normal",snpinfo =try)

try <- scan1snps(genoprob, map=map1,pheno=cross2$pheno, keep_all_snps=T, model = "normal")







• Call ‘query_func()’ to grab all SNPs over a region.

     • Use ‘index_snps()’ to group SNPs into equivalence classes.

     • Use ‘genoprob_to_snpprob()’ to convert ‘genoprobs’ to SNP
       probabilities.

     • Use ‘scan1()’ to do a single-QTL scan at the SNPs.

index_snps(map, snpinfo, tol = 0.00000001)



try <- scan1snps(genoprob, map = cross2$gmap, pheno=cross2$pheno, kinship = NULL, addcovar = NULL,
       Xcovar = NULL, model = "normal", query_func=query_variants)

       start = NULL, end = NULL, snpinfo = NULL, batch_length = 20,
       keep_all_snps = FALSE, cores = slurmcore)


try <- scan1snps(genoprob, map = cross2$gmap, pheno=cross2$pheno)

####NOT WORKING
z <- qb.BestPattern(mc, epistasis=F,category = "nqtl",score.type ="distance",include = "all", center = "median", level = 2)

qb.BestPattern(mc)



qb.BayesFactor(mc,nmax=4,cutoff.pattern=0.25, cutoff.pairs = 0.25, epistasis = F)
qb.BayesFactor(mc, items = "pairs",cutoff.pairs = 1,cutoff.pattern=0.25, nmax = 15, epistasis = TRUE)
qb.BayesFactor(mc, items = "pattern",nmax = 4, epistasis = F)
qb.best(mc)
summary(qb.BayesFactor(mc, item = "pattern"))
###########
temp <- summary(qb.scanone(mc, type = "2logBF"), threshold = c(upper = 5))
cross.arch <- qb.arch(temp, chr = c(1,1,2,3), pos = c(15,45,12,15))
try <- step.fitqtl(mc)
plot(slice, figs = c("effects", "cellmean", "effectplot"))
best <- qb.best(mc,category = "nqtl")
