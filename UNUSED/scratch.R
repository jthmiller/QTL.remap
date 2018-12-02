################################
rabbit hole
i <- 1
disto <- sapply(1:24,function(X){
length(which(all.rank$chr==X & all.rank$rank<(max(all.rank$rank)*i)))
})

for (i in seq.int(0.01,0.99,.025)){

  bindo <- sapply(1:24,function(X){
  length(which(all.rank$chr==X & all.rank$rank<(max(all.rank$rank)*i)))
})
  disto <- rbind(disto,bindo)
}

disto <- t(disto)
colnames(disto) <- c(1,seq.int(0.01,0.99,.025))
rownames(disto) <- 1:24

png('/home/jmiller1/public_html/rabithole.png')
plot(as.numeric(colnames(disto)),disto[1,])
sapply(2:24, function(X){ points(as.numeric(colnames(disto)),disto[X,])})
dev.off()

### ln model between position and bp for each chrm
models <- lapply(cross2.map, function(x) {
  pos <- as.numeric(x)
  bp <- as.numeric(gsub("*.:", "", names(x)))
  chr <- data.frame(pos, bp)
  mod <- lm(pos ~ bp, data = chr)
})

## positions to predict CM position for
positions <- lapply(cross2.map, function(x) {
  data.frame(bp = as.numeric(unlist(gsub("*.:", "", names(x)))))
})

predicted <- mapply("predict", models, positions)

AHR.bed$cm.str <- NA
AHR.bed$cm.stp <- NA

for (i in seq_along(AHR.bed$str)) {
  chr <- AHR.bed$chrom[i]
  pos1 <- AHR.bed$str[i]
  pos2 <- AHR.bed$str[i]
  AHR.bed$cm.str[i] <- predict(models[[as.character(chr)]],data.frame(bp=pos1))
  AHR.bed$cm.stp[i] <- predict(models[[as.character(chr)]],data.frame(bp=pos2))
}


  chr$pos <- as.numeric(chr$pos)
  chr$bp <- as.numeric(chr$bp)
  mod <- lm(pos ~ bp, data = chr)
  new <- data.frame(bp = as.numeric(AHR.bed$str[i]))
  AHR.bed$str[i] <- predict(mod, newdata = new)
  new <- data.frame(bp = as.numeric(AHR.bed$stp[i]))
  AHR.bed$stp[i] <- predict(mod, newdata = new)
}



map <- pull.map(cross.nbh)
pos <- gsub('*.\\.*.:','',names(unlist(map)))
elp <- gsub('*.\\.*.:','',names(unlist(pull.map(cross.elr))))
elc <- gsub('\\..*','',names(unlist(pull.map(cross.elr))))
elo <- data.frame(chr=elc,pos=elp)

nep <- gsub('*.\\.*.:','',names(unlist(pull.map(cross.elr))))

elm <-  interpPositions(elo,pull.map(cross.elr),pull.map(cross.nbh))



all <- c(cross.NBH,cross.NEW,cross.ELR)
all <- convert2cross2(all)


all.prob <- calc_genoprob(all, map = NULL, error_prob = 0.1,
       map_function =  "kosambi", lowmem = FALSE,
       quiet = TRUE, cores = 12)

all.max <- maxmarg(all.prob, map = NULL, minprob = 0.95, chr = NULL, pos = NULL,
       return_char = FALSE, quiet = TRUE, cores = 12, tol = 0.0000000000001)


all.max <- est_map(all,maxit=100)
  , map.function="kosambi",n.cluster=12)

try <- scan1(all.max,all$pheno,model = "normal")


sim_geno()


insert_pseudomarkers(map, step = 0, off_end = 0, stepwidth = c("fixed",
     "max"), pseudomarker_map = NULL, tol = 0.01, cores = 1)

nbh <- convert2cross2(cross.NBH)
new <- convert2cross2(cross.NEW)
elr <- convert2cross2(cross.ELR)

all <- c(nbh,new,elr)



probs_to_grid - subset genotype probabilities to a grid of pseudomarkers


map <- insert_pseudomarkers(iron$gmap, step=1)
probs_map <- interp_genoprob(probs, map)











## Scratch

##### Phenotype effect plots

nbh2 <- convert2cross2(NBH.x)
# scan1.nbh <- scan1(genoprobs, pheno, kinship = NULL, addcovar = NULL, Xcovar =
# NULL, intcovar = NULL, weights = NULL, reml = TRUE, model = 'normal', cores =
# slurmcore) model_source_file
out.nbh <- fitqtl(NBH$cross.18, qtl = NBH$qtl.uns)
out.new <- fitqtl(NEW$cross.18, qtl = NEW$qtl.uns)
out.elr <- fitqtl(ELR$cross.18, qtl = ELR$qtl.uns)

nbh2 <- convert2cross2(NBH$cross.18)
elr2 <- convert2cross2(ELR$cross.18)
new2 <- convert2cross2(NEW$cross.18)

est_map


write_control_file



cross.18 <- read.cross.jm(file = file.path(indpops, paste(pop, ".unphased.f2.csvr",
  sep = "")), format = "csvr", geno = c(1:3), estimate.map = FALSE)

### Pull names from plinkfile
path <- file.path(indpops, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross.18$pheno$ID <- paste(popname, indname, sep = "_")

## Subset and drop parents
cross.pars <- subset(cross.18, ind = is.na(cross.18$pheno$Phen))


library(qtl2)
grav2 <- read_cross2("~/my_data/grav2.yaml")




nbh2.map <- est_map(nbh2, error_prob = 0.01, map_function = "kosambi", lowmem = FALSE,
  maxit = 10000, tol = 1e-04, quiet = TRUE, save_rf = FALSE, cores = slurmcore)

nbh_qtl_1 <- find_marker(nbh2.map, NBH$qtl.uns$chr[1], NBH$qtl.uns$pos[1])
nbh_qtl_3 <- find_marker(nbh2.map, NBH$qtl.uns$chr[3], NBH$qtl.uns$pos[3])

# map <- insert_pseudomarkers(cross2$gmap, step = 1)

g <- calc_genoprob(nbh2, map = nbh2.map, error_prob = 0.01, map_function = "kosambi",
  lowmem = FALSE, quiet = TRUE, cores = slurmcore)

q1 <- maxmarg(g, nbh2.map, chr = NBH$qtl.uns$chr[1], pos = NBH$qtl.uns$pos[1], return_char = TRUE)
q2 <- maxmarg(g, nbh2.map, chr = NBH$qtl.uns$chr[3], pos = NBH$qtl.uns$pos[3], return_char = TRUE)


m <- maxmarg(g, minprob = 0.5)
inferg <- predict_snpgeno(g, m)


png("/home/jmiller1/public_html/phenotypes.png", width = 600)
plot_pxg(q1, nbh2$pheno[, 2], pch = 19, cex = 2)
dev.off()

nbh2.imp <- sim_geno(nbh2, n_draws = 200, error_prob = 0.01, map_function = "kosambi",
  lowmem = FALSE, quiet = TRUE, cores = 12)



full <- stepwiseqtl(cross.18, additive.only = T, method = "imp", pheno.col = 6, scan.pairs = T)
out.fq <- fitqtl(cross.18, qtl = qtl.uns, formula = y ~ Q1 + Q2 + Q3 * Q4)
out.fq <- fitqtl(NBH$cross.18, qtl = NBH$qtl.uns)

##### Not needed>?

nbh.grid <- sim.geno(NBH.x, n.draws = 500, step = 5, off.end = 10, error.prob = 0.01,
  map.function = "kosambi", stepwidth = "fixed")
new.grid <- sim.geno(NEW.x, n.draws = 500, step = 5, off.end = 10, error.prob = 0.01,
  map.function = "kosambi", stepwidth = "fixed")
elr.grid <- sim.geno(ELR.x, n.draws = 500, step = 5, off.end = 10, error.prob = 0.01,
  map.function = "kosambi", stepwidth = "fixed")

### Saved to home save.image(file.path('~/NEW.NBH.ELR.QTLmap.Rsave'))

cross.nbh <- subset(cross, ind = cross$pheno$cross1 == 1)
cross.new <- subset(cross, ind = cross$pheno$cross2 == 1)
cross.elr <- subset(cross, ind = cross$pheno$cross3 == 1)

nbh.grid <- reduce2grid(cross.nbh)
new.grid <- reduce2grid(cross.new)
elr.grid <- reduce2grid(cross.elr)

scan.norm.imp.NBH <- scanone(cross.nbh, method = "imp", model = "normal", pheno.col = 6)
scan.norm.imp.NEW <- scanone(cross.new, method = "imp", model = "normal", pheno.col = 6)
scan.norm.imp.ELR <- scanone(cross.elr, method = "imp", model = "normal", pheno.col = 6)








################## not needed? ###
nbh.grid <- reduce2grid(nbh.grid)
new.grid <- reduce2grid(new.grid)
elr.grid <- reduce2grid(elr.grid)

scan.NBH <- scanone(nbh.grid, method = "imp", model = "normal", pheno.col = 6)
scan.NEW <- scanone(new.grid, method = "imp", model = "normal", pheno.col = 6)
scan.ELR <- scanone(elr.grid, method = "imp", model = "normal", pheno.col = 6)

melted.nbh <- data.frame(pop = "NBH", chr = scan.NBH$chr, pos = scan.NBH$pos, lod = scan.NBH$lod)
melted.new <- data.frame(pop = "NEW", chr = scan.NEW$chr, pos = scan.NEW$pos, lod = scan.NEW$lod)
melted.elr <- data.frame(pop = "ELR", chr = scan.ELR$chr, pos = scan.ELR$pos, lod = scan.ELR$lod)
melted <- rbind(melted.nbh, melted.new, melted.elr)

melted2 <- melted[which(!melted$pop == "ELR"), ]
melted2$pop <- factor(melted2$pop, levels = c("NEW", "NBH"))











png("/home/jmiller1/public_html/ggplot2.qtl.png", width = 3000)
plotPheno(nbh2, pheno.col = 1, ...)
dev.off()





cross2 <- convert2cross2(cross2)



### rQTL
cross.18 <- sim.geno(cross.18, error.prob = 0.1, step = 5, n.draws = 250)
scon <- scanone(cross.18, model = "normal", pheno.col = 1, method = "imp", addcovar = cross.18$pheno$sex)

fake.f2 <- argmax.geno(crOb, step = 5, off.end = 5, err = 0.1)
scon <- scanone(fake.f2, model = "binary", pheno.col = 5, method = "imp", addcovar = fake.f2$pheno$sex)

sims <- pull.draws(cross.18)


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

bf <- qb.BayesFactor(mc, items = c("nqtl","chrom","pattern"))
