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
