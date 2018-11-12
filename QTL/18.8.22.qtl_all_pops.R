#!/bin/R

debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library("ggridges")

# load(file.path('~/NEW.NBH.ELR.QTLmap.Rsave'))

NBH <- new.env()
load("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NBH/REMAPS/QTLmap.Rsave", envir = NBH)
ELR <- new.env()
load("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/ELR/REMAPS/QTLmap.Rsave", envir = ELR)
NEW <- new.env()
load("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NEW/REMAPS/QTLmap.Rsave", envir = NEW)

ELR$cross.18$pheno$ID_pop <- "ELR"
NEW$cross.18$pheno$ID_pop <- "NEW"
NBH$cross.18$pheno$ID_pop <- "NBH"

#### Sim and reduce to grid ####
NBH.x <- qtl::clean(NBH$cross.18)
NEW.x <- qtl::clean(NEW$cross.18)
ELR.x <- qtl::clean(ELR$cross.18)

nbh.grid <- sim.geno(NBH.x, n.draws = 500, step = 5, off.end = 10, error.prob = 0.01,
  map.function = "kosambi", stepwidth = "fixed")
new.grid <- sim.geno(NEW.x, n.draws = 500, step = 5, off.end = 10, error.prob = 0.01,
  map.function = "kosambi", stepwidth = "fixed")
elr.grid <- sim.geno(ELR.x, n.draws = 500, step = 5, off.end = 10, error.prob = 0.01,
  map.function = "kosambi", stepwidth = "fixed")

cross <- c(NBH.x, NEW.x, ELR.x)
cross <- sim.geno(cross, n.draws = 20, step = 1, off.end = 10, error.prob = 0.001,
  map.function = "kosambi", stepwidth = "fixed")

POS.map.18 <- est.map(cross, error.prob = 0.05, map.function = "kosambi", maxit = 1000)

### Saved to home {}save.image(file.path('~/NEW.NBH.ELR.QTLmap.Rsave'))
### load(file.path('~/NEW.NBH.ELR.QTLmap.Rsave'))

cross.nbh <- subset(cross, ind = cross$pheno$cross1 == 1)
cross.new <- subset(cross, ind = cross$pheno$cross2 == 1)
cross.elr <- subset(cross, ind = cross$pheno$cross3 == 1)

nbh.grid <- reduce2grid(cross.nbh)
new.grid <- reduce2grid(cross.new)
elr.grid <- reduce2grid(cross.elr)

scan.norm.imp.NBH <- scanone(cross.nbh, method = "imp", model = "normal", pheno.col = 6)
scan.norm.imp.NEW <- scanone(cross.new, method = "imp", model = "normal", pheno.col = 6)
scan.norm.imp.ELR <- scanone(cross.elr, method = "imp", model = "normal", pheno.col = 6)

melted.nbh <- data.frame(pop = "NBH", chr = scan.norm.imp.NBH$chr, pos = scan.norm.imp.NBH$pos,
  lod = scan.norm.imp.NBH$lod)
melted.new <- data.frame(pop = "NEW", chr = scan.norm.imp.NEW$chr, pos = scan.norm.imp.NEW$pos,
  lod = scan.norm.imp.NEW$lod)
melted.elr <- data.frame(pop = "ELR", chr = scan.norm.imp.ELR$chr, pos = scan.norm.imp.ELR$pos,
  lod = scan.norm.imp.ELR$lod)
melted <- rbind(melted.nbh, melted.new, melted.elr)

melted$pop <- factor(melted$pop, levels = rev(c("NBH", "NEW", "ELR")))

# All chr
png("/home/jmiller1/public_html/ggplot2.qtl.png", width = 3000)
ggplot(melted, aes(y = as.factor(pop), x = pos, height = lod, fill = pop)) + ylab(NULL) +
  geom_ridgeline(stat = "identity", scale = 0.2) + facet_wrap(~as.factor(chr),
  scales = "free_x", nrow = 1, ncol = 24, strip.position = NULL) + scale_y_discrete("lod") +
  scale_fill_manual(name = "Populations", guide = F, values = popcol) + theme(axis.text = element_text(size = 10),
  axis.text.y = element_blank())
dev.off()


# model_source_file
out.nbh <- fitqtl(NBH$cross.18, qtl = NBH$qtl.uns)
out.new <- fitqtl(NEW$cross.18, qtl = NEW$qtl.uns)
out.elr <- fitqtl(ELR$cross.18, qtl = ELR$qtl.uns)


# Phenotype plot
pheno.all <- phen <- read.csv("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/metadata/ALL_phenotype_Dist.csv",
  header = T)
phen <- phen[!is.na(phen$pheno_all),]
phen$Pheno_05 <- phen$pheno_all
#index <- which(phen$pop_all == "BRP")

phen$pop_all <- factor(phen$pop_all, levels = c("NBH", "BP", "NEW", "ELR"))

# Phenotype figure
png("/home/jmiller1/public_html/phenotypes.png", width = 600)
ggplot(phen, aes(Pheno_05, y = pop_all, fill = pop_all, group = pop_all, height = ..density..,
  show.legend = F, order = as.numeric(pop_all))) + geom_density_ridges(scale = 1.33,
  bandwidth = 0.5) + xlab("Phenotype Score") + ylab("Density") + scale_alpha(guide = "none") +
  scale_y_discrete(limits = rev(levels(phen$pop_all)), expand = c(0.01, 0.01)) +
  scale_x_continuous(expand = c(0.02, 0.02), limits = c(0, 5)) + scale_fill_manual(name = "Populations",
  guide = "legend", values = popcol) + theme(text = element_text(size = 20), axis.text.y = element_blank(),
  panel.background = element_rect(fill = "white", colour = "black", size = 0.1,
    linetype = 1), )
dev.off()


`?`(plot_pxg)








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


# All chr
png("/home/jmiller1/public_html/ggplot2.qtl.png", width = 3000)
ggplot(melted, aes(y = as.factor(pop), x = pos, height = lod, fill = pop)) + geom_ridgeline(stat = "identity",
  alpha = 0.5, scale = 0.25) + facet_wrap(~as.factor(chr), scales = "free_x", nrow = 1,
  ncol = 24, strip.position = NULL) + # scale_x_continuous('pos') + scale_y_continuous('lod') +
scale_y_discrete("lod") + theme(axis.text = element_text(size = 10))
dev.off()


# again with the qtl2

cross2 <- convert2cross2(cross.18)
map <- insert_pseudomarkers(cross2$gmap, step = 1)
sim_geno <-
nbh2 <- convert2cross2(NBH$cross.18)


plot_pxg(nbh2$geno$"1"[, 1], nbh2$pheno[, 2])

nbh2.imp <- sim_geno(nbh2, n_draws = 200, error_prob = 0.01, map_function = "kosambi",
  lowmem = FALSE, quiet = TRUE, cores = 12)



full <- stepwiseqtl(cross.18, additive.only = T, method = "imp", pheno.col = 6, scan.pairs = T)
out.fq <- fitqtl(cross.18, qtl = qtl.uns, formula = y ~ Q1 + Q2 + Q3 * Q4)
out.fq <- fitqtl(NBH$cross.18, qtl = NBH$qtl.uns)





png("/home/jmiller1/public_html/ggplot2.qtl.png", width = 3000)
plotPheno(nbh2, pheno.col = 1, ...)
dev.off()

# guess_phase - turn imputed genotypes into phased genotypes along chromosomes
# sim_geno - multiple imputations of underlying genotypes given marker data
# predict_snpgeno - predict SNP genotypes in a multiparent population from
# inferred genotypes plus founder strains’ SNP alleles.

perms.unstrat <- scan1perm(pr, cross2$pheno, model = "binary", cores = 12, n_perm = 200)

par(mar = c(4.1, 4.1, 0.6, 0.6))
plot_pxg(nbh2, iron$pheno[, "liver"], ylab = "Liver phenotype")


pxg1 <- find_marker(nbh2, chr = "2", 198)

png("/home/jmiller1/public_html/Q1.png", width = 300)
plotPXG(nbh2, pheno.col = 2)
