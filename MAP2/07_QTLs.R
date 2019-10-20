## NBH

################################################################################
## step-wise
full.norm.add_only <- stepwiseqtl(cross, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=8)
################################################################################

################################################################################
qtl.2 <- makeqtl(cross, chr=c(2), pos=c(80),  what="draws")
out.a.2 <- addqtl(cross, qtl=qtl.2, formula=y~Q1+Q2, method="imp",model="binary",pheno.col=4)
out.i.2 <- addqtl(cross, qtl=qtl.2, formula=y~Q1*Q2, method="imp",model="binary",pheno.col=4)
## indicates possible interaction on 13 and 18

qtl.8 <- makeqtl(cross, chr=c(8), pos=c(88),  what="draws")
out.i.8 <- addqtl(cross, qtl=qtl.8, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=4)
out.a.8 <- addqtl(cross, qtl=qtl.8, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=4)

qtl.18 <- makeqtl(cross, chr=c(18), pos=c(26.7),  what="draws")
out.a.18 <- addqtl(cross, qtl=qtl.18, formula=y~Q1+Q2, method="imp",model="binary",pheno.col=4)
out.i.18 <- addqtl(cross, qtl=qtl.18, formula=y~Q1*Q2, method="imp",model="binary",pheno.col=4)
################################################################################

################################################################################
## manual at stepwise positions
qtl.18 <- makeqtl(cross, chr=c(18), pos=c(27),  what="draws")
out.i.18.b <- addqtl(cross, qtl=qtl.18, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=5)
out.a.18.b <- addqtl(cross, qtl=qtl.18, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=5)

qtl.2 <- makeqtl(cross, chr=c(2), pos=c(80),  what="draws")
out.i.2.b <- addqtl(cross, qtl=qtl.2, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=5)
out.a.2.b <- addqtl(cross, qtl=qtl.2, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=5)

qtl.8 <- makeqtl(cross, chr=c(8), pos=c(67),  what="draws")
out.i.8.b <- addqtl(cross, qtl=qtl.8, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=5)
out.a.8.b <- addqtl(cross, qtl=qtl.8, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=5)
################################################################################

####################################################################################
## Multi/interacting QTL

################################################################################
qtl <- makeqtl(cross, chr=c(2,8,18), pos=c(80,67,27),  what="draws")
fitted <- fitqtl(cross,qtl=qtl, formula=y~Q1+Q2+Q3, pheno.col=4)

qtl <- makeqtl(cross, chr=c(2,8,18), pos=c(80,67,27),  what="draws")
fitted <- fitqtl(cross, qtl=qtl,  formula=y~Q1+Q2+Q3, pheno.col=5, method="imp",model="normal")
fitted.noint.add <- addqtl(cross, qtl=qtl,  formula=y~Q1+Q2+Q3, pheno.col=5, method="imp",model="normal")

qtl <- makeqtl(cross, chr=c(2,18), pos=c(80,27))
fitted.int.only <- fitqtl(cross, qtl=qtl,  formula=y~Q1+Q2+Q1:Q2, pheno.col=5, method="imp",model="normal")
out.i <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2+Q3, pheno.col=5, method="imp",model="normal")

qtl <- makeqtl(cross, chr=c(2,18), pos=c(80,27))
out.ap13 <- addqtl(cross, qtl=qtl, formula = y~Q1 + Q1*Q3 + Q2,  model='normal', method = "imp", pheno.col = 5)
out.ap23 <- addqtl(cross, qtl=qtl, formula = y~Q1 + Q2 + Q2*Q3,  model='normal', method = "imp", pheno.col = 5)

qtl <- makeqtl(cross, chr=c(2,3,13,18,19), pos=c(80,11.18,38.09,27,41.3))
out.ap13 <- fitqtl(cross, qtl=qtl, formula = y~Q1 + Q1*Q2 + Q1*Q3 + Q4 + Q4*Q5,  model='normal', method = "imp", pheno.col = 5)

### 19 interacts w 2 (from picture)
### 18 and 23
####################################################################################
