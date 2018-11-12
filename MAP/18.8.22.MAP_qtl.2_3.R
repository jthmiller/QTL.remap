#!/bin/R

### Map QTLs 2 of 3
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")

cross.18 <- read.cross(format = "csv", dir = popdir, file = paste("chr", X, "_", 
  outname, ".manymarks.QTLmap.csv", sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", 
  "B"))

zero.map <- shiftmap(pull.map(cross.18))
cross.18 <- replacemap(cross.18, zero.map)

print(summary(pull.map(cross.18))[as.character(X), ])

marker.warning()
print("Dropping ~3 markers that inflate the map. Takes a long time...")

# cross.18 <- dropone.par(cross.18, X, drop.its = 3, maxit = 5, map.function =
# 'kosambi', error.prob = 0.03, sex.sp = F, verbose = F, parallel = T, cores =
# slurmcore)


for (i in 1:droppo) {
  
  dropone <- parallel.droponemarker(cross.18, chr = X, error.prob = 0.03, map.function = "kosambi", 
    m = 0, p = 0, maxit = 5, cores = slurmcore, tol = 1e-06, sex.sp = FALSE, 
    verbose = F, parallel = T)
  
  dropone$Ldiff <- abs(as.numeric(dropone$Ldiff))
  drops <- rownames(summary(dropone)[1, ])
  
  cross.18 <- drop.markers(cross.18, drops)
  
  # cross.18 <- qtlTools::dropByDropone(cross = cross.18, droponeRes = dropone,
  # midMarkerThresh = 15, endMarkerThresh = 15, re.est.map = F)
  
  marker.warning()
  
  print(i)
  
}

marker.warning()

if (reorder.marks == T) {
  print("Re-order markers")
  cross.18 <- orderMarkers(cross.18, chr = X, window = 5, use.ripple = T, error.prob = ers, 
    map.function = "kosambi", sex.sp = F, maxit = 3000, tol = 0.001)
}

print("Re-estimating the map")
POS.map.18 <- est.map(cross.18, error.prob = 0.05, map.function = "kosambi", chr = X, 
  maxit = 3000)

cross.18 <- replace.map(cross.18, POS.map.18)

print(summary(pull.map(cross.18))[as.character(X), ])

## Error rate
print("determine error rate for last round of mapping")
ers <- er.rate(cross = cross.18, cpus = slurmcore, maxit = 2000)

print(paste(pop, "error rate for chromosome", X, "is", ers))
system(paste("echo", pop, X, "_", outname, ers, ">> /home/jmiller1/QTL_Map_Raw/popgen/rQTL/remap_out/genotyping_error_rate.txt"))

print("saving...")
save.image(paste(popdir, "/chr", X, "_", outname, ".QTLmap.Rsave", sep = ""))

### Write Map information
system(paste("echo", pop, X, outname, " >> /home/jmiller1/QTL_Map_Raw/popgen/rQTL/remap_out/map.txt", 
  sep = "   "))
line <- unlist(summary(pull.map(cross.18))[as.character(X), ])
write("mar  length  avesp  max", file = "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/remap_out/map.txt", 
  append = TRUE)
write(line, file = "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/remap_out/map.txt", append = TRUE)

### Write cross to file
print("Writing the markers to rQTL format")
write.cross(cross.18, filestem = paste(popdir, "/chr", X, "_", outname, "_2.QTLmap", 
  sep = ""), format = "csv", chr = X)
