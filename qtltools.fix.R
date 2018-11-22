#!/bin/R
findGenecM <- function(cross, marker.info, gff, gffCols = NULL, gene.id.set = NULL, 
  attributeParse = c("ID="), seqnameParse = c("Chr", "scaffold_"), dropNonColinearMarkers = TRUE, 
  verbose = TRUE, ...) {
  
  dropNonColMar <- function(map) {
    tdiff <- function(y) {
      d1 <- diff(c(y, y[length(y)]))
      d2 <- diff(c(0, y))
      out <- d1 < 0 | d2 < 0
      out[is.na(out)] <- FALSE
      return(out)
    }
    tmp <- map
    tmp$ord <- 1:nrow(tmp)
    good <- unlist(lapply(unique(map$chr), function(x) {
      tc <- tmp[tmp$chr == x, ]
      tc <- tc[order(tc$pos), ]
      d <- tdiff(tc$bp)
      bads <- numeric()
      while (any(d)) {
        bads <- c(bads, tc$ord[d])
        tc <- tc[!d, ]
        d <- tdiff(tc$bp)
      }
      return(tc$ord)
    }))
    return(map[tmp$ord %in% good, ])
  }
  if (dropNonColinearMarkers) {
    marker.info <- dropNonColMar(marker.info)
  }
  
  if (is.null(gffCols) & ncol(gff) != 9) 
    stop("gff file must be of standard format, see details")
  if (!all(c("marker.name", "chr", "pos", "bp") %in% colnames(marker.info))) 
    stop("marker.info must be of standard format, see details")
  
  if (!is.null(gffCols)) {
    if (length(gffCols) != 6) 
      stop("if supplied, gffCols must be a vector length 6")
    if (!all(gffCols %in% colnames(gff))) 
      stop("gffCols must be a vector that matches column names in the gff file")
    gff <- data.frame(chr = gff[, gffCols[1]], source = NA, feature = gff[, gffCols[2]], 
      start = gff[, gffCols[3]], end = gff[, gffCols[4]], score = NA, strand = gff[, 
        gffCols[5]], frame = NA, attribute = gff[, gffCols[6]], stringsAsFactors = F)
  }
  
  gff$geneID <- gene.id.set
  
  
  if (verbose) 
    cat("culling chromosomes to those in the cross\n")
  colnames(gff) <- c("chr", "source", "feature", "start", "end", "score", "strand", 
    "frame", "attribute", "geneID")
  for (i in seqnameParse) {
    gff$chr <- gsub(i, "", gff$chr, fixed = T)
  }
  
  if (is.factor(gff$chr)) {
    gff$chr <- as.character(gff$chr)
  }
  
  gff <- gff[gff$chr %in% as.character(chrnames(cross)), ]
  
  gff$bp <- (gff[, 4] + gff[, 5])/2
  if (verbose) 
    cat("inferring mapping position for:\n")
  out <- lapply(unique(gff$chr), function(i) {
    tgff <- gff[gff$chr == i, ]
    if (verbose) 
      cat("chr ", i, " (n. features = ", nrow(tgff), ")\n", sep = "")
    tmap <- marker.info[marker.info$chr == i, ]
    
    outint <- lapply(2:nrow(tmap), function(j) {
      bpcm <- tmap[(j - 1):j, c("pos", "bp")]
      gffbp <- tgff[tgff$bp >= min(bpcm$bp) & tgff$bp < max(bpcm$bp), ]
      if (nrow(gffbp) >= 1) {
        mod <- lm(pos ~ bp, data = bpcm)
        gffbp$pos = predict(mod, newdata = gffbp)
      }
      return(gffbp)
    })
    outchr <- do.call(rbind, outint)
    if (any(tgff$bp < min(tmap$bp))) {
      outlow <- tgff[tgff$bp < min(tmap$bp), ]
      outlow$pos <- 0
      outchr <- rbind(outlow, outchr)
    }
    if (any(tgff$bp > max(tmap$bp))) {
      outhi <- tgff[tgff$bp >= max(tmap$bp), ]
      outhi$pos <- max(tmap$pos)
      outchr <- rbind(outchr, outhi)
    }
    return(outchr)
  })
  
  return(do.call(rbind, out))
}
