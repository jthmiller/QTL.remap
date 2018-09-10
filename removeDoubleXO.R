
#####

# markersInInterval(cross, chr, min, max)

# returns a list of marker names that fall within a range of locations

# (in cM) along a linkage group (chr). Requires a single linkage group

# output is a character vector

#####



markersInInterval <- function(cross, chr, min, max) {

 names(which(pull.map(cross=cross, chr=chr)[[chr]] < max &

             pull.map(cross=cross, chr=chr)[[chr]] > min))

}



#####

# singleMarkerInInterval(cross, chr, min, max)

# returns the marker name if only a single marker falls in a given

# interval (including the possible case of NAs surrounding marker

# on a linkage group, else returns value FALSE

#####



singleMarkerInInterval <- function(cross, chr, min, max) {

  tmp <- markersInInterval(cross,chr,min,max)  

  val <- ifelse(sum(!is.na(tmp) == 1 | is.na(tmp) == 2), tmp[!is.na(tmp)], FALSE)

}

#####

# removeDoubleXO(cross,chr,verbose)

# Takes a cross object, returns a cross object with all double-recombinant

# genotypes removed.  This is useful to determine if the genotyping errors

# from stacks influence the genetic maps.  if verbose = 1 (TRUE), minimal

# reporting happens, if verbose > 1, all changed genotypes are written to

# STDOUT

#####



removeDoubleXO <- function(cross, chr, verbose=TRUE) { 

  if (!missing(chr)) 

      chr <- matchchr(chr, names(cross$geno))

  else chr <- names(cross$geno)

  for (ch in chr) {  # loop through all linkage groups

    if (verbose)

      cat("Starting linkage group",ch,"\n")



    # find all recombination events in cross. xo is a list spanning

    # individuals.  each element of xo is either NULL or a numeric

    # vector with estimated crossover locations

    

    xo <- locateXO(cross, ch)



    # initialize some variables to keep track of the number

    # of changes

    

    total_removed <- 0

    tot_genotypes <- length(cross$geno[[ch]]$data[,]) -

                          sum(is.na(cross$geno[[ch]]$data[,]))



    for (ind in 1:length(xo)) {            # loop through individuals

      if (length(xo[[ind]]) <= 1) 

        next  # skip individuals with one or fewer recombination events

      # walk along the linkage groups recombination events

      for (location in 3:length(xo[[ind]])-1) {



        # Determine if there is a single marker between each recombination

        # event and the next (location -1 thru location)

        

        sMar <- singleMarkerInInterval(cross,ch,xo[[ind]][location-2],xo[[ind]][location+1])

        if (sMar!=FALSE) { # if there are double recombination events

          oldValue <- cross$geno[[ch]]$data[ind,sMar] # original genotype call

          cross$geno[[ch]]$data[ind,sMar] <- NA # assign genotype to NA

          total_removed <- total_removed + 1 # count removed genotypes

          if (verbose>1)

            cat("individual", ind, "marker", sMar, "oldvalue=", oldValue,

                "newvalue = ", cross$geno[[ch]]$data[ind,sMar], "\n")

        }

      }

    }

    if (verbose) {

      cat("Removed",total_removed,"of",tot_genotypes,

          "genotyped markers on linkage group",ch,"\n")

    }

  }

  cross

}



