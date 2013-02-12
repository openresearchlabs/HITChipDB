# Copyright (C) 2011-2013 Leo Lahti and Jarkko Salojarvi 
# Contact: <microbiome-admin@googlegroups.com>. All rights reserved.

# This file is a part of the microbiome R package
# http://microbiome.github.com/

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' Description: Profiling preprocessing script
#'
#' Arguments:
#'   @param dbuser MySQL username
#'   @param dbpwd  MySQL password
#'   @param dbname MySQL database name
#'   @param verbose monitor processing through intermediate messages
#'   @param host host; needed with FTP connections
#'   @param port port; needed with FTP connections
#'   @param use.precalculated.phylogeny use precalculated phylogeny?
#'                                        
#' Returns:
#'   @return Preprocessed data and parameters
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

preprocess.chipdata <- function (dbuser, dbpwd, dbname, verbose = TRUE, host = NULL, port = NULL, use.precalculated.phylogeny = TRUE) {

  # library(HITChipDB); library(microbiome); fs <- list.files("~/Rpackages/microbiome/HITChipDB/R/", full.names = T); for (f in fs) {source(f)}

  microbiome::InstallMarginal("RMySQL")

  ## ask parameters or read from R-file
  drv <- dbDriver("MySQL")
  if (!(is.null(host) && is.null(port))) {    
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname, host = host, port = port)
  } else { 
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname)
  }  
  
  params <- ReadParameters(con)  
  params$chip <- detect.chip(dbname)
  params$rm.phylotypes <- phylotype.rm.list(params$chip) 
  # List oligos and phylotypes to remove by default

  # Use precalculated phylogeny with HITChip to speed up
  if (is.null(use.precalculated.phylogeny)) {
    if (params$chip == "HITChip")  {use.precalculated.phylogeny <- TRUE}
    if (!params$chip == "HITChip") {use.precalculated.phylogeny <- FALSE}
  }

  # Minmax parameters hard-coded to standardize normalization;
  # Using the parameters from HITChip atlas with 3200 samples
  # params$minmax.points <- c(30.02459, 132616.91371)
  params$minmax.points <- c(30, 133000) 
  # TODO pre-calculated quantiles encoded here, too?

  # Get sample information matrix for the selected projects	
  project.info <- fetch.sample.info(params$prj$projectName, chiptype = NULL, 
  	       	  		    dbuser = dbuser, 
				    dbpwd = dbpwd, 
				    dbname = dbname, 
				    host = host, 
				    port = port,
  	       	  		    selected.samples = params$samples$sampleID)

  message("Get probe-level data for the selected hybridisations")
  tmp <- get.probedata(unique(project.info[["hybridisationID"]]), params$rm.phylotypes$oligos, dbuser, dbpwd, dbname, host = host, port = port)
  fdat.orig <- tmp$data       # features x hybs, original non-log scale
  fdat.oligoinfo <- tmp$info  # oligoinfo

  # Annotations for selected hybridisations
  fdat.hybinfo <- project.info[match(colnames(fdat.orig), project.info$hybridisationID), ]
  rownames(fdat.hybinfo) <- colnames(fdat.orig)

  ## Discard the hybs that contain only NAs
  onlyNA <- colMeans(is.na(fdat.orig)) == 1
  naHybs <- colnames(fdat.orig)[onlyNA]
  if(sum(onlyNA) > 0){
    message("Removing the following hybs, because they contain only NAs:\n")
    message(naHybs,"\n\n")
    fdat.orig <- fdat.orig[, !onlyNA]
    fdat.hybinfo <- fdat.hybinfo[, !onlyNA]
  }
  
  ##############################
  ## Between-array normalization
  ##############################

  # selected scaling for featurelevel data
  # Background correction after this step, if any. 
  # Order of normalization / bg correction was validated empirically.
  # bg.adjust intentionally set to NULL 
  # bg correction done _after_ oligo summarization, if any (see next steps)
  d.scaled <- ScaleProfile(fdat.orig, params$normalization, bg.adjust = NULL, minmax.points = params$minmax.points) 

  ##################################
  ## GET OLIGO-PHYLOTYPE MAPPINGS
  ##################################

  if (!use.precalculated.phylogeny || !chip == "HITChip") {

    message("Fetching Phylogeny from the database")
    phylogeny.full <- get.phylogeny.info(params$phylogeny, 
	    		     dbuser = dbuser, 
			     dbpwd = dbpwd, 
			     dbname = dbname, 
			     host = host, 
			     port = port,
			     verbose = verbose, 
			     chip = params$chip)

    # This handles also pmTm, complement and mismatch filtering
    # This is the phylogeny used in probe summarization into taxonomic levels
    rm.oligos <- sync.rm.phylotypes(params$rm.phylotypes, phylogeny.full)$oligos
    phylogeny.filtered <- prune16S(phylogeny.full, pmTm.margin = 2.5, complement = 1, mismatch = 0, rmoligos = params$rm.phylotypes$oligos, remove.nonspecific.oligos = params$remove.nonspecific.oligos)

    # Keep only relevant cols
    phylogeny.full <- phylogeny.full[, 1:5]; 
    phylogeny.filtered <- phylogeny.filtered[, 1:5]; 

    # Remove duplicate rows
    phylogeny.full <- phylogeny.full[!duplicated(phylogeny.info.full),]
    phylogeny.filtered <- phylogeny.filtered[!duplicated(phylogeny.info.filtered),]

    #write.table(phylogeny.full, file = "~/Rpackages/microbiome/microbiome/inst/extdata/phylogeny.full.tab", sep = "\t", quote = F, row.names = F)
    #write.table(phylogeny.filtered, file = "~/Rpackages/microbiome/microbiome/inst/extdata/phylogeny.filtered.tab", sep = "\t", quote = F, row.names = F)
    #write.table(phylogeny.filtered, file = "~/Rpackages/microbiome/microbiome/inst/extdata/phylogeny.tab", sep = "\t", quote = F, row.names = F)

  } else {

    message("Using pre-calculated phylogeny")
    data.directory <- system.file("extdata/", package = "microbiome")
    phylogeny.filtered <- microbiome::read.profiling(level = "phylogeny", data.dir = data.directory)
    phylogeny.full <- microbiome::read.profiling(level = "phylogeny.full", data.dir = data.directory)
    if (!chip == "HITChip") { stop("Pre-calculated phylogeny available only for HITChip") }
    
  }
  phylogeny.info <- phylogeny.filtered

  ####################
  ## COMPUTE SUMMARIES
  ####################

  # Summarize probes into oligos and hybridisations into samples
  oligo.log10 <- summarize.rawdata(log10(d.scaled), 
  	      			   fdat.hybinfo, 
				   fdat.oligoinfo = fdat.oligoinfo, 
				   oligo.ids = sort(unique(phylogeny.info$oligoID)))

  # Return to the original scale
  oligo.abs <- matrix(10^oligo.log10, nrow = nrow(oligo.log10)) # - 1  
  rownames( oligo.abs ) <- rownames( oligo.log10 )
  colnames( oligo.abs ) <- colnames( oligo.log10 )

  # fs <- list.files("~/Rpackages/microbiome/microbiome/R/", full.names = T); for (f in fs) {source(f)}; fs <- list.files("~/Rpackages/microbiome/HITChipDB/R/", full.names = T); for (f in fs) {source(f)}

  # Oligo summarization
  finaldata <- list()
  finaldata[["oligo"]] <- oligo.abs
  levels <- c("species", "L2", "L1")
  if (params$chip == "MITChip" || params$chip == "PITChip") {levels <- c(levels, "L0")}
  for (level in levels) {
    finaldata[[level]] <- list()
    for (method in c("sum", "rpa", "frpa", "nmf", "ave")) {

        message(paste(level, method))
    	summarized.log10 <- summarize.probesets(
					phylogeny.info = phylogeny.info,		
			    		  oligo.data = oligo.log10, 
      			       	          method = method, 
					   level = level)$summarized.matrix

        # Store the data in absolute scale					
        finaldata[[level]][[method]] <- 10^summarized.log10

    }
  }

  list(data = finaldata, phylogeny.info = phylogeny.info, phylogeny.info.full = phylogeny.full, naHybs = naHybs, params = params)

}
