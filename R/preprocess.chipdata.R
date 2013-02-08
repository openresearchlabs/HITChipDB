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
#'   @param probe.parameters Optional. If probe.parameters are given,
#'          the summarization is based on these and model parameters are not
#' 	    estimated. A list. One element for each probeset with the following probe vectors: 
#'	    affinities, variances
#'                                        
#' Returns:
#'   @return Preprocessed data and parameters
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

preprocess.chipdata <- function (dbuser, dbpwd, dbname, verbose = TRUE, host = NULL, port = NULL, probe.parameters = list()) {

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

  # Minmax parameters hard-coded to standardize normalization;
  # Using the parameters from HITChip atlas with 3200 samples
  # params$minmax.points <- c(30.02459, 132616.91371)
  params$minmax.points <- c(30, 133000) 
  # TODO pre-calculated quantiles encoded here, too?

  # Get sample information matrix for the selected projects	
  project.info <- fetch.sample.info(params$prj$projectName, chiptype = NULL, 
  	       	  		    dbuser, dbpwd, dbname, 
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

  # This handles also pmTm, complement and mismatch filtering
  phylogeny.info.full <- get.phylogeny.info(params$phylogeny, 
	    		     dbuser = dbuser, dbpwd = dbpwd, dbname = dbname, 
			     verbose = verbose, 
			     chip = params$chip)


  # TODO add this in output, too
  phylogeny.info.filtered <- prune16S(phylogeny.info.full, pmTm.margin = 2.5, complement = 1, mismatch = 0, rmoligos = params$rm.phylotypes$oligos, remove.nonspecific.oligos = params$remove.nonspecific.oligos)
  
  phylogeny.info <- phylogeny.info.filtered

			     
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

  # Oligo summarization
  finaldata <- list()
  finaldata[["oligo"]] <- oligo.abs
  levels <- c("species", "L2", "L1")
  if (params$chip == "MITChip" || params$chip == "PITChip") {levels <- c(levels, "L0")}
  for (level in levels) {
    finaldata[[level]] <- list()
    for (method in c("sum", "rpa", "nmf")) {
        message(paste(level, method))
    	summarized.log10 <- summarize.probesets(phylogeny.info,		
			    		  oligo.log10, 
      			       	          method = method, 
					   level = level, 	
				   rm.phylotypes = params$rm.phylotypes,
				probe.parameters = probe.parameters)$summarized.matrix

        # Store the data in absolute scale					
        finaldata[[level]][[method]] <- 10^summarized.log10

    }
  }

  list(data = finaldata, phylogeny.info = phylogeny.info.filtered, phylogeny.info.full = phylogeny.info.full, naHybs = naHybs, params = params)

}

