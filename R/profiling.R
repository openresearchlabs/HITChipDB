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

#' run.profiling.script
#' 
#' Description: Profiling main script
#'
#' Arguments:
#'   @param dbuser MySQL username
#'   @param dbpwd  MySQL password
#'   @param dbname MySQL database name (HITChip: "Phyloarray"; MITChip: "Phyloarray_MIT"; PITChip old: "Phyloarray_PIT"; PITChip new: "pitchipdb")
#'   @param verbose verbose
#'   @param host host; needed with FTP connections
#'   @param port port; needed with FTP connections
#'   @param summarization.methods List summarization methods to be included in output
#'   @param which.projects Optionally specify the projects to extract. All samples from these projects will be included.
#'
#' Returns:
#'   @return Profiling parameters. Also writes output to the user-specified directory.
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

run.profiling.script <- function (dbuser, dbpwd, dbname, verbose = TRUE, host = NULL, port = NULL, summarization.methods = c("sum", "ave", "nmf", "rpa", "frpa"), which.projects = NULL) {

  # dbuser = "root"; dbpwd = "fidipro"; dbname = "phyloarray"; host = '127.0.0.1'; port = 3307; verbose = T
  
  # Fetch and preprocess the data		     
  chipdata  <- preprocess.chipdata(dbuser, dbpwd, dbname, 
  	       				verbose = verbose,
					   host = host,
					   port = port, 
		 	  summarization.methods = summarization.methods, 
			         which.projects = which.projects)

  finaldata <- chipdata$data
  params    <- chipdata$params
  phylogeny.info <- chipdata$phylogeny.info
  phylogeny.info.full <- chipdata$phylogeny.info.full

  ## Write preprocessed data in tab delimited file
  outd <- WriteChipData(finaldata, params$wdir, phylogeny.info, phylogeny.info.full, verbose = verbose)

  # Add oligo heatmap into output directory
  # Provide oligodata in the _original (non-log) domain_
  hc.params <- add.heatmap(finaldata[["oligo"]], output.dir = params$wdir, phylogeny.info = phylogeny.info)

  # Plot hierachical clustering trees into the output directory
  if (ncol(finaldata[["oligo"]]) > 2) { 

    method <- "complete"
    dat <- finaldata[["oligo"]]

    png(paste(params$wdir, "hclust_oligo_spearman_", method, ".png", sep = ""), width = 480, height = 480 * ncol(dat)/20)
    hc <- hclust(as.dist(1 - cor(dat, use = "pairwise.complete.obs", method = "spearman")), method = method)
    plot(hc, hang = -1, main = "hclust/spearman/oligo/complete", xlab = "Samples")
    dev.off()

  }

  # Plot hclust trees on screen
  tmp <- htree.plot(finaldata[["oligo"]])

  # Write parameters into log file
  tmp <- WriteLog(chipdata$naHybs, params)
  params$logfilename <- tmp$log.file
  params$paramfilename <- tmp$parameter.file

  params

}


