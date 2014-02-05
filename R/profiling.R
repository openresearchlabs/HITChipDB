# Copyright (C) 2011-2014 Leo Lahti and Jarkko Salojarvi 
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
#'   @param summarization.methods List summarization methods to be included in output. For HITChip frpa always used; for other chips, rpa always used. Other options: sum, ave, nmf.
#'   @param which.projects Optionally specify the projects to extract. All samples from these projects will be included.
#'
#' Returns:
#'   @return Profiling parameters. Also writes output to the user-specified directory.
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

run.profiling.script <- function (dbuser, dbpwd, dbname, verbose = TRUE, host = NULL, port = NULL, summarization.methods = c("frpa", "sum"), which.projects = NULL) {

  # Fetch and preprocess the data		     
  chipdata  <- preprocess.chipdata(dbuser, dbpwd, dbname, 
  	       				verbose = verbose,
					   host = host,
					   port = port, 
		 	  summarization.methods = summarization.methods, 
			         which.projects = which.projects)

  finaldata <- chipdata$data
  params    <- chipdata$params

  # Phylogeny used for L1/L2/species summarization
  phylogeny.info <- chipdata$phylogeny.info
  # Complete phylogeny before melting temperature etc. filters
  phylogeny.info.full <- chipdata$phylogeny.full

  ## Write preprocessed data in tab delimited file
  outd <- WriteChipData(finaldata, params$wdir, phylogeny.info, phylogeny.info.full, verbose = verbose)

  # Add oligo heatmap into output directory
  # Provide oligodata in the _original (non-log) domain_
  hc.params <- add.heatmap(log10(finaldata[["oligo"]]), output.dir = params$wdir, phylogeny.info = phylogeny.info)

  # Plot hierachical clustering trees into the output directory
  if (ncol(finaldata[["oligo"]]) > 2) { 

    method <- "complete"
    dat <- finaldata[["oligo"]]

    # Clustering
    hc <- hclust(as.dist(1 - cor(log10(dat), use = "pairwise.complete.obs", method = "pearson")), method = method)

    # Save into file
    png(paste(params$wdir, "/hclust_oligo_pearson_", method, ".png", sep = ""), height = 480, width = 480 * ncol(dat)/20)
    plot(hc, hang = -1, main = "hclust/pearson/oligo/log10/complete", xlab = "Samples", ylab = "1 - Correlation")
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









#' Description: Default list of removed phylotypes and oligos
#'
#' Arguments:
#'  @param chip Chip name (HIT/MIT/PIT/Chick)Chip
#' Returns:
#'   @return List of removed oligos and phylotypes
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

phylotype.rm.list <- function (chip) {

  rm.phylotypes <- list()

  if (chip == "HITChip") {
    
    rm.phylotypes[["oligos"]] <- c("UNI 515", "HIT 5658", "HIT 1503", "HIT 1505", "HIT 1506")
    rm.phylotypes[["species"]] <- c("Victivallis vadensis")
    rm.phylotypes[["L1"]] <- c("Lentisphaerae")
    rm.phylotypes[["L2"]] <- c("Victivallis")

  } else if (chip == "MITChip") {

    rm.phylotypes[["oligos"]] <- c("Bacteria", "DHC_1", "DHC_2", "DHC_3", "DHC_4", "DHC_5", "DHC_6", "Univ_1492")
    rm.phylotypes[["species"]] <- c()
    rm.phylotypes[["L1"]] <- c()
    rm.phylotypes[["L2"]] <- c()

  } else if (chip == "PITChip") {

    # Based on JZ mail 9/2012; LL

    rm.old.oligos <- c("Bacteria", "DHC_1", "DHC_2", "DHC_3", "DHC_4", "DHC_5", "DHC_6", "Univ_1492")
    rm.new.oligos <- c("PIT_1083", "PIT_1022", "PIT_1057", "PIT_1023", "PIT_1118", "PIT_1040", "PIT_1058", "PIT_1119", "PIT_122", "PIT_1221", "PIT_1322", "PIT_1367", "PIT_1489", "PIT_160", "PIT_1628", "PIT_1829", "PIT_1855", "PIT_1963", "PIT_1976", "PIT_1988", "PIT_2002", "PIT_2027", "PIT_2034", "PIT_2101", "PIT_2196", "PIT_2209", "PIT_2281", "PIT_2391", "PIT_2392", "PIT_2418", "PIT_2425", "PIT_2426", "PIT_2498", "PIT_2555", "PIT_2563", "PIT_2651", "PIT_2654", "PIT_2699", "PIT_2741", "PIT_2777", "PIT_2786", "PIT_2936", "PIT_35", "PIT_425", "PIT_427", "PIT_428", "PIT_429", "PIT_435", "PIT_481", "PIT_605", "PIT_7", "PIT_733", "PIT_734", "PIT_892")
    rm.phylotypes[["oligos"]] <- c(rm.old.oligos, rm.new.oligos)
    rm.phylotypes[["species"]] <- c()
    rm.phylotypes[["L0"]] <- c("Nematoda", "Apicomplexa", "Euryarchaeota", "Ascomycota", "Parabasalidea", "Chordata")
    rm.phylotypes[["L1"]] <- c("Chromadorea", "Coccidia", "Methanobacteria", "Saccharomycetales", "Trichomonada", "Mammalia")
    rm.phylotypes[["L2"]] <- c("Ascaris suum et rel.", "Eimeria  et rel.", "Methanobrevibacter et rel.", "Saccharomyces et rel.", "Trichomonas et rel.", "Uncultured Mammalia", "Uncultured methanobacteria")

  } else if (chip == "ChickChip") {
    warning("No universal probes excluded from ChichChip yet!")
  }

  rm.phylotypes

}

