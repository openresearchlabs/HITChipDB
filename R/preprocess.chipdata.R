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
#'   @param summarization.methods List summarization methods to be included in output. With HITChip frpa always used; with other chips rpa always used. Other options: "sum", "ave"
#'   @param which.projects Optionally specify the projects to extract. All samples from these projects will be included.
#'   @param all.samples Use all samples from the selected project by default? TRUE / FALSE
#'                                        
#' Returns:
#'   @return Preprocessed data and parameters
#'
#' @export
#' @importFrom DBI dbDriver
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

preprocess.chipdata <- function (dbuser, dbpwd, dbname, verbose = TRUE, host = NULL, port = NULL, use.precalculated.phylogeny = NULL, summarization.methods = c("frpa", "sum"), which.projects = NULL, all.samples = TRUE) {

  ## ask parameters or read from R-file
  drv <- dbDriver("MySQL")
  if (!(is.null(host) && is.null(port))) {    
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname, host = host, port = port)
  } else { 
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname)
  }
  
  chip <- detect.chip(dbname)
  params <- ReadParameters(con, which.projects = which.projects, all.samples = all.samples, chip = chip)  

  if (params$chip == "HITChip" && "rpa" %in% summarization.methods)  {
    warning("Frozen-RPA (fRPA) used instead of RPA for HITChip.")
    summarization.methods <- unique(gsub("^rpa", "frpa", summarization.methods))
  } else if (!params$chip == "HITChip" && "frpa" %in% summarization.methods)  {
    warning(paste("RPA used instead of frozen-RPA (fRPA) for", params$chip))
    summarization.methods <- unique(gsub("frpa", "rpa", summarization.methods))

  }

  message("Get sample information matrix for the selected projects")
  project.info <- fetch.sample.info(params$prj$projectName, chiptype = NULL, 
  	       	  		    dbuser = dbuser, 
				    dbpwd = dbpwd, 
				    dbname = dbname, 
				    host = host, 
				    port = port,
  	       	  		    selected.samples = params$samples$sampleID)


  # Let us require that all data is from a unique chip design; otherwise stop				    
  if (length(unique(project.info$designID)) > 1) {
    # message(table(project.info$designID, project.info$projectName))
    stop("The selected projects are from different array versions or missing the array design info! Combining these is not allowed.")
  }

  message("Get probe-level data for the selected hybridisations")
  tmp <- get.probedata(unique(project.info[["hybridisationID"]]), params$rm.phylotypes$oligos, dbuser, dbpwd, dbname, host = host, port = port)
  fdat.orig <- tmp$data       # features x hybs, original non-log scale
  fdat.oligoinfo <- tmp$info  # oligoinfo

  # Annotations for selected hybridisations
  fdat.hybinfo <- project.info[match(colnames(fdat.orig), project.info$hybridisationID), ]
  rownames(fdat.hybinfo) <- colnames(fdat.orig)

  message("Discard the hybs that contain only NAs")
  onlyNA <- colMeans(is.na(fdat.orig)) == 1
  naHybs <- colnames(fdat.orig)[onlyNA]
  if(sum(onlyNA) > 0){
    message("Removing the following hybs, because they contain only NAs:\n")
    message(naHybs,"\n")
    fdat.orig <- fdat.orig[, !onlyNA]
    fdat.hybinfo <- fdat.hybinfo[, !onlyNA]
  }

  ##############################
  ## Between-array normalization
  ##############################

  message("Minmax parameters hard-coded to standardize normalization")
  # Using the parameters from HITChip atlas 
  # params$minmax.points <- c(30.02459, 132616.91371)
  params$minmax.points <- c(30, 133000) 

  # "selected scaling for featurelevel data"
  # Background correction after this step, if any. 
  # Order of normalization / bg correction was validated empirically.
  # bg.adjust intentionally set to NULL 
  # bg correction done _after_ oligo summarization, if any (see next steps)
  d.scaled <- ScaleProfile(fdat.orig, params$normalization, bg.adjust = NULL, minmax.points = params$minmax.points) 

  ##################################
  ## GET OLIGO-PHYLOTYPE MAPPINGS
  ##################################

  message("ReadPhylogeny")
  ph <- ReadPhylogeny(params$phylogeny, params$rm.phylotypes, 
     			params$remove.nonspecific.oligos,
	    		     dbuser = dbuser, 
			     dbpwd = dbpwd, 
			     dbname = dbname, 
			     host = host, 
			     port = port,
			     verbose = verbose, 
			     chip = params$chip)

  taxonomy.filtered <- ph$filtered
  taxonomy.full <- ph$full

  ####################
  ## COMPUTE SUMMARIES
  ####################

  message("Summarize probes into oligos and hybridisations into samples")
  oligo.log10 <- summarize.rawdata(log10(d.scaled), 
  	      			   fdat.hybinfo, 
				   fdat.oligoinfo = fdat.oligoinfo)

  # Return to the original scale
  oligo.abs <- matrix(10^oligo.log10, nrow = nrow(oligo.log10))
  rownames( oligo.abs ) <- rownames( oligo.log10 )
  colnames( oligo.abs ) <- colnames( oligo.log10 )

  levels <- c("species", "L2", "L1")
  if (params$chip %in% c("MITChip", "PITChip", "ChickChip")) {
    levels <- c(levels, "L0")
    if ("frpa" %in% summarization.methods) { 
      warning("fRPA summarization not implemented for PITChip or MITChip - using RPA instead!")
      summarization.methods[[which(summarization.methods == "frpa")]] <- "rpa" 
      summarization.methods <- unique(summarization.methods)
   }
  }

  params$summarization.methods <- summarization.methods

  list(probedata = oligo.abs, taxonomy = taxonomy.filtered, taxonomy.full = taxonomy.full, naHybs = naHybs, params = params)

}


