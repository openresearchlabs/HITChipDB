# Copyright (C) 2011-2012 Leo Lahti and Jarkko Salojarvi 
# Contact: <microbiome-admin@googlegroups.com>. All rights reserved.

# This file is a part of the microbiome R package
# http://microbiome.github.com/

# This program is open source software; you can redistribute it and/or
# modify it under the terms of the FreeBSD License (keep this notice):
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' fetch.sample.info
#'
#' Description: Fetch sample information from HITChip atlas
#'
#' Arguments:
#'   @param allowed.projects list projects for which to fetch the data
#'   @param chiptype chiptype (eg. new.chip)
#'   @param dbuser MySQL user
#'   @param dbpwd MySQL password
#'   @param dbname MySqL database name
#'   @param selected.samples Sample to investigate. By default all.
#'   @param host host; needed with FTP connections
#'   @param port port; needed with FTP connections
#' Returns:
#'   @return project.info data.frame
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

fetch.sample.info <- function (allowed.projects, chiptype = NULL, 
		  dbuser, dbpwd, dbname, selected.samples = NULL, host = NULL, port = NULL) { 

  if (!require(RMySQL)) {
    install.packages("RMySQL")
    require(RMySQL)
  }

  drv <- dbDriver("MySQL")
  if (!(is.null(host) && is.null(port))) {
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname, host = host, port = port)
  } else { 
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname)
  }  


  # Fetch all data from the database
   rs <- dbSendQuery(con, paste("SELECT p.projectName,p.projectID,s.subjectID,s.sampleID,s.samplingDate,s.normAlgVersion,h.hybridisationID,h.dye,a.arrayID,a.barcode,sl.designID,s.reproducibility,s.normalisationFinished,s.imageID,fe.extractionID,fe.extractionName,fe.noSampleNormalisation,h.isDiscarded,fe.hasReproCheck
     FROM sample s               
     JOIN hybridisation h USING (sampleID) JOIN featureextraction fe USING (hybridisationID)
     JOIN project p USING (projectID)
     JOIN array a USING (arrayID)
     JOIN slide sl USING (barcode)
     ORDER BY s.projectID, s.sampleID, h.hybridisationID, fe.extractionID"))

  message("Fetch selected projects and samples")
  project.info.all <- fetch(rs, n = -1)

  # If no chiptype specified, use all
  if (is.null(chiptype)) {chiptype <- unique(project.info.all$designID)}
  if (is.null(selected.samples)) {selected.samples <- unique(project.info.all$sampleID)}

  # Close MySQL connection
  dbDisconnect(con) 

  # Filter out samples based on predefined criteria
  rkeep <- project.info.all$projectName %in% allowed.projects &
           !as.logical(project.info.all$isDiscarded) &
           !as.logical(project.info.all$noSampleNormalisation) &
           as.logical(project.info.all$normalisationFinished) &
           project.info.all$hasReproCheck & 
           project.info.all$designID %in% chiptype &
           project.info.all$normAlgVersion == 1.1 &
    	   project.info.all$sampleID %in% selected.samples
 
  # Remove annotations which are identical for all samples
  ckeep <- sapply(project.info.all, function (x) {!length(unique(x)) == 1})

  message("Filter the data")
  project.info <- project.info.all[rkeep, ckeep]

  project.info                  
     
}

#' list.mysql.projects
#' 
#' Description: List projects in MySQL database
#'
#' Arguments:
#'   @param dbuser MySQL user
#'   @param dbpwd  MySQL password
#'   @param dbname MySqL database name
#'   @param host host; needed with FTP connections
#'   @param port port; needed with FTP connections

#' Returns:
#'   @return project names vector
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

list.mysql.projects <- function (dbuser, dbpwd, dbname, host = NULL, port = NULL) { 

  if (!require(RMySQL)) {
    install.packages("RMySQL")
    require(RMySQL)
  }

  drv <- dbDriver("MySQL")

  if (!(is.null(host) && is.null(port))) {
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname, host = host, port = port)
  } else { 
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname)
  }
  
  # Fetch all data from the database
  rs <- dbSendQuery(con, paste("SELECT p.projectName FROM project p"))

  project.info <- fetch(rs, n = -1) 

  unique(project.info$projectName)

}



#' get.phylogeny.info
#' 
#' Description: Get phylogeny
#' 
#' Arguments:
#'   @param phylogeny phylogeny (default: 16S)
#'   @param rmoligos oligos to exclude
#'   @param dbuser MySQL user
#'   @param dbpwd MySQL password
#'   @param dbname MySqL database name
#'   @param verbose verbose
#'   @param remove.nonspecific.oligos Logical. Remove oligos with multiple targets.
#'   @param chip chip type
#'   @param host host; needed with FTP connections
#'   @param port port; needed with FTP connections
#' Returns:
#'   @return phylogeny.info
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

get.phylogeny.info <- function (phylogeny = "16S", rmoligos = NULL, dbuser, dbpwd, dbname, verbose = TRUE, remove.nonspecific.oligos = FALSE, chip = "HITChip", host = NULL, port = NULL) {   

  # phylogeny = "16S"; rmoligos = NULL; verbose = TRUE; remove.nonspecific.oligos = FALSE; chip = "HITChip"

  if (!require(RMySQL)) {
    install.packages("RMySQL")
    require(RMySQL)
  }

  if (verbose) { message("Load phylogeny.info info") }

  require(RMySQL)
  drv <- dbDriver("MySQL")
  if (!(is.null(host) && is.null(port))) {
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname, host = host, port = port)
  } else { 
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname)
  }
  
  ## Collect the full phylogenetic information for oligos
  message("Collect the full phylogeny")
  excloligos <- ifelse(length(rmoligos>0),
                       paste('AND NOT (', paste("o.oligoID='",rmoligos,"'",sep="",collapse=" OR "), ') ', sep=''), '')
                       

  if (chip == "MITChip" || chip == "PITChip") {

    # Also get level0

    full16Squery <- paste("SELECT 		   		 
   		   		 l0.name AS 'L0', 
   		                 l1.name AS 'L1', ", "                        	 
   		                 l2.name AS 'L2', ", "                        	 
				 species.name AS 'species', 
				 specimen.name AS 'specimen', 
				 ot.oligoID AS 'oligoID', ","
                        	 o.pmTm, 
				 ot.Tm, 
				 ot.mismatch, 
				 ot.complement ",

                        "FROM phylogeny ph ",

                        "JOIN taxon l0 USING (phylogenyID) ",
                        "JOIN taxtotax tt0 ON (tt0.parentID=l0.taxonID) ",

                        "JOIN taxon l1 ON (tt0.childID=l1.taxonID) ",
                        "JOIN taxtotax tt1 ON (tt1.parentID=l1.taxonID) ",

                        "JOIN taxon l2 ON (tt1.childID=l2.taxonID) ",
                        "JOIN taxtotax tt2 ON (tt2.parentID=l2.taxonID) ",

                        "JOIN taxon species ON (tt2.childID=species.taxonID) ",
                        "JOIN taxtotax tt3 ON (tt3.parentID=species.taxonID) ",
                        "JOIN taxon specimen ON (tt3.childID=specimen.taxonID) ",
                        "JOIN oligotargetpair ot ON (ot.targetID=specimen.targetID) ",
                        "JOIN oligo o ON (ot.oligoID=o.oligoID) ",

                        "AND l0.taxonLevel='level 0' ",

                        "AND tt0.nodeDistance=1 ",
                        "AND tt1.nodeDistance=1 ",
                        "AND tt2.nodeDistance=1 ",
                        "AND tt3.nodeDistance=1 ",

                        excloligos,
                        "ORDER BY l0.name, l1.name, l2.name, species.name, specimen.name, o.oligoID;", sep="")
  } else {

    # get levels 1-3
    full16Squery <- paste("SELECT l1.name AS 'L1', l2.name AS 'L2', ", 
                        "species.name AS 'species', specimen.name AS 'specimen', ot.oligoID AS 'oligoID', ",
                        "o.pmTm, ot.Tm, ot.mismatch, ot.complement ",
                        "FROM phylogeny ph ",
                        "JOIN taxon l1 USING (phylogenyID) ",
                        "JOIN taxtotax tt1 ON (tt1.parentID=l1.taxonID) ",
                        "JOIN taxon l2 ON (tt1.childID=l2.taxonID) ",
                        "JOIN taxtotax tt2 ON (tt2.parentID=l2.taxonID) ",
                        "JOIN taxon species ON (tt2.childID=species.taxonID) ",
                        "JOIN taxtotax tt3 ON (tt3.parentID=species.taxonID) ",
                        "JOIN taxon specimen ON (tt3.childID=specimen.taxonID) ",
                        "JOIN oligotargetpair ot ON (ot.targetID=specimen.targetID) ",
                        "JOIN oligo o ON (ot.oligoID=o.oligoID) ",
                        "AND l1.taxonLevel='level 1' ",
                        "AND tt1.nodeDistance=1 ",
                        "AND tt2.nodeDistance=1 ",
                        excloligos,
                        "ORDER BY l1.name, l2.name, species.name, specimen.name, o.oligoID;", sep="")

  }

  rs <- dbSendQuery(con, full16Squery)
  full16S <- fetch(rs, n = -1)
  
  # Close MySQL connection
  dbDisconnect(con)

  # Apply the standard filters
  if (verbose) {
    message("Hybridisation temperature: Tm >= pmTm - 2.5\n")
    message("Number of mismatches: 0\n")
    message("Must be a complement sequence\n")
    message("No requirement for a full-length hybridisation\n\n")
  }

  phylogeny.info <- prune16S(full16S, pmTm.margin = 2.5, complement = 1, mismatch = 0)

  rmoligos2 <- rmoligos
  if (remove.nonspecific.oligos) {
    if (verbose) {message("Removing oligos that have multiple targets at L2 level")}
    nPhylotypesPerOligo <- n.phylotypes.per.oligo(phylogeny.info, "L2") 
    nonspecific.oligos <- setdiff(phylogeny.info$oligoID, names(which(nPhylotypesPerOligo == 1)))
    rmoligos2 <- c(rmoligos, nonspecific.oligos)
  } 

  phylogeny.info <- phylogeny.info[!phylogeny.info$oligoID %in% rmoligos2, ]

  phylogeny.info        

}



#' fetch.projects
#' 
#' Fetch projects from the phyloarray database
#' @param con MySQL connection
#' @param condition list of lists with field-value pairs
#'
#' Fetch complete records from the \emph{projects}
#' table. If \code{condition=NULL}; will fetch all records and if
#' \code{condition} is defined it will fetch only those records that
#' comply with the condition.  Condition must be a list of lists, each
#' having at least the fields \emph{field} and \emph{value}.  The
#' \emph{field} field must be a character vector of length 1 and
#' \emph{value} must be a vector.  Each of the values will be evaluated
#' as a optional value for the \emph{field} field.  This is an example:
#' "\code{list(list(field='projectName',value=c(A,B)))}" that will be
#' evaluated to the SQL condition "\code{WHERE (projectName='A' OR
#' projectName='B')}".  
#' @return A dataframe with the selected records from the projects table
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples # TBA
#' @keywords utilities

fetch.projects <- function (con, condition = NULL) {

   if (phyloarrayConnection(con)) {

    if (!require(RMySQL)) {
      install.packages("RMySQL")
      require(RMySQL)
    }

      stm <- paste("SELECT * FROM project", expandCondition(condition), sep='')
      rs <- dbSendQuery(con, stm)
      prjs <- fetch(rs, n=-1)
      return(prjs)
   } else {
     stop("Provide proper connection for fetch.projects")
   }
}





#' Fetch samples from a phyloarray database
#'
#' The function fetches complete records from the \emph{samples}
#' table. If \code{condition=NULL} it will fetch all records and if
#' \code{condition} is defined it will fetch only those records that
#' comply with the condition.  Condition must be a list of lists, each
#' having at least the fields \emph{field} and \emph{value}.  The
#' \emph{field} field must be a character vector of length 1 and
#' \emph{value} must be a vector.  Each of the values will be evaluated
#' as a optional value for the \emph{field} field.  This is an example:
#' "\code{list(list(field='projectName',value=c(A,B)))}" that will be
#' evaluated to the SQL condition "\code{WHERE (projectName='A' OR
#' projectName='B')}".
#'
#' @param con MySQL connection
#' @param condition list of lists with field-value pairs
#'
#' @return A dataframe with the selected records from the samples table
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples # TBA
#' @keywords utilities

fetch.samples <- function (con, condition = NULL) {
   if (phyloarrayConnection(con)) {

      if (!require(RMySQL)) {
        install.packages("RMySQL")
        require(RMySQL)
      }

      stm <- paste("SELECT * FROM sample", expandCondition(condition), sep='')
      rs <- dbSendQuery(con, stm)
      smps <- fetch(rs, n=-1)
      return(smps)
   } else {
     stop("Provide proper connection for fetch.samples")
   }
}
