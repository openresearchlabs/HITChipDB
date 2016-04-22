

#' list.mysql.projects
#' 
#' Description: List projects in MySQL database
#'
#' Arguments:
#' @param dbuser MySQL user
#' @param dbpwd  MySQL password
#' @param dbname MySqL database name
#' @param host host; needed with FTP connections
#' @param port port; needed with FTP connections

#' Returns:
#'   @return project names vector
#'
#' @export
#' @importFrom DBI dbDriver
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

list.mysql.projects <- function (dbuser, dbpwd, dbname, host = NULL, port = NULL) { 

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
#' @param phylogeny phylogeny (default: 16S)
#' @param dbuser MySQL user
#' @param dbpwd MySQL password
#' @param dbname MySqL database name
#' @param verbose verbose
#' @param chip chip type
#' @param host host; needed with FTP connections
#' @param port port; needed with FTP connections
#' @param rmoligos oligos to exclude
#' @return phylogeny.info
#'
#' @export
#' @importFrom DBI dbDriver
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

get.phylogeny.info <- function (phylogeny = "16S", dbuser, dbpwd, dbname, verbose = TRUE, chip = "HITChip", host = NULL, port = NULL, rmoligos = NULL) {   

  if (verbose) { message("Load phylogeny.info info") }

  drv <- dbDriver("MySQL")
  if (!(is.null(host) && is.null(port))) {
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname, host = host, port = port)
  } else { 
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname)
  }
  
  ## Collect the full phylogenetic information for oligos
  message("Collect the full phylogeny")

  #excloligos <- ifelse(length(rmoligos>0),
  #                     paste('AND NOT (', paste("o.oligoID='",rmoligos,"'",sep="",collapse=" OR "), ') ', sep=''), '')

  if (chip %in% c("MITChip", "PITChip", "ChickChip")) {

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

                        #excloligos,
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
                        #excloligos,
                        "ORDER BY l1.name, l2.name, species.name, specimen.name, o.oligoID;", sep="")

  }

  rs <- dbSendQuery(con, full16Squery)
  full16S <- fetch(rs, n = -1)
  
  # Close MySQL connection
  dbDisconnect(con)

  # Return full 16S phylogeny
  full16S

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

      stm <- paste("SELECT * FROM sample", expandCondition(condition), sep='')
      rs <- dbSendQuery(con, stm)
      smps <- fetch(rs, n=-1)
      return(smps)

   } else {

     stop("Provide proper connection for fetch.samples")

   }
}
