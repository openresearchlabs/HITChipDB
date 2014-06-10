#' Description: Get probedata
#' 
#' Arguments:
#'   @param hybridization.ids Specify the hybridizations to retrieve
#'   @param rmoligos oligos to exclude
#'   @param dbuser MySQL user
#'   @param dbpwd  MySQL password
#'   @param dbname MySqL database name
#'   @param host host; needed with FTP connections
#'   @param port port; needed with FTP connections
#' Returns:                                        
#'   @return list with data (features x hybridizations matrix) and info (features x info) fields 
#'
#' @export
#' @import RMySQL
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

get.probedata <- function (hybridization.ids, rmoligos, dbuser, dbpwd, dbname, host = NULL, port = NULL) {

  #hybridization.ids <- unique(project.info[["hybridisationID"]]); rmoligos <- params$rm.phylotypes$oligos

  microbiome::InstallMarginal("RMySQL")

  # List unique hybridisations for the selected samples
  hids <- mysql.format(hybridization.ids)
                      
  drv <- dbDriver("MySQL")
  if (!(is.null(host) && is.null(port))) {
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname, host = host, port = port)
  } else { 
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname)
  }  
  
  # JS 16.5. correction for duplicate extractionIDs under same hybridisationID:
  # choose the one which gets picked to sample normalisation (set during manual quality control of the data):
  rs <- dbSendQuery(con, statement = paste("SELECT featureID,extractionID,fe.hybridisationID,spatNormSignal,isOutlier
      		FROM featuremeasurement 
		JOIN featureextraction fe USING (extractionID)
		JOIN hybridisation h USING (hybridisationID)
                JOIN arrayfeature af USE INDEX (PRIMARY) USING (featureID)
		WHERE fe.hybridisationID IN", hids,"and NOT fe.noSampleNormalisation"))
  rawdata <- fetch(rs, n = -1)

  ## Check if there is any data
  if(nrow(rawdata) == 0) {
    stop("No data found for these samples (perhaps they are not normalized yet?).\n\n")
  }

  message("Remove outliers")
  rawdata$spatNormSignal[as.logical(rawdata$isOutlier)] <- NA

  message("Split data into arrays")
  rawdata.esplit <- split(rawdata, rawdata$hybridisationID)

  message("Remove NAs")
  na.inds <- sapply(rawdata.esplit, function (x) all(is.na(x$spatNormSignal)))
  rawdata.esplit <- rawdata.esplit[!na.inds]

  # Get probeID - featureID - oligoID mappings
  rs <- dbSendQuery(con, "SELECT fe.featureID,p.probeID,p.oligoID,fe.arrayCol,fe.arrayRow FROM arrayfeature fe JOIN probe p USING (probeID)")
  probes <- fetch(rs, n = -1) 
  
  # Remove specified oligos
  probes <- probes[!probes$oligoID %in% rmoligos,]

  ftab.info <- data.frame(list(featureID = unique(rawdata$featureID)))
  ftab.info[["probeID"]] <- probes$probeID[match(ftab.info$featureID, probes$featureID)]
  ftab.info[["oligoID"]] <- probes$oligoID[match(ftab.info$probeID, probes$probeID)]

  message("Remove NA oligos")
  keep <- !is.na(ftab.info$oligoID)
  ftab.info <- ftab.info[keep, ]
  rownames(ftab.info) <- ftab.info$featureID

  # LL 4.4.2012. With HITChip atlas we encountered some cases where the arrays had different number of entries
  # due to duplicates on some arrays. Now added automated handling here to avoid problems with cases
  # that may have different natural number of elements on the array.
  # JS 16.5.2013 NOT duplicates, see above. 
  # LL: 27.1.2014 fixed this part; in earlier version the probe matching was incorrect but this occurred only in few situations 

  # Form features x hybridizations matrix
  ftab <- matrix(NA, nrow = nrow(ftab.info), ncol = length(rawdata.esplit))
  rownames(ftab) <- rownames(ftab.info)
  colnames(ftab) <- names(rawdata.esplit)
  for (hid in names(rawdata.esplit)) { 
    inds <- match(ftab.info$featureID, rawdata.esplit[[hid]]$featureID)
    ftab[, hid] <- I(rawdata.esplit[[hid]][inds, "spatNormSignal"]) 
  }

  # Close MySQL connection
  dbDisconnect(con)

  # Clean up memory
  gc()

  list(data = ftab, info = ftab.info)
}
