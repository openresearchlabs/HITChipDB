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
#'   @return project.info data.frame
#'
#' @examples # info <- fetch.sample.info(allowed.projects, dbuser, dbpwd, 
#'          # dbname, selected.samples = NULL, host = NULL, port = NULL)
#' @export
#' @importFrom DBI dbDriver
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
fetch.sample.info <- function (allowed.projects, chiptype = NULL, 
		  dbuser, dbpwd, dbname, 
		  selected.samples = NULL, 
		  host = NULL, port = NULL) { 

 # allowed.projects <- params$prj$projectName; chiptype = NULL; selected.samples = params$samples$sampleID
 # selected.samples = NULL

  drv <- dbDriver("MySQL")
  if (!(is.null(host) && is.null(port))) {
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname, host = host, port = port)
  } else { 
    con <- dbConnect(drv, username = dbuser, password = dbpwd, dbname = dbname)
  }  

  # Fetch all data from the database
  # Main info
   rs <- dbSendQuery(con, paste("SELECT p.projectName,p.projectID,s.subjectID,s.sampleID,s.samplingDate,s.normAlgVersion,h.hybridisationID,h.dye, h.arrayID ,s.reproducibility,s.normalisationFinished,s.imageID,fe.extractionID,fe.extractionName,fe.noSampleNormalisation,h.isDiscarded,fe.hasReproCheck
     FROM sample s               
     JOIN hybridisation h USING (sampleID) JOIN featureextraction fe USING (hybridisationID)
     JOIN project p USING (projectID)", 
     paste("WHERE projectName in ('", paste(unique(allowed.projects),collapse="','"), "')", sep = ""),
     "ORDER BY s.projectID, s.sampleID, h.hybridisationID, fe.extractionID", sep = " "))
     #paste("WHERE sampleID in ('", paste(unique(selected.samples),collapse="','"), "')", sep = ""),

  message("Fetch selected projects and samples")
  project.info.all <- fetch(rs, n = -1)

  # arrayID and barcode
  rs <- dbSendQuery(con, paste("SELECT a.arrayID,a.barcode,sl.designID 
     FROM array a               
     JOIN slide sl USING (barcode)
     WHERE arrayID in ('",paste(unique(project.info.all$arrayID),collapse="','"),"')",
     sep=""))
  project.info.arrays <- fetch(rs, n = -1)
  #combine 
  project.info.all <- cbind(project.info.all,project.info.arrays[match(project.info.all$arrayID,project.info.arrays$arrayID),c("barcode","designID")])
   
  # if no chiptype specified, use all
  if (is.null(chiptype)) {chiptype <- unique(project.info.all$designID)}
  if (is.null(selected.samples)) {selected.samples <- unique(project.info.all$sampleID)}

  # Pick selected samples only
  project.info.all <- project.info.all[project.info.all$sampleID %in% selected.samples,]

  # Close MySQL connection
  dbDisconnect(con) 

  # Filter out samples based on predefined criteria
  filter.table <- cbind(allowed.project = (project.info.all$projectName %in% allowed.projects), 
               	        notDiscarded = (!as.logical(project.info.all$isDiscarded)),
           		sampleNormalized = (!as.logical(project.info.all$noSampleNormalisation)),
           		normalisationFinished = (as.logical(project.info.all$normalisationFinished)),
			hasReproCheck = (project.info.all$hasReproCheck),
           		correctChip = (project.info.all$designID %in% chiptype),
           		correctNormAlgVersion = (project.info.all$normAlgVersion == 1.1),
    	   		selected.sample = (project.info.all$sampleID %in% selected.samples))

  filter.table[is.na(filter.table)] <- 0
  rkeep <- (rowMeans(filter.table == 1) == 1)

  # Remove annotations which are identical for all samples
  ckeep <- sapply(project.info.all, function (x) {!length(unique(x)) == 1})

  message("Filter the data")
  project.info <- project.info.all[rkeep, ckeep]

  project.info                  
     
}

