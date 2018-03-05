#' @title Fetch data from the database
#' @description Fetch data from the database.
#' @param params params 
#' @param con con
#' @param scriptVersion scriptVersion
#' @param save.data save.data
#' @param scaling scaling
#' @param cmetrics cmetrics
#' @return data 
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
FetchData <- function (params, con, scriptVersion, save.data, scaling, cmetrics) {

  ## COLLECTING DATA FROM THE DATABASE
  message("Collecting data from the database\n")

  ## Collecting data for ALL probes (used to be: oligo's) JN09032011
  inclsamples <- paste("sampleID='", params$samples$sampleID, "'", sep = "", collapse = " OR ")

  rs <- dbSendQuery(con,"DROP TABLE IF EXISTS tmp1")
  query <- paste("CREATE TEMPORARY TABLE tmp1 ",
                 "SELECT s.sampleID, s.projectID, h.hybridisationID, fe.extractionID, h.dye, p.oligoID, p.probeID, af.featureID, ",
                 " fm.isOutlier, spatNormSignal AS featureSignal ",
                 "FROM sample s USE INDEX (PRIMARY) ",
                 "JOIN hybridisation h USING (sampleID) ",
                 "JOIN featureextraction fe USING (hybridisationID) ",
                 "JOIN featuremeasurement fm USING (extractionID) ",
                 "JOIN arrayfeature af USE INDEX (PRIMARY) USING (featureID) ",
                 "JOIN probe p USING (probeID) ",
                 "WHERE (",inclsamples,") ",
                 "AND normalisationFinished ",
                 "AND NOT (oligoID IS NULL) ",
                 "AND NOT isDiscarded ",
                 "AND NOT noSampleNormalisation ",
                 "ORDER BY s.sampleID, h.hybridisationID,  p.probeID, af.featureID",
                 sep ="")

  rs <- dbSendQuery(con, query)
  rs <- dbSendQuery(con, "ALTER TABLE tmp1 ADD INDEX (featureID)")
  rs <- dbSendQuery(con, 'SELECT * FROM tmp1')

  rawdata <- fetch(rs, n = -1)

  ## Check if there is any data
  if(nrow(rawdata) == 0)
    stop("No data found for these samples (perhaps they are not normalized yet?).\n\n")

  ## Create the data matrix (featuretab) for clustering based on all array features, 
  ## each hybridisation having one column in this table and each feature having one row. 
  ## The first column contains the oligoID's
  message("Create the FULL data matrix (featuretab) for clustering.\n")

  ## Change outlier values to NAs
  rawdata$featureSignal[as.logical(rawdata$isOutlier)] <- NA

  ## Get the hybIDs and initialize featuretab with the first hyb
  hybIDs <- unique(rawdata$hybridisationID)
  featuretab <- rawdata[rawdata$hybridisationID==hybIDs[1], c("oligoID","probeID","featureID", "featureSignal")]

  ## Name the column with sampleID, hybID, and dye
  samplename <- unique(rawdata[rawdata$hybridisationID==hybIDs[1],"sampleID"])
  dye <- unique(rawdata[rawdata$hybridisationID==hybIDs[1],"dye"])
  projectname <- unique(rawdata[rawdata$hybridisationID==hybIDs[1],"projectID"])
  columnnames <- paste(samplename, hybIDs[1], dye, projectname, sep=".")

  message(str(featuretab))
  if (length(hybIDs)>1) {
    for (i in hybIDs[2:length(hybIDs)]) {
      addcol <- rawdata[rawdata$hybridisationID==i,"featureSignal"]
      featuretab <- cbind(featuretab, addcol)

      ## Name the columns with sampleID, hybID, and dye
      samplename <- unique(rawdata[rawdata$hybridisationID==i,"sampleID"])
      dye <- unique(rawdata[rawdata$hybridisationID==i,"dye"])
      projectname <- unique(rawdata[rawdata$hybridisationID==i,"projectID"])
      columnnames <- c(columnnames, paste(samplename, i, dye, projectname, sep="."))
    }
  }

  colnames(featuretab) = c("oligoID","probeID","featureID", columnnames)
  rownames(featuretab) <- featuretab$featureID 

  ## Discard the hybs contains only NAs
  onlyNA <- colSums(is.na(featuretab))==dim(featuretab)[1]
  naHybs <- names(onlyNA)[onlyNA]
  if(sum(onlyNA)>0){
    cat("Removing the following hybs, because they contain only NAs:\n")
    cat(naHybs,"\n\n")
    featuretab <- featuretab[,!onlyNA]
  }

  # Remove rmoligos
  featuretab <- featuretab[!featuretab$oligoID %in% params$rm.phylotypes$oligos, ]

  ## Write log of parameters used in profiling in to the file
  tmpTime <- strsplit(as.character(Sys.time()), split=" ")[[1]]
  tmpDate <- tmpTime[1]
  tmpTime <- paste(strsplit(tmpTime[2], split=":")[[1]], collapse=".")
  profTime <- paste(tmpDate,tmpTime,sep="_")
  logfilename <- paste(params$wdir,"/",profTime,"_profiling_log.txt", sep="")

  cat("Log of profiling script\n", "\n", file=logfilename)
  cat("profiling date: ",profTime, "\n", file=logfilename, append=T)
  cat("script version: ",scriptVersion,  "\n",file=logfilename, append=T)
  cat("data retrieved from db: ",params$useDB,  "\n", file=logfilename, append=T)
  cat("project IDs: ",params$prj$projectID,  "\n", file=logfilename, append=T)
  cat("sample IDs: ",params$samples$sampleID,  "\n", file=logfilename, append=T)
  cat("excluded oligos: ",params$rmoligos,  "\n", file=logfilename, append=T)
  cat("excluded hybridisations: ",naHybs,  "\n", file=logfilename, append=T)
  cat("phylogeny: ",params$phylogeny,  "\n", file=logfilename, append=T)
  cat("scaling: ",params$scal,  "\n", file=logfilename, append=T)
  cat("clustering tree in: ",params$clusterGraphFile,  "\n", file=logfilename, append=T)
  cat("tree ratio: ",params$figureratio, "\n",file=logfilename, append=T)
  cat("clustering metric: ",params$clmet, "\n",file=logfilename, append=T)
  cat("phylogeny level in figure: ",params$lev, "\n",file=logfilename, append=T)
  cat("figure coloring: ", params$pal, "\n",file=logfilename, append=T)
  cat("figure fontsize: ", params$fontsize, "\n",file=logfilename, append=T)
  cat("data saved: ", save.data, "\n",file=logfilename, append=T)
  cat("data in directory: ",params$wdir, "\n",file=logfilename, append=T)

  ## Write parameters used in profiling to the file
  paramfilename <- paste(params$wdir,"/",profTime,"_profiling_params.Rdata", sep="")

  save(logfilename, profTime, scriptVersion, params, naHybs, cmetrics, save.data, file=paramfilename)
  
  ## Collect the full phylogenetic information for oligos  
  message("Collect the full 16S phylogeny\n")

  full16Squery <- paste("SELECT l1.name AS 'level 1', l2.name AS 'level 2', ", 
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
                        "WHERE ph.name='",params$phylogeny,"' ",
                        "AND l1.taxonLevel='level 1' ",
                        "AND tt1.nodeDistance=1 ",
                        "AND tt2.nodeDistance=1 ",
                        "ORDER BY l1.name, l2.name, species.name, specimen.name, o.oligoID;", sep="")
  rs <- dbSendQuery(con, full16Squery)
  full16S <- fetch(rs, n = -1)
  full16S <- full16S[full16S$oligoID %in% params$rm.phylotypes$oligos,]          
	  
  message("FINISHED COLLECTING THE DATA\n")

  list(full16S = full16S, featuretab = featuretab, logfilename = logfilename)

}
