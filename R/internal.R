#' Description: Calculate species summaries and possibly update d.oligo2
#'
#' Arguments:
#' @param d.oligo2 d.oligo2
#' @param bgc.method background correction method
#' Returns:
#' @return Background-corrected data matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords internal
oligo.bg.correction <- function (d.oligo2, bgc.method) {

  if ( bgc.method == "6*sd bkg intensity" ){ bgth <- 6 }

  d.oligo2 <- threshold.data(d.oligo2, bgth)
  d.oligo2 <- apply(d.oligo2, c(1,2), function(x) max(0, x))
  
  d.oligo2

}

# Database utilities for package-internal use only

#' Tests whether the database connection is a phyloarray connection.
#' Expands one element (one "field", "value" pair list) from a list 
#' of "field", "value" pair lists
#'
#' @param con a MySQL database connection.
#' @return TRUE when the test succeeds. Otherwise a program halt.
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

phyloarrayConnection <- function (con) {

   # microbiome::InstallMarginal("RMySQL")

   if (!(class(con)=='MySQLConnection')) {
      stop('Input must be a DBI connection to a phyloarray database')
   }
   essential <- c("array",
                  "arraydesign",
                  "arrayfeature",
                  "arrayhybridisations",
                  "featureextraction",
                  "featuremeasurement",
                  "hybridisation",
                  "image",
                  "oligo",
                  "oligoclass",
                  "oligotargetpair",
                  "phylogeny",
                  "project",
                  "sample",
                  "slide",
                  "target",
                  "taxon",
                  "taxonlevel",
                  "taxtotax")
   if (!(length(intersect(dbListTables(con),essential))==length(essential))) {
      stop('Essential tables missing in the connected database. Not a phyloarray database?')
   }
   return(TRUE)
}



#' Expands one element (one "field", "value" pair list) from a list of "field", "value" pair lists
#' @param elm TBA
#'
#' @return TBA
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

expandElement <- function (elm) {
   if (is.list(elm)) {
      if (is.null(elm$field)) {
         stop("Missing 'field' field in condition")
      }
      if (!is.vector(elm$field)) {
         stop("Field 'field' must be vector with length 1")
      }
      if (length(elm$field)>1) {
         stop("Field 'field' must be vector with length 1")
      }
      if (!is.vector(elm$value)) {
         stop("'value' field must be vector")
      }
      if (is.null(elm$value)) {
         stop("Missing 'value' field in condition")
      }
      if (is.character(elm$value)) {
         elm$value <- paste("'",elm$value,"'",sep='')
      }
       return(paste("(",paste(elm$field,elm$value,sep='=',collapse=' OR '),")",sep=''))
   }
   else {
      stop("Argument 'condition' must be a list of lists with 'field' and 'value' pairs")
   }
}  

#' Expands a list of "field", "value" pair lists into an SQL condition
#' 
#' Arguments:
#'@param condition condition
#'
#' Returns:
#'  @return TBA
#'
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples # TBA
#' @keywords utilities
expandCondition <- function (condition) {
   if (is.null(condition)) {
      return(NULL)
   } else {
      if (is.list(condition)) {
         return(paste(" WHERE",paste(lapply(condition, expandElement),collapse=' AND ')))
      } else {
         stop("Argument 'condition' must be a list of lists with 'field' and 'value' pairs")
      }
   }  
}


#' Description: populate radiobuttons 
#'
#' Arguments:
#' @param tt TBA
#' @param title TBA
#' @param var.names TBA 
#' @param var.values TBA
#' @param var.init TBA
#'
#' Returns:
#'   @return A list.
#'
#' @references
#' See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
populate.radiobuttons <- function(tt, title, var.names, var.values, var.init) {

  title.font <- tkfont.create(weight = "bold", size = 10)

  frm <- tkframe(tt, relief = "groove", borderwidth = 2)

  label.widget <- tklabel(frm, text = title,font = title.font)

  tkpack(label.widget, side = "top")

  for (i in 1:length(var.values)){
    button.widget <- tkradiobutton(frm, text = var.names[i], 
                                  variable = var.init, value = var.values[i])
    tkpack(button.widget,side = "top")
  }

  tkpack(frm,side = "top")

  return(list(frame = frm, var = var.init))

}



#' Description: Select projects to analyze
#' 
#' Arguments:
#' @param con valid MySQL connection
#' @param multi enable selection of multiple options
#' @param condition TBA
#' Returns:
#'   @return vector of project names 
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
choose.projects <- function (con, multi = TRUE, condition = NULL) {
   prjs <- fetch.projects(con, condition = condition)
   projects <- select.list(sort(prjs$projectName), multiple = multi, title = "Select studies:")
   prjs <- fetch.projects(con, condition = list(list(field = 'projectName', value = projects)))
   return(prjs)
}


#' choose.samples
#'
#' @param con MySQL connection
#' @param multi multiple selections allowed
#' @param title title
#' @param condition TBA
#'
#' @return sample vector
#' @export 
#' @references See citation("microbiome")
#' @author Douwe Molenaar. Maintainer: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @examples # TBA
#' @keywords utilities
choose.samples <- function (con, multi = TRUE, title = 'Select samples:', condition = NULL) {

    smps <- fetch.samples(con, condition = condition)
    samples <- select.list(smps$sampleID, multiple = multi, title = title)
    smps <- fetch.samples(con, condition = list(list(field = 'sampleID', value = samples)))
    return(smps)

}



#' Description: Define parameters in select box
#'
#' Arguments:
#' @param con Output from dbConnect(dbDriver("MySQL"), username = dbuser, password = dbpwd, dbname = 'PhyloArray')
#' @param which.projects Optionally list which projects to take. All samples returned.
#' @param all.samples Return all samples from the selected projects: TRUE/FALSE
#' @param chip chip
#' 
#' Returns:
#'   @return list with defined parameters
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
ReadParameters <- function (con, which.projects = NULL, all.samples = TRUE, chip = NULL) {

  ## Determine the working directory
  wdir <- tclvalue(tkchooseDirectory(title = "Save output files into directory:")) 
        
  ## Choose samples to display
  if (is.null(which.projects)) {
    prj <- choose.projects(con, multi = TRUE, condition = NULL)
    if(nrow(prj) < 1) { stop("Choose at least 1 project") }
  } else {    
    prjs <- fetch.projects(con)
    prj <- prjs[match(which.projects, prjs$projectName), ]
  }

  defaults <- list(phylogeny = "16S", remove.nonspecific.oligos = FALSE, normalization = "minmax", all.samples = "Use all samples")
  s <- NULL; for (nam in names(defaults)) {s <- paste(s, paste(nam, ":", defaults[[nam]], sep = ""), "; ", sep = "")}

  if (chip == "MITChip") {
    defaults[["remove.nonspecific.oligos"]] <- TRUE
  }

  use.default.parameters <- tk_select.list(c(paste("Yes, use the defaults:", s), "No, proceed to parameter selection"), preselect = paste("Yes, use the defaults:", s), multiple = FALSE, title = paste('Use default parameters?'))

  rs <- dbSendQuery(con, "SELECT phylogenyID, name FROM phylogeny WHERE NOT name='ROOT'")
  phylogenies <- fetch(rs, n = -1)
  phylogenies <- unique(phylogenies$name)
  defaults$phylogeny <- phylogenies[grep(defaults$phylogeny, phylogenies)][[1]] 

  if (substr(use.default.parameters, 1, 2) == "No") {    

    ## Choose the phylogeny and the (lowest) summary taxonomic level to use
    if (length(phylogenies) > 1) {
      phylogeny <- tk_select.list(phylogenies, preselect = defaults$phylogeny, multiple = FALSE, title = 'Select phylogeny for profiling')
    } else {
      phylogeny <- phylogenies[[1]]
    }

    # Exclude non-specific oligos?
    remove.nonspecific.oligos <- tk_select.list(c("Yes", "No"), multiple = FALSE, preselect = defaults$remove.nonspecific.oligos, title = "Remove non-specific oligos?")
    if (remove.nonspecific.oligos == "Yes") {remove.nonspecific.oligos <- TRUE}
    if (remove.nonspecific.oligos == "No")  {remove.nonspecific.oligos <- FALSE}
   
    ## Normalization method
    scal <- tk_select.list(c("none", "minmax", "quantile"), preselect = defaults$normalization, multiple = FALSE, title = "Select normalization method")

    ## Sample selection
    all.samples <- tk_select.list(c("Select samples manually", "Use all samples"), preselect = defaults$all.samples, multiple = FALSE, title = "Select samples manually?")

  } else {

    phylogeny <- defaults$phylogeny
    remove.nonspecific.oligos <- defaults$remove.nonspecific.oligos
    scal <- defaults$normalization
    all.samples <- defaults$all.samples

  }

  # Extract samples
  if (all.samples == "Select samples manually") {
    # Select samples manually
    samples <- choose.samples(con, multi=TRUE, title='Select samples', condition=list(list(field='projectID', value=prj$projectID)))
  } else {
    # Take all samples in the project
    samples <- fetch.samples(con, condition = list(list(field='projectID', value=prj$projectID)))
  }

  # LL + JS decided to remove default BG correction 29.3.2012
  # based on empirical validations
  # bgc.method <- select.list(c("2*sd bkg intensity", "6*sd bkg intensity"), multiple = FALSE, preselect = "6*sd bkg intensity", title = "Select background correction method:")
  bgc.method <- NULL # Intentional

  # List oligos and phylotypes to remove by default
  rm.phylotypes <- phylotype.rm.list(chip) 

  list(wdir = wdir, prj = prj, samples = samples, phylogeny = phylogeny, normalization = scal, bgc.method = bgc.method, remove.nonspecific.oligos = remove.nonspecific.oligos, chip = chip, rm.phylotypes = rm.phylotypes)

}

#' Description: Calculate d.oligo2
#'
#' Arguments:
#' @param featuretab featuretab
#' @param d.scaled d.scaled
#' @param oligo.ids oligo.ids
#' Returns:
#'   @return d.oligo2
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
get.doligo2 <- function (featuretab, d.scaled, oligo.ids) {

  d.oSplit <- split(cbind(featuretab[,1:3],d.scaled), featuretab$oligoID)
  d.oSplit.pruned <- d.oSplit[oligo.ids]
  d.oligo <- t(sapply(d.oSplit, function(x) apply((x[,4:dim(x)[2]]), 2, mean, na.rm=T))) # hybs separate
  sampleID <- get.sampleid(d.oligo)

  d.oligo2 <- t(sapply(d.oSplit,
                     function(x){
                       temp <- apply((x[,4:dim(x)[2]]), 2, mean, na.rm=T)
                       temp2 <- sapply(split(temp, sampleID), mean, na.rm=T)
                       return(temp2)
                     }
                     ))# hybs averaged
  
  d.oligo2
}




#' Description: Pick sampleIDs from d.oligo column names
#'
#' Arguments:
#' @param d.oligo d.oligo matrix
#' Returns:
#'   @return sampleID vector
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
get.sampleid <- function (d.oligo) {

  sampleID <- sapply(colnames(d.oligo), function(z) { # edit by S.Tims
    # allows samples with "." in the samplename
    s <- strsplit(z, split="\\.")[[1]] 
    if(length(s) > 4){
       s <- head(s,-3)
       s <- paste(s, collapse = ".")
    } else { s <- s[1] }
   return(s)})

   sampleID
}



#' Description: Minmax scaling. 
#'
#' Arguments:
#' @param dat data matrix in original absolute scale
#' @param quantile.points quantiles for minmax
#' @param minmax.points pre-calculated quantiles for minmax in log10 scale; overrides quantile.points argument
#' @param robust Select minmax version. 
#' 
#' Returns:
#'   @return normalized data matrix in absolute scale
#'
#' @note With robust = FALSE, the standard minmax is carried out. This
#'   shifts and scales each array such that their min and max values are
#'   identical across arrays. The robust = TRUE will perform the scaling
#'   such that the upper quantiles and minimum values of the data match 
#'   (instead of maximum values).
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
scaling.minmax <- function (dat, quantile.points = NULL, minmax.points = NULL, robust = FALSE) {

  if ( is.null(minmax.points) ) {

    maxabs <- mean(apply(dat, 2, quantile, max(quantile.points), na.rm = TRUE))
    minabs <- mean(apply(dat, 2, quantile, min(quantile.points), na.rm = TRUE))

  } else {

    maxabs <- max(minmax.points)
    minabs <- min(minmax.points)

  }

  if (!robust) {

    r <- apply(dat, 2, function (x) { 
      x = (((x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T)))*(maxabs - minabs)) + minabs;
    return(x)})

  } else {
    
    r <- apply(dat, 2, function (x) {

    # Shift data to start from zero
    xz <- x - min(x, na.rm = T);

    # Check the quantile points
    maxq <- quantile(xz, max(quantile.points), na.rm = TRUE);   

    # Determine the scaling factor such that the max quantiles will 
    # match between arrays
    k <- maxabs/maxq;

    # Scale the data to match max quantiles
    xs <- k * xz + min(dat, na.rm = TRUE);
    xs})
  }

  r
}





#' Description: Write log file
#'
#' Arguments:
#' @param naHybs hybridisation that were removed due to NAs  
#' @param params parameters
#'
#' Returns:
#'   @return List of scaling methods
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

WriteLog <- function (naHybs, params) {

  scriptVersion1 <- sessionInfo()$otherPkgs$microbiome$Version # microbiome package number
  scriptVersion2 <- sessionInfo()$otherPkgs$HITChipDB$Version # microbiome package number

  ## Write log of parameters used in profiling in to the file
  tmpTime <- strsplit(as.character(Sys.time()), split=" ")[[1]]
  tmpDate <- tmpTime[1]
  tmpTime <- paste(strsplit(tmpTime[2], split=":")[[1]], collapse=".")
  profTime <- paste(tmpDate,tmpTime,sep="_")
  logfilename <- paste(params$wdir,"/",profTime,"_profiling_log.txt", sep="")

  cat("Log of profiling script\n", "\n", file=logfilename)
  cat("profiling date: ",profTime, "\n", file=logfilename, append=T)
  cat("script version microbiome: ", scriptVersion1,  "\n",file=logfilename, append=T)
  cat("script version HITChipDB: ", scriptVersion2,  "\n",file=logfilename, append=T)
  cat("data retrieved from db: ",params$useDB,  "\n", file=logfilename, append=T)
  cat("project IDs: ",params$prj$projectID,  "\n", file=logfilename, append=T)
  cat("sample IDs: ",params$samples$sampleID,  "\n", file=logfilename, append=T)
  cat("excluded oligos: ",params$rm.phylotypes$oligos,  "\n", file=logfilename, append=T)
  cat("excluded species: ",params$rm.phylotypes[["species"]], "\n", file=logfilename, append=T)
  cat("excluded level 1: ",params$rm.phylotypes[["L1"]], "\n", file=logfilename, append=T)
  cat("excluded level 2: ",params$rm.phylotypes[["L2"]], "\n", file=logfilename, append=T)
  cat("excluded hybridisations: ",naHybs,  "\n", file=logfilename, append=T)
  cat("remove non-specific oligos: ",params$remove.nonspecific.oligos, "\n",file=logfilename, append=T)
  cat("phylogeny: ",params$phylogeny,  "\n", file=logfilename, append=T)
  cat("scaling: ",params$scal,  "\n", file=logfilename, append=T)
  cat("data in directory: ",params$wdir, "\n",file=logfilename, append=T)

  ## Save profiling parameters 
  paramfilename <- paste(params$wdir,"/",profTime,"_profiling_params.Rdata", sep="")
  save(logfilename, profTime, scriptVersion1, scriptVersion2, params, naHybs, file = paramfilename)  

  list(log.file = logfilename, parameter.file = paramfilename)
  
}


#' Description: Writed data into the output directory
#' @param finaldata preprocessed data matrices in absolute scale (from the chipdata function)
#' @param output.dir output directory
#' @param tax.table tax.table used in summarization
#' @param tax.table.full tax.table.full unfiltered phylogenyinfo
#' @param meta sample metadata samples x features
#' @param verbose verbose
#' @return Preprocessed data in absolute scale, tax.table, and parameters
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
WriteChipData <- function (finaldata, output.dir, tax.table, tax.table.full, meta, verbose = TRUE) {

  ## Write oligoprofile in original (non-log) domain
  fname <- paste(output.dir, "/oligoprofile.tab", sep = "")
  mydat <- finaldata[["oligo"]]
  WriteMatrix(cbind(rownames(mydat), mydat), fname, verbose)
    
  ## Write the other levels in log domain
  for (level in setdiff(names(finaldata), "oligo")) {
    for (method in names(finaldata[[level]])) {
        fname <- paste(output.dir, "/", level, "-", method, ".tab", sep = "")
        mydat <- finaldata[[level]][[method]]
        WriteMatrix(cbind(rownames(mydat), mydat), fname, verbose)
    }
  }
  
  # Write metadata template
  fname <- paste(output.dir, "/meta.tab", sep = "")
  WriteMatrix(meta, fname, verbose)  

  # Write tax.table that shall be used for probe summarization
  fname <- paste(output.dir, "/taxonomy.tab", sep = "")
  WriteMatrix(tax.table, fname, verbose)

  # Write filtered tax.table 
  fname <- paste(output.dir, "/taxonomy.filtered.tab", sep = "")
  WriteMatrix(tax.table, fname, verbose)

  # Write unfiltered tax.table
  fname <- paste(output.dir, "/taxonomy.full.tab", sep = "")
  WriteMatrix(tax.table.full, fname, verbose)

   # Write metadata
  fname <- paste(output.dir, "/meta.tab", sep = "")
  WriteMatrix(meta, fname, verbose)

  # Return path to the output directory 
  output.dir
 
}




#' Description: Format string vector to mysql query format
#' 
#' Arguments:
#' @param s string vector
#'
#' Returns:
#'   @return mysql query version
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
mysql.format <- function (s) {
  paste("('",paste(as.character(s),collapse="','",sep="'"),"')",sep="")
}






#' Description: Calculate species summaries and possibly update d.oligo2
#'
#' Arguments:
#' @param d.oligo2 d.oligo2
#' @param bgc.method background correction method
#' Returns:
#'   @return Background-corrected data matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

oligo.bg.correction <- function (d.oligo2, bgc.method) {

  if ( bgc.method == "6*sd bkg intensity" ){ bgth <- 6 }

  d.oligo2 <- threshold.data(d.oligo2, bgth)
  d.oligo2 <- apply(d.oligo2, c(1,2), function(x) max(0, x))
  
  d.oligo2

}

#' Description: Between-arrays normalization 
#'
#' Arguments:
#' @param dat data matrix in original absolute scale
#' @param method normalization method
#' @param bg.adjust background adjustment 
#' @param minmax.quantiles quantiles for minmax
#' @param minmax.points minmax end points
#'
#'
#' Returns:
#'   @return Normalized data matrix in absolute scale
#'
#' @export
#' @importFrom preprocessCore normalize.quantiles
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

ScaleProfile <- function (dat, method = 'minmax', bg.adjust = NULL, minmax.quantiles = c(0.005, 0.995), minmax.points = NULL) {

  # d.scaled <- ScaleProfile(fdat.orig, params$normalization, bg.adjust = NULL, minmax.points = params$minmax.points) 
  # dat <- fdat.orig; method = params$normalization; bg.adjust = NULL; minmax.quantiles = c(0.005, 0.995); minmax.points = NULL

  message(paste("Normalizing with", method))
  
  ## Table dat is a copy of featuretab containing 
  ## logarithms of the values in 
  ## featuretab

  if (method=='minmax') {
    r <- scaling.minmax(dat, quantile.points = minmax.quantiles, minmax.points = minmax.points, robust = FALSE)
  } else if (method=='minmax.robust') {
    r <- scaling.minmax(dat, quantile.points = minmax.quantiles, minmax.points = minmax.points, robust = TRUE)
  } else if (method=='quantile') {
    dn <- dimnames(dat)
    r <- normalize.quantiles(dat)
    dimnames(r) <- dn
  } else if (method=='normExpQuant') {
    ## Impute NA's with sample medians
    na.inds <- which(is.na(dat), arr.ind=T)
    r <- apply(dat, 2, function(x){x[is.na(x)] <- median(x, na.rm=T); return(x)})
    dn <- dimnames(dat)
    rc <- apply(10^(r), 2, bg.adjust)
    r <- normalize.quantiles(log10(rc+1))
    dimnames(r) <- dn
    r[na.inds] <- NA
  } else {
    stop("No between-array normalization recognized!!")
  }
 
  return(r)
}


#' Description: determine detection threshold for the data
#'
#' Arguments:
#' @param dat data
#' @param sd.times standard deviation threshold
#'
#' Returns:
#'   @return thresholded data matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

threshold.data <- function(dat, sd.times = 6){

  thr <- apply(dat, 2, function(x){
      DD <- density(as.numeric(x), adjust = 1.2, na.rm = T);
      noise_mode <- DD$x[which(DD$y==max(DD$y))[1]];
      noise_sd   <- sd(x[x < noise_mode], na.rm = T);
      low.thresh <- noise_mode + sd.times*noise_sd;
      low.thresh 
    })

  # Subtract background from signal intensities in each sample
  data.mat <- t(apply(dat, 1, function(Tr){ Tr-thr })) 
  return(data.mat)
}


#' Description: List oligos associated with removed phylotypes from different levels
#'
#' Arguments:
#' @param rm.phylotypes rm.phylotypes
#' @param tax.table tax.table
#'
#' Returns:
#'   @return probe name vector
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

sync.rm.phylotypes <- function (rm.phylotypes, tax.table) {

  # If remove L0 is not NULL, then 
  # add L1 groups under this group to removal list
  if (!is.null(rm.phylotypes$L0)) {
    rm.phylotypes$L1 <- c(rm.phylotypes$L1, tax.table[which(tax.table[["L0"]] == rm.phylotypes$L0), "L1"])
    rm.phylotypes$L1 <- unique(rm.phylotypes$L1)
  }

  # If remove L1 is not NULL, then add L2 groups under this group to removal list
  if (!is.null(rm.phylotypes$L1)) {
    rm.phylotypes$L2 <- c(rm.phylotypes$L2, tax.table[which(tax.table[["L1"]] == rm.phylotypes$L0), "L2"])
    rm.phylotypes$L2 <- unique(rm.phylotypes$L2)
  }

  # If remove L2 is not NULL, then 
  # add species groups under this group to removal list
  if (!is.null(rm.phylotypes$L2)) {
    rm.phylotypes$species <- c(rm.phylotypes$species, tax.table[which(tax.table[["L2"]] == rm.phylotypes$L0), "species"])
    rm.phylotypes$species <- unique(rm.phylotypes$species)
  }

  rm.phylotypes
}





