#' read.profiling.010
#' 
#' Description: read profiling script output from profiling script v. 010 into R 
#'
#' Arguments:
#'   @param level phylogenetic level ("oligo" / "species" / "L1" / "L2" / "L0") or "phylogeny.info"
#'   @param method ("rpa" / "sum" / "ave" / "nmf")
#'   @param data.dir Profiling script output directory for reading the data. If not given, GUI will ask to specify the file and overruns the possible level / method arguments in the function call.
#'   @param log10 Logical. Logarithmize the data TRUE/FALSE. By default, the data is in original non-log scale.
#'   @param impute impute missing oligo signals
#'
#' Returns:
#'   @return data matrix (phylo x samples)
#'
#' @export
#' @importFrom svDialogs tk_choose.files
#' @importFrom svDialogs tclServiceMode
#'
#' @examples # params <- read.profiling.010(...); 
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

read.profiling.010 <- function(level = NULL, method = "rpa", data.dir = NULL, log10 = TRUE, impute = TRUE){

  if (level == "L0") { level <- "level0"}
  if (level == "L1") { level <- "level1"}
  if (level == "L2") { level <- "level2"}
  if (method == "rpa") {method <- "RPA"}
  if (method == "sum") {method <- "Sum"}
  if (method == "ave") {method <- "log10Ave"}
  if (method == "nmf") {method <- "NMF"}

  InstallMarginal("svDialogs")

  ##  Select file
  if (is.null(data.dir)) {

    f <- tk_choose.files(multi = F)

    # Recognize level and method from the file name 
    level <- NULL; method <- NULL
    if (!length(grep("oligo", f)) == 0) { level <- "oligo"}
    if (!length(grep("species", f)) == 0) { level <- "species"}
    if (!length(grep("level0", f)) == 0) { level <- "level0"}
    if (!length(grep("level1", f)) == 0) { level <- "level1"}
    if (!length(grep("level2", f)) == 0) { level <- "level2"}
    if (!length(grep("phylogeny", f)) == 0) { level <- "phylogeny.info"}
    if (!length(grep("RPA", f)) == 0) { method <- "RPA"}
    if (!length(grep("Sum", f)) == 0) { method <- "Sum"}
    if (!length(grep("Ave", f)) == 0) { method <- "log10Ave"}
    if (!length(grep("NMF", f)) == 0) { method <- "NMF"}

    tclServiceMode(FALSE)

  } else {
    if (level %in% c("level0", "level1", "level2", "species")) {
      f <- paste(data.dir, "/", level, "_", method, "_010.tab", sep = "")
    } else if (level == "oligo") {
      f <- paste(data.dir, "/oligoprofile_010.tab", sep = "")
    } else if (level == "phylogeny.info") {
      f <- paste(data.dir, "/phylogenyinfo_010.tab", sep = "")
    }
  }

  # level2_Sum_010.tab  
  # oligoprofile_010.tab  
  # species_log10Ave_010.tab

  # Read the data
  if (level == "phylogeny.info") {
    tab <- read.table(f, sep = "\t", header = T, row.names = NULL)
  } else {
    tab <- read.table(f, sep = "\t", header = T, row.names = 1)
  }

  # Check that the data is logarithmized as required in the arguments
  if (level == "phylogeny.info") {
    tab <- tab
  } else if (!length(grep("log10", f)) == 0 && !log10) { 
    tab <- 10^tab
  } else if (length(grep("log10", f)) == 0 && log10) {
    tab <- log10(tab)
  } else {
    tab <- tab
  }

  # Always impute by rows and for log10 data
  if (impute && any(is.na(tab))) {
    warning(paste("The matrix has ", sum(is.na(tab)), " missing values - imputing."))
    if (!log10) {
      tab <- 10^t(impute(t(log10(tab))))
    } else {
      tab <- t(impute(t(tab)))
    }
  }

  as.matrix(tab)

}



#' Phylogeneticenrichments
#'
#' Description: Computes the dependence of the higher phylogenetic
#' level groups of the given list x, which can be a vector/data frame of either
#' oligos, species, or level2 groups. Computes only the enrichments
#' if onlyEnrich=T (neglects the unexpected disappearences of the
#' groups in the given list). Outputs only the enrichments below the
#' given p-value threshold p.th.
#' Requires phylogeny.info=the first 4 columns from phylogenyprofile,
#' origlevel=the level from which the list is, and maplevel= the
#' level to which the mapping should be made and the enrichments of
#' which computed.
#' Outputs a four component list: p-values from Fisher's tests,
#' input tables for tests, the actual tests outputs, and the full
#' phylogenetic mapping information for the given list.
#' Changes: Version 1 computed two-tailed p-values also for single
#' occurrences on origlevel.
#' v2: two-tailed p-values, but not for single occurrences
#' v3: one-tailed p-values (enrichments), if wanted, outputting only
#' enrichments under the given p-value
#'
#' Arguments:
#'   @param x x
#'   @param phylogeny.info phylogeny.info
#'   @param origlevel origlevel
#'   @param maplevel maplevel 
#'   @param onlyEnriched onlyEnriched
#'   @param p.th p-value threshold
#'
#' Returns:
#'   @return enrichments
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

Phylogeneticenrichments <- function(x, phylogeny.info, origlevel = colnames(phylogeny.info)[3], maplevel = "L1", onlyEnriched = T, p.th = 0.05)
{
  
  ## Convert to character vector
  x <- as.character(x)
    
  ## collect origlevel groups for all x items
  origlevel.split <- split(phylogeny.info,phylogeny.info[[origlevel]])
  x.origlevelGroups <- origlevel.split[x]
  
  ## vector (inX and inXind) showing which items from the origlevel of
  ## phylogeny.info are in the x list
  origlevel.ugroups <- names(origlevel.split)
  inXind <- match(x, origlevel.ugroups)
  inX <- rep(F,length(origlevel.ugroups))
  inX[inXind] <- T

  ## Collect the full phyloinfo
  phyloM <-  x.origlevelGroups[[1]]
  nolGroups <- length(x.origlevelGroups)
  if(nolGroups>1)
    for(i in 2:nolGroups){
      phyloM <- rbind(phyloM,x.origlevelGroups[[i]])
    }
  
  if(length(x)>2){
    
    ## compute enrichments of maplevel groups in x
    e <- c()
    pvals <- c()
    estimates <- c()
    tables <- list()
    maplevel.ugroups <- as.character(unique(phyloM[[maplevel]]))

    for(g in maplevel.ugroups){

      ## compute in which maplevel groups the given item occurs
      inMaplevelGroup <- unlist(lapply(origlevel.split,
                                         function(y){is.element(g,y[[maplevel]])
                                                   }))

      if(sum(inMaplevelGroup)>0){
        if(onlyEnriched)
          tmp <- try(fisher.test(inX, inMaplevelGroup, alternative="g"))
        else
          tmp <- try(fisher.test(inX, inMaplevelGroup, alternative="t"))

          e[g] <- list(tmp)
          pvals[g] <- tmp$p.value          
          tables[g] <- list(table(inX,inMaplevelGroup))
          estimates[g] <- tmp$estimate

      }

      pvals.adjusted <- p.adjust(pvals, method = "BH")
      inds <- (pvals.adjusted < p.th)
      e <- e[inds]
      pvals <- pvals.adjusted[inds]
      tables <- tables[inds]
      estimates <- estimates[inds]

    }

    ## Return enrichments i.e. Fisher's exact tests in table form

    if (length(pvals)>0)
       list(pval.adj=t(rbind(pvals,estimates)), 
           tables=tables, tests=e, phylogeny.info=phyloM) 
    else
       list(pvalues=pvals, tables=tables, tests=e, phylogeny.info=phyloM) 

  } else
    list(pvalues=1, tables=NULL, tests=NULL, phylogeny.info=phyloM) 
}

#' Description: Check number of matching phylotypes for each probe
#' 
#' Arguments:
#'   @param phylogeny.info oligo - phylotype matching data.frame
#'   @param level phylotype level
#'
#' Returns:
#'   @return number of matching phylotypes for each probe
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

n.phylotypes.per.oligo <- function (phylogeny.info, level) {
  sapply(split(phylogeny.info[, c("oligoID", level)], phylogeny.info$oligoID), function(x) length(unique(x[[level]])))
}

#' Description: Write matrix in tab file
#'
#' Arguments:
#'   @param dat data matrix
#'   @param filename output file
#'   @param verbose verbose
#' Returns:
#'   @return output file location
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

WriteMatrix <- function (dat, filename, verbose = FALSE) { 

  if (verbose) { message(paste("Writing output in ", filename)) }
  write.table(dat, file = filename, quote = FALSE, sep = "\t", row.names = FALSE)
  filename

}



#' Description: determine threshold for bg correction
#'
#' Arguments:
#'   @param dat data matrix (in approximately normal scale ie. logged)
#'
#' Returns:
#'   @return threshold value
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

estimate.min.threshold <- function (dat) {

  #estimate min threshold 
  DD <- density(as.numeric(unlist(dat)))

  #find mode
  noise_mode <- DD$x[[which.max(DD$y)]] # DD$x[which(diff(DD$y)<0)[1]]

  #compute sd of noise
  noise_sd <- sd(dat[dat < noise_mode])

  #threshold
  low.thresh <- noise_mode + 6*noise_sd

  low.thresh
}


