#' @title Define parameters in select box
#' @description Define parameters in select box.
#' @param con Output from dbConnect(dbDriver("MySQL"), username = dbuser, password = dbpwd, dbname = 'PhyloArray')
#' @param which.projects Optionally list which projects to take. All samples returned.
#' @param all.samples Return all samples from the selected projects: TRUE/FALSE
#' @param chip chip
#' @param save.dir save.dir
#' @param use.default.parameters use.default.parameters
#' @return list with defined parameters
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
ReadParameters <- function (con, which.projects = NULL, all.samples = TRUE, chip = NULL, save.dir = NULL, use.default.parameters = FALSE) {

  ## Determine the working directory
  if (is.null(save.dir)) {
    wdir <- tclvalue(tkchooseDirectory(title = "Save output files into directory:")) 
  } else {
    wdir <- save.dir
  }
    
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

  rs <- dbSendQuery(con, "SELECT phylogenyID, name FROM phylogeny WHERE NOT name='ROOT'")
  phylogenies <- fetch(rs, n = -1)
  phylogenies <- unique(phylogenies$name)
  defaults$phylogeny <- phylogenies[grep(defaults$phylogeny, phylogenies)][[1]] 

  if (is.null(use.default.parameters)) {
    use.default.parameters <- tk_select.list(c(paste("Yes, use the defaults:", s), "No, proceed to parameter selection"), preselect = paste("Yes, use the defaults:", s), multiple = FALSE, title = paste('Use default parameters?'))
  } 

  if (!use.default.parameters || substr(use.default.parameters, 1, 2) == "No") {    

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
