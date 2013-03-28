#' Description: Compress normalized raw data matrix into final probe-level matrix:
#'              summarize oligos into probes and hybridisations into samples
#'
#' Arguments:
#'   @param fdat.log10 normalized raw data matrix oligos x hybridisations in log10 scale
#'   @param fdat.hybinfo hybridization info table
#'   @param fdat.oligoinfo oligo info table
#'   @param oligo.ids oligo.ids
#' 
#' Returns:
#'   @return probes x samples matrix in log10 scale
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

summarize.rawdata <- function (fdat.log10, fdat.hybinfo, fdat.oligoinfo, oligo.ids) {

  # List rows for each oligo (each oligo has multiple features which will be averaged)
  d.oSplit <- split(1:nrow(fdat.log10), fdat.oligoinfo$oligoID)[as.character(oligo.ids)] 

  # Remove oligos with no data (probes discarded earlier from the data)
  d.oSplit <- d.oSplit[!sapply(d.oSplit, is.null)]

  # Probes x hybs: oligo summary as means of log feature signals per oligo, hybs separate
  message("Probe summarization with mean of log feature signals per oligo, hybs separate")
  oligo.data  <- t(sapply(d.oSplit, function(x) colMeans(fdat.log10[x,], na.rm = TRUE)))

  ## Average over all hybridisations/extractions associated with this sample
  # List hybridisations associated with the same sample
  indlist <- split(1:ncol(oligo.data), fdat.hybinfo$sampleID)

  # Keep only cases with multiple (typically 2) hybs per sample 
  # (see table(sapply(indlist, length)))
  message("Removing samples with only one hybridisation")
  indlist <- indlist[sapply(indlist, length) > 1]

  message("Forming probes x samples matrix (average over hybridisations for each sample)")
  oligo.data <- sapply(indlist, function(inds) { rowMeans(oligo.data[, inds]) } )

  oligo.data
}