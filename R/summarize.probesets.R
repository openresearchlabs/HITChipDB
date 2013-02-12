#' Description: summarize.probesets
#'
#' Arguments:
#'   @param phylogeny.info oligo - phylotype matching data.frame
#'   @param oligo.data preprocessed probes x samples data matrix in log10 domain
#'   @param method summarization method
#'   @param level summarization level
#'   @param verbose print intermediate messages
#'   @param species.matrix Optional. Provide pre-calculated species-level summaries to speed up computation.
#'
#' Returns:
#'   @return List with two elements: summarized.matrix (summarized data matrix in log10 scale) and probe.parameters (only used with rpa, probe-level parameter estimates)
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
summarize.probesets <- function (phylogeny.info, oligo.data, method, level, verbose = TRUE, species.matrix = NULL) {

  # Option 1: Summarize probes through species level (default with RPA)
  if (method %in% c("rpa", "frpa", "sum.spe", "ave.spe")) {

    res <- summarize.probesets.through.species(level = level, phylogeny.info = phylogeny.info, oligo.data = oligo.data, method = gsub(".spe", "", method))

  # Option 2: Summarize from oligos to all levels directly 
  # (default with SUM and AVE)
  } else if (method %in% c("rpa.experimental", "frpa.experimental", "sum", "ave")) {

    res <- summarize.probesets.directly(level, phylogeny.info, oligo.data, gsub(".experimental", "", method))

  # NMF is always straight from oligos to levels
  } else if (method == "nmf") {
  
    if (level == "species") {
      warning("nmf oligo summarization method not implemented at species level");     return(NULL)
    }
  
    # Add +1 to avoid taking log10 for 0
    summarized.matrix <- 1 + deconvolution.nonneg(10^oligo.data, phylogeny.info, level)
    res <- list(summarized.matrix = summarized.matrix, probe.parameters = list())
  
  }

  res

}

