#' Description: summarize.probesets.through.species
#'
#' Arguments:
#'   @param level summarization level
#'   @param phylogeny.info oligo - phylotype matching data.frame
#'   @param oligo.data preprocessed probes x samples data matrix in log10 domain
#'   @param method summarization method
#'
#' Returns:
#'   @return List with two elements: summarized.matrix (summarized data matrix in log10 scale) and probe.parameters (only used with rpa, probe-level parameter estimates)
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

summarize.probesets.through.species <- function (level, phylogeny.info, oligo.data, method) {

  # level = level; phylogeny.info = phylogeny.info; oligo.data = oligo.data; method = gsub(".spe", "", method); probe.parameters = probe.parameters

  # Get species-level summary first
  probeset.summaries <- summarize.probesets.species(phylogeny.info, oligo.data, method, verbose = FALSE)

  species.matrix   <- 10^probeset.summaries$summarized.matrix
  probe.parameters <- probeset.summaries$probe.parameters 

  if (level == "species") {
     res <- list(summarized.matrix = log10(species.matrix), probe.parameters = probe.parameters)
     return(res)
  }

  # List all species for the given level (L0 / L1 / L2)")
  phylogroups <- levelmap(phylotypes = NULL, level.from = level, level.to = "species", phylogeny.info)

  # Remove specified phylogroups
  # phylogroups <- phylogroups[setdiff(names(phylogroups), rm.phylotypes[[level]])]

  summarized.matrix <- matrix(NA, nrow = length(phylogroups), ncol = ncol(oligo.data))
  rownames(summarized.matrix) <- sort(names(phylogroups))
  colnames(summarized.matrix) <- colnames(oligo.data)

  # Go through each phylogroup and summarize from species level
  for (pg in names(phylogroups)) {

    specs <- unique(phylogroups[[pg]])
    mat <- matrix(species.matrix[specs,], nrow = length(specs))

    if (method == "ave") { vec <- colMeans(mat) }
    if (method == "sum") { vec <- colSums(mat)  } 
    if (method %in% c("rpa", "frpa")) { vec <- colSums(mat) } # For RPA, use the sum for L1/L2

    summarized.matrix[pg, ] <- vec

  }

  list(summarized.matrix = log10(summarized.matrix), probe.parameters = probe.parameters)

}