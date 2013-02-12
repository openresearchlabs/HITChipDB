
#' Description: filter 16S data
#' 
#' Arguments:
#'   @param full16S full16S
#'   @param pmTm.margin pmTm margin
#'   @param complement logical
#'   @param mismatch logical
#'   @param rmoligos oligos to exclude
#'   @param remove.nonspecific.oligos Logical. Remove oligos with multiple targets.
#'
#' Returns:
#'   @return filtered 16S data
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

prune16S <- function (full16S, pmTm.margin = 2.5, complement = 1, mismatch = 0, rmoligos = NULL, remove.nonspecific.oligos = FALSE) {

<<<<<<< HEAD
  message("Applying probe filters..")
  message(paste("Hybridisation temperature: Tm >= pmTm - ", pmTm.margin, "\n", sep = ""))
  message(paste("Number of mismatches: ", mismatch, "\n", sep = ""))
  message(paste("Must be a complement sequence (complement = ", complement, ")\n", sep = ""))
  message("No requirement for a full-length hybridisation\n\n")
=======
    message("Applying probe filters..")
    message(paste("Hybridisation temperature: Tm >= pmTm - ", pmTm.margin, "\n", sep = ""))
    message(paste("Number of mismatches: ", mismatch, "\n", sep = ""))
    message(paste("Must be a complement sequence (complement = ", complement, ")\n", sep = ""))
    message("No requirement for a full-length hybridisation\n\n")
>>>>>>> 388dc8fb9ee265a39825c8979f8f011f1e6e6cc6

  keep <- ((full16S$Tm >= (full16S$pmTm-pmTm.margin)) & (full16S$complement == complement)) & (full16S$mismatch == mismatch)

  phylogeny.info <- full16S[keep, ]

  rmoligos2 <- rmoligos
  if (remove.nonspecific.oligos) {
    # if (verbose) {message("Removing oligos that have multiple targets at L2 level")}
    nPhylotypesPerOligo <- n.phylotypes.per.oligo(phylogeny.info, "L2") 
    nonspecific.oligos <- setdiff(phylogeny.info$oligoID, names(which(nPhylotypesPerOligo == 1)))
    rmoligos2 <- c(rmoligos, nonspecific.oligos)
  } 

  phylogeny.info <- phylogeny.info[!phylogeny.info$oligoID %in% rmoligos2, ]

  phylogeny.info
}