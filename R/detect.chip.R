#' Description: Detect the chip type (H/M/PITChip) from database name
#'
#' Arguments:
#'   @param dbname MySQL database name
#' 
#' Returns:
#'   @return chip name
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

detect.chip <- function (dbname) {

  if (dbname %in% c("PhyloArray_MIT", "phyloarray_mit", "Phyloarray_MIT")) {
    chip <- "MITChip"
  } else if (dbname %in% c("PhyloArray_PIT", "phyloarray_pit", "Phyloarray_PIT", "pitchipdb")) {
    chip <- "PITChip"
  } else if (dbname %in% c("PhyloArray_HIT", "phyloarray_hit", "Phyloarray_HIT", "Phyloarray", "phyloarray")) {
    chip <- "HITChip"
  } else if (dbname %in% c("chickchipdb")) {
    chip <- "ChickChip"
  } else {
    stop("Check database name (dbname)!")
  }

  chip

}
