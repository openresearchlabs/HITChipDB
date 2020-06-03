#' @import microbiome
#' @import parallel
#' @importFrom RPA summarize_probedata
#' @importFrom graphics image
#' @importFrom graphics layout
#' @importFrom graphics lcm
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics plot.new
#' @importFrom graphics rect
#' @importFrom graphics strheight
#' @importFrom graphics strwidth
#' @importFrom graphics text
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom stats as.dist
#' @importFrom stats cor
#' @importFrom stats density
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom stats quantile
#' @importFrom stats sd
#' @importFrom utils head
#' @importFrom utils read.csv
#' @importFrom utils select.list
#' @importFrom utils sessionInfo
#' @importFrom utils str
#' @importFrom utils write.table
.onAttach <- function(lib, pkg)
{
  packageStartupMessage('\nHITChipDB R package (microbiome.github.com)\n(C) 2011-2016 Leo Lahti and Jarkko Salojarvi <microbiome-admin@googlegroups.com>\n')
}
