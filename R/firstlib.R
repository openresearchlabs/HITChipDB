#' @import microbiome
#' @import parallel
#' @import tcltk
#' @importFrom RPA summarize_probedata

.onAttach <- function(lib, pkg)
{
  packageStartupMessage('\nHITChipDB R package (microbiome.github.com)\n(C) 2011-2016 Leo Lahti and Jarkko Salojarvi <microbiome-admin@googlegroups.com>\n')
}
