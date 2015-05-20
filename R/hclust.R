

#' calculate.hclust
#' 
#' Calculate hierarchical clustering for standard selections in 
#' profiling script
#'
#' 
#'   @param dat data matrix (use log10 with pearson!)
#'   @param method hierarchical clustering method (see ?hclust)
#'   @param metric clustering metrics (spearman / pearson / euclidean)
#'
#' 
#'   @return hclust object for log10 and for absolute scale data
#'
#' @export
#' @examples 
#' \dontrun{  
#'   data(peerj32)
#'   dat <- peerj32$microbes
#'   hc <- calculate.hclust(dat, 'complete', 'pearson')
#' }
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
calculate.hclust <- function(dat, method = "complete", metric = "pearson") {
    
    if (metric == "euclidean") {
        hc <- hclust(dist(t(dat)), method = method)
    } else if (metric %in% c("spearman", "pearson")) {
        hc <- hclust(as.dist(1 - cor(dat, use = "complete.obs", 
                     method = metric)), 
            method = method)
    } else {
        stop("Provide proper metric for calculate.hclust!")
    }
    
    hc
    
}
