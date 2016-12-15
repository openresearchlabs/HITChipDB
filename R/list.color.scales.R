#' @title List color scales
#' @description List color scales
#' @param ... Arguments to be passed
#' @return list of color scales
#' @export
#' @examples list.color.scales()
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
list.color.scales <- function(...) {
    ## Different colour scales
    list(`white/blue` = colorRampPalette(c("white", "darkblue"), 
         interpolate = "linear")(100), 
         `white/black` = colorRampPalette(c("white", "black"), 
         interpolate = "linear")(100), 
        `black/yellow/white` = colorRampPalette(c("black", "yellow", "white"), 
         bias = 0.5, interpolate = "linear")(100))
}

