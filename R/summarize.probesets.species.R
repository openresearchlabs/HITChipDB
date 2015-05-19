#' Description: Probeset summarization with various methods.
#' 
#' Arguments:
#'   @param phylogeny.info oligo - phylotype matching data.frame
#'   @param oligo.data preprocessed probes x samples data matrix in log10 domain
#'   @param method summarization method
#'   @param verbose print intermediate messages
#'   @param rm.species Species to exclude
#'   @param probe.parameters Optional. If probe.parameters are given,
#'          the summarization is based on these and model parameters are not
#' 	    estimated. A list. One element for each probeset with the following probe vectors: 
#'	    affinities, variances
#' Returns:
#'   @return List with two elements: summarized.matrix (summarized data matrix in log10 scale) and probe.parameters (only used with rpa, probe-level parameter estimates)
#'
#' @export
#' @importFrom RPA d.update.fast 
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

summarize.probesets.species <- function (phylogeny.info, oligo.data, method, verbose = TRUE, rm.species = c("Victivallis vadensis"), probe.parameters = list()) {

  level <- "species"			    
  probesets <- retrieve.probesets(phylogeny.info, level = level)
  probesets <- probesets[setdiff(names(probesets), rm.species)]
  nPhylotypesPerOligo <- n.phylotypes.per.oligo(phylogeny.info, level) 

  # initialize
  summarized.matrix <- array(NA, dim = c(length(probesets), ncol(oligo.data)), 
  		    	      dimnames = list(names(probesets), colnames(oligo.data))) 

  probe.parameters <- list()
  if (method == "frpa") {
    if (verbose) {message("Loading pre-calculated preprocessing parameters")}
    rpa.hitchip.species.probe.parameters <- list()
    load(system.file("extdata/probe.parameters.rda", package = "HITChipDB"))
    probe.parameters <- rpa.hitchip.species.probe.parameters

    # Ensure we use only those parameters that are in the filtered phylogeny
    for (bac in names(probe.parameters)) {
      probe.parameters[[bac]] <- probe.parameters[[bac]][intersect(names(probe.parameters[[bac]]), probesets[[bac]])]
    }
  }

  for (set in names(probesets)) {

    # if (verbose) { message(set) }

    # Pick expression for particular probes
    probes <- probesets[[set]]

    # Pick probe data for the probeset: probes x samples
    # oligo.data assumed to be already in log10
    dat <- matrix(oligo.data[probes,], length(probes)) 
    rownames(dat) <- probes
    colnames(dat) <- colnames(oligo.data)

    if (method == "frpa") {

      # Summarize with pre-calculated variances
      vec <- d.update.fast(dat, probe.parameters[[set]])

    } else if (method == "rpa") {

      # RPA is calculated in log domain
      # Downweigh non-specific probes with priors with 10% of virtual data and
      # variances set according to number of matching probes
      # This will provide slight emphasis to downweigh potentially
      # cross-hybridizing probes
      res <- rpa.fit(dat, 
      	     		  alpha = 1 + 0.1*ncol(oligo.data)/2, 
			  beta  = 1 + 0.1*ncol(oligo.data)*nPhylotypesPerOligo[probes]^2)

      vec <- res$mu
      probe.parameters[[set]] <- res$tau2

    } else if (method == "ave") {

      vec <- log10(colMeans((10^dat), na.rm = T))

    } else if (method == "sum") {

      # Weight each probe by the inverse of the number of matching phylotypes
      # Then calculate sum -> less specific probes are downweighted
      # However, set the minimum signal to 0 in log10 scale (1 in original scale)!
      dat2 <- (10^dat) / nPhylotypesPerOligo[rownames(dat)]
      dat2[dat2 < 1] <- 1
      vec <- log10(colSums(dat2, na.rm = T))
      vec[which(vec == -Inf)] <- 0

    }
    
    summarized.matrix[set, ] <- vec 

  }

  list(summarized.matrix = summarized.matrix, probe.parameters = probe.parameters)
  
}
