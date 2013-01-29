#' summarize.probesets
#'
#' Arguments:
#'   @param phylogeny.info oligo - phylotype matching data.frame
#'   @param oligo.data preprocessed probes x samples data matrix in log10 domain
#'   @param method summarization method
#'   @param level summarization level
#'   @param verbose print intermediate messages
#'   @param rm.phylotypes Phylotypes to exclude (a list with fields species, L1, L2)
#'   @param species.matrix Optional. Provide pre-calculated species-level summaries to speed up computation.
#' Returns:
#'   @return List with two elements: summarized.matrix (summarized data matrix in log10 scale) and probe.parameters (only used with rpa, probe-level parameter estimates)
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

summarize.probesets <- function (phylogeny.info, oligo.data, method, level, verbose = TRUE, rm.phylotypes = NULL, species.matrix = NULL) {

  # oligo.data <- oligo.log10; method <- method.name; verbose = TRUE; species.matrix = NULL

  rm.phylotypes <- sync.rm.phylotypes(rm.phylotypes, phylogeny.info)

  # print("Remove specified oligos")
  rm.oligos <- rm.phylotypes$oligos
  if (!is.null(rm.oligos)) { oligo.data <- oligo.data[setdiff(rownames(oligo.data), rm.oligos), ]}
  phylogeny.info <- phylogeny.info[!phylogeny.info$oligoID %in% rm.oligos, ]

  if (verbose) message("Get species matrix in original scale") # if not pre-calculated version available
  if (is.null(species.matrix)) {
    probeset.summaries <- summarize.probesets.species(phylogeny.info, oligo.data, method, verbose = FALSE, rm.phylotypes$species)
    species.matrix   <- 10^probeset.summaries$summarized.matrix
    probe.parameters <- probeset.summaries$probe.parameters # optinally with rpa
  } else {
    probe.parameters <- list()
  }

  if (level == "species") {

    summarized.matrix <- species.matrix

  } else if (level %in% c("L0", "L1", "L2")) {

    if (method %in% c("rpa", "ave", "sum", "rpa.full")) {

      if (verbose) {message(paste(level, method))}

      # List all species for the given level (L0 / L1 / L2)")
      phylogroups <- levelmap(phylotypes = NULL, level.from = level, level.to = "species", phylogeny.info)

      # Remove specified phylogroups
      phylogroups <- phylogroups[setdiff(names(phylogroups), rm.phylotypes[[level]])]

      summarized.matrix <- matrix(NA, nrow = length(phylogroups), ncol = ncol(oligo.data))
      rownames(summarized.matrix) <- sort(names(phylogroups))
      colnames(summarized.matrix) <- colnames(oligo.data)

      for (pg in names(phylogroups)) {
        specs <- unique(phylogroups[[pg]])
        mat <- matrix(species.matrix[specs,], nrow = length(specs))

        if (method == "ave") { vec <- colMeans(mat) }
        if (method == "sum") { vec <- colSums(mat)  } 
        if (method %in% c("rpa", "rpa.full")) { vec <- colSums(mat)  } # For RPA, use the sum for L1/L2

        summarized.matrix[pg, ] <- vec
      }

    } else if (method == "nmf") {

      # Add +1 to avoid taking log10 for 0
      summarized.matrix <- 1 + deconvolution.nonneg(10^oligo.data, phylogeny.info, level)

    }

  } else {

    message(level)
    message(nchar(level))
    message(colnames(phylogeny.info))
    stop("Provide proper level!")

  }

  # Return in the original log10 domain    
  list(summarized.matrix = log10(summarized.matrix), probe.parameters = probe.parameters)

}

#' Description: Probeset summarization with various methods.
#' 
#' Arguments:
#'   @param phylogeny.info oligo - phylotype matching data.frame
#'   @param oligo.data preprocessed probes x samples data matrix in log10 domain
#'   @param method summarization method
#'   @param verbose print intermediate messages
#'   @param rm.species Species to exclude
#' Returns:
#'   @return List with two elements: summarized.matrix (summarized data matrix in log10 scale) and probe.parameters (only used with rpa, probe-level parameter estimates)
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

summarize.probesets.species <- function (phylogeny.info, oligo.data, method, verbose = TRUE, rm.species = c("Victivallis vadensis")) {

  level <- "species"			    

  if (method == "nmf") {warning("nmf oligo summarization method not implemented at species level"); return(NULL)}
  probesets <- retrieve.probesets(phylogeny.info, level = level)
  probesets <- probesets[setdiff(names(probesets), rm.species)]
  nPhylotypesPerOligo <- n.phylotypes.per.oligo(phylogeny.info, level) 

  # initialize
  summarized.matrix <- array(NA, dim = c(length(probesets), ncol(oligo.data)), 
  		    	      dimnames = list(names(probesets), colnames(oligo.data))) 

  probe.parameters <- list()

  for (set in names(probesets)) {

    if (verbose) { message(set) }

    # Pick expression for particular probes
    probes <- probesets[[set]]

    # Pick probe data for the probeset: probes x samples
    # oligo.data assumed to be already in log10
    dat <- matrix(oligo.data[probes,], length(probes)) 
    rownames(dat) <- probes
    colnames(dat) <- colnames(oligo.data)

    if (method == "rpa") {

      # RPA is calculated in log domain
      # Downweigh non-specific probes with priors with 10% of virtual data and
      # variances set according to number of matching probes
      # This will provide slight emphasis to downweigh potentially
      # cross-hybridizing probes
      vec <- rpa.fit(dat, tau2.method = "robust", 
      	     		  alpha = 1 + 0.1*ncol(oligo.data)/2, 
			  beta  = 1 + 0.1*ncol(oligo.data)*nPhylotypesPerOligo[probes]^2)$mu

      # Switch to this when the probe variances estimated from large atlas become available.		  
      # Remember to compare performance
      # vec <- d.update.fast(dat, variances)
      # vec <- d.update.fast(dat - affinities, variances)}, mc.cores = mc.cores), identity))

    } else if (method == "rpa.full") {

      # RPA is calculated in log domain
      # Downweigh non-specific probes with priors with 10% of virtual data and
      # variances set according to number of matching probes
      # This will provide slight emphasis to downweigh potentially
      # cross-hybridizing probes
      res <- rpa.fit(dat, tau2.method = "robust", 
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