probe.summarization <- function (probedata, taxonomy, levels, summarization.methods) {

  # Ensure we have all levels
  levels <- intersect(levels, colnames(taxonomy))

  # Oligo summarization
  finaldata <- list()
  finaldata$oligo <- probedata

  for (level in levels) {

    finaldata[[level]] <- list()

    for (method in summarization.methods) {

        message(paste(level, method))
	# For species/L1/L2 summarization use the filtered phylogeny: phylogeny.filtered!
    	summarized.log10 <- summarize.probesets(
					taxonomy = taxonomy,		
			    		  oligo.data = log10(probedata), 
      			       	          method = method, 
					   level = level)$summarized.matrix

        # Store the data in absolute scale					
        finaldata[[level]][[method]] <- 10^summarized.log10

    }
  }

  finaldata

}



#' Description: summarize.probesets.through.species
#'
#' Arguments:
#'   @param level summarization level
#'   @param taxonomy oligo - phylotype matching data.frame
#'   @param oligo.data preprocessed probes x samples data matrix in log10 domain
#'   @param method summarization method
#'   @param verbose verbose
#'
#' Returns:
#'   @return List with two elements: summarized.matrix (summarized data matrix in log10 scale) and probe.parameters (only used with rpa, probe-level parameter estimates)
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

summarize.probesets.through.species <- function (level, taxonomy, oligo.data, method, verbose = FALSE) {

  # Get species-level summary first
  message("Summarizing through species...")
  probeset.summaries <- summarize.probesets.species(taxonomy, oligo.data, method, verbose = verbose)

  species.matrix   <- 10^probeset.summaries$summarized.matrix
  probe.parameters <- probeset.summaries$probe.parameters 

  if (level == "species") {
    res <- list(summarized.matrix = log10(species.matrix), probe.parameters = probe.parameters)
    return(res)
  }

  # List all species for the given level (L0 / L1 / L2)")
  phylogroups <- levelmap(phylotypes = NULL, from = level, to = "species", taxonomy)

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
    if (length(grep("rpa", method)) > 0) { vec <- colSums(mat) } # For RPA, use the sum for L1/L2

    summarized.matrix[pg, ] <- vec

  }

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


#' Description: summarize.probesets.directly
#'
#' Arguments:
#'   @param level summarization level
#'   @param phylogeny.info oligo - phylotype matching data.frame
#'   @param oligo.data preprocessed probes x samples data matrix in log10 domain
#'   @param method summarization method
#'   @param verbose verbose
#'
#' Returns:
#'   @return List with two elements: summarized.matrix (summarized data matrix in log10 scale) and probe.parameters (only used with rpa, probe-level parameter estimates)
#'
#' @export
#' @importFrom RPA d.update.fast rpa.fit
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

summarize.probesets.directly <- function (level, phylogeny.info, oligo.data, method, verbose = TRUE) {

  message("Direct summarization..")			     

  # Retrieve oligos for each taxonomic group
  probesets <- retrieve.probesets(phylogeny.info, level = level)

  # Number of phylotypes per oligo
  nPhylotypesPerOligo <- n.phylotypes.per.oligo(phylogeny.info, level) 

  # initialize
  summarized.matrix <- array(NA, dim = c(length(probesets),   ncol(oligo.data)), 
  		    	    dimnames = list(names(probesets), colnames(oligo.data))) 

  probe.parameters <- list()
  if (method == "frpa") {
    stop("fRPA not yet implement for direct summarization")
    if (verbose) {message("Loading pre-calculated preprocessing parameters")}
    rpa.hitchip.species.probe.parameters <- list()
    load(system.file("extdata/probe.parameters.rda", package = "HITChipDB"))
    probe.parameters <- rpa.hitchip.species.probe.parameters
  }

  for (set in names(probesets)) {
   
    # Pick data for particular probes
    probes <- probesets[[set]]

    # Pick probe data for the probeset: probes x samples
    # oligo.data assumed to be already in log10
    dat <- matrix(oligo.data[probes,], length(probes)) 
    rownames(dat) <- probes
    colnames(dat) <- colnames(oligo.data)

    if (method == "frpa") {

      message(method)

      # Summarize with pre-calculated variances
      vec <- d.update.fast(dat, probe.parameters[[set]])

    } else if (method == "rpa") {

      message(method)

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

   } else if (method == "rpa.with.affinities") {

      message(method)

     # Also include affinities in summarization

      	    res <- rpa.fit(dat, 
      	     		  alpha = 1 + 0.1*ncol(oligo.data)/2, 
			  beta  = 1 + 0.1*ncol(oligo.data)*nPhylotypesPerOligo[probes]^2, 
			  summarize.with.affinities = TRUE)

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

#' Description: Generate toydata for the package
#'
#' Arguments:
#'   @param output.dir output directory name
#'
#' Returns:
#'   @return output file name
#'
#' @export
#' @importFrom microbiome read.profiling
#' @importFrom microbiome GetPhylogeny
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

GenerateSimulatedData <- function (output.dir) {

  data.directory <- system.file("extdata/", package = "microbiome")

  #phylogeny.info <- read.profiling(level = "phylogeny.info", data.dir = data.directory)
  phylogeny.info <- as.data.frame(GetPhylogeny("HITChip", "filtered"))

  oligo.matrix.nolog.simulated <- read.profiling(data.dir = data.directory, log10 = FALSE)
  N <- ncol(oligo.matrix.nolog.simulated)
  colnames(oligo.matrix.nolog.simulated) <- paste("Sample.", 1:N, sep = "")

  # Oligo summarization
  finaldata <- list()
  finaldata[["oligo"]] <- oligo.matrix.nolog.simulated
  levels <- c("species", "L2", "L1")
  for (level in levels) {
    finaldata[[level]] <- list()
    for (method in c("sum", "rpa")) {

        message(paste(level, method))
    	summarized.log10 <- summarize.probesets(phylogeny.info, log10(oligo.matrix.nolog.simulated), method = method, level = level)$summarized.matrix
      			       	          
        # Store the data in absolute scale					
        finaldata[[level]][[method]] <- 10^summarized.log10

    }
  }

  # Write summary matrices into the output directory
  outd <- WriteChipData(finaldata, output.dir, phylogeny.info)

  set.seed(344)
  metadata.simulated <- data.frame(list(
              sampleID = I(colnames(oligo.matrix.nolog.simulated)),
	      time = rep(1:4, 5),
              age = runif(N, 0, 100),
              bmi = runif(N, 20, 40),
	      subjectID = I(paste("subjectID", rep(1:4, 5), sep = "-")),
	      group = I(sample(paste("group", rep(1:4, 5), sep = "-"))),
              gender = I(sample(c("M", "F"), N, replace = TRUE)),
              diet = I(sample(c("Apricots", "Beverages", "Carrots"), N, replace = TRUE))))
  
  write.table(metadata.simulated, file = paste(output.dir, "/metadata.tab", sep = ""), quote = FALSE, sep = "\t")

  output.dir

}



#' Description: count
#'
#' For cross-hyb control
#'
#' Arguments:
#'   @param d TBA
#'
#' Returns:
#'   @return list
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

count <- function(d){
   tabulate(d)
}



#' Description: summarize.oligos 
#'
#' Arguments:
#'   @param oligo.matrix oligo.matrix
#'   @param phylogeny.info phylogeny.info
#'   @param level taxonomic level
#'
#' Returns:
#'   @return list
#'
#' @export
#' @importFrom microbiome polish.phylogeny.info
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

summarize.oligos <- function(oligo.matrix, phylogeny.info, level = "L2"){

  phylogeny.info <- polish.phylogeny.info(phylogeny.info)
		 
  oligo <- split(phylogeny.info$oligoID,phylogeny.info[,level])
  oligo <- lapply(oligo, function(x) unique(as.character(x)))

  Sim <- sapply(oligo,function(x){  
    if (length(x) > 1)
      return(colSums(oligo.matrix[x, ]))
    else
      return(oligo.matrix[x, ])
  })
  t(Sim)
}


#' Description: mixingMatrix
#'
#' For cross-hyb control
#'
#' Arguments:
#'   @param phylogeny.info phylogeny.info
#'   @param level taxonomic level
#'
#' Returns:
#'   @return oligos x phylotypes mixing matrix
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

mixingMatrix <- function(phylogeny.info, level){

  M <- matrix(0,length(unique(phylogeny.info$oligoID)),length(unique(phylogeny.info[,level])),dimnames=list(sort(as.character(unique(phylogeny.info$oligoID))),sort(as.character(unique(phylogeny.info[,level])))))

  for (i in 1:nrow(phylogeny.info))
    M[as.character(phylogeny.info$oligoID[i]),as.character(phylogeny.info[i,level])]=1

  M <- apply(M, 2, function(x) x/sum(x))

  return(M)

}

#' Description: ngp: Non-negative probabilistic solution to probe mixing problem.
#'
#' Uses gamma prior to solve known cross-hybridisation problems. 
#' Perfect probe targets are given in phylogeny.info data frame. 
#'
#' Arguments:
#'   @param oligo.data oligo.data in absolute domain
#'   @param phylogeny.info oligo - phylotype mapping data frame
#'   @param level taxonomic level
#'   @param lambda - stregth of gamma prior. Default: 0.001.
#'   @param alpha - alpha parameter of gamma prior. Default: 1.
#'   @param beta - beta parameter of gamma prior. Default: 1.
#'
#' Returns:
#'   @return list
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

ngp <- function(oligo.data, phylogeny.info, level, lambda=0.001, alpha=1,beta=1){

   # Avoid NOTE in buildcheck
   ginv <- NULL

   # oligos x phylotypes mixing matrix for the given level
   M <- mixingMatrix(phylogeny.info, level)
   coms <- intersect(rownames(oligo.data), rownames(M))
   M <- M[coms, ]  
   oligo.data <- oligo.data[coms, ]

   # starting guess using pseudoinverse
   W <- t(M)%*%M
   X <- t(oligo.data) %*% M
   A <- ginv(W)%*%t(X)
   A[which(A<0)] <- 0

   for (i in 1:ncol(oligo.data)){
      Acol=A[,i]
      for (j in 1:ncol(M)){  
         a=2*W[j,j]
         b=sum((W[j,-j]+W[-j,j])*Acol[-j])-2*X[i,j]-beta*lambda
         c=lambda*(alpha-1)
         if (a!=0){
           D=b^2-4*a*c
           if (D>0)
             A[j,i]=max((-b+sqrt(D)),(-b-sqrt(D)))/(2*a)
         }else
            A[j,i]=-c/b
     }
   }
   colnames(A) <- colnames(oligo.data)
   rownames(A) <- colnames(M)
   return(A)
}

#' Description: Deconvolution
#'
#' !OBSOLETE! Used nmf package (currently removed from CRAN) for solving cross-hyb problem.
#' Calls function ngp with standard arguments (SEE: ngp).
#'
#' Arguments:
#'   @param oligo.data oligo.data in absolute domain
#'   @param phylogeny.info oligo - phylotype mapping data frame
#'   @param level taxonomic level
#'   @param block.solution block.solution
#'   @param verbose verbose
#'   @param ... parameters to be passed to ngp function
#'
#'
#' Returns:
#'   @return list
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Jarkko Salojarvi \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

deconvolution.nonneg <- function(oligo.data, phylogeny.info, level, block.solution = T,verbose=F,...){
   ngp(oligo.data, phylogeny.info, level,...)
}

# --------------------------------------------------------------------



