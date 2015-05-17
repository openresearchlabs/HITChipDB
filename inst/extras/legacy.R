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



