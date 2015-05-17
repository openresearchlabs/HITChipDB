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
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

GenerateSimulatedData <- function (output.dir) {

  data.directory <- system.file("extdata/", package = "microbiome")

  #phylogeny.info <- read.profiling(level = "phylogeny.info", data.dir = data.directory)
  phylogeny.info <- as.data.frame(GetPhylogeny("HITChip", "filtered"))

  oligo.matrix.nolog.simulated <- read.profiling(level = "oligo", data.dir = data.directory, log10 = FALSE)
  N <- ncol(oligo.matrix.nolog.simulated)
  colnames(oligo.matrix.nolog.simulated) <- paste("Sample.", 1:N, sep = "")

  # Oligo summarization
  finaldata <- list()
  finaldata[["oligo"]] <- oligo.matrix.nolog.simulated
  levels <- c("species", "L2", "L1")
  for (level in levels) {
    finaldata[[level]] <- list()
    for (method in c("sum", "rpa", "nmf")) {

        message(paste(level, method))
    	summarized.log10 <- summarize.probesets(phylogeny.info, log10(oligo.matrix.nolog.simulated), method = method, level = level)$summarized.matrix
      			       	          
        # Store the data in absolute scale					
        finaldata[[level]][[method]] <- 10^summarized.log10

    }
  }

  # Write summary matrices into the output directory
  outd <- HITChipDB::WriteChipData(finaldata, output.dir, phylogeny.info)

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


