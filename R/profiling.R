#' run.profiling.script
#' 
#' Description: Profiling main script
#'
#' Arguments:
#'   @param dbuser MySQL username
#'   @param dbpwd  MySQL password
#'   @param dbname MySQL database name (HITChip: "Phyloarray"; MITChip: "Phyloarray_MIT";
#'                                PITChip old: "Phyloarray_PIT"; PITChip new: "pitchipdb")
#'   @param verbose verbose
#'   @param host host; needed with FTP connections
#'   @param port port; needed with FTP connections
#'   @param summarization.methods List summarization methods to be included in output. For HITChip frpa always used; for other chips, rpa always used. Other options: sum, ave
#'   @param which.projects Optionally specify the projects to extract. All samples from these projects will be included.
#'
#' Returns:
#'   @return Profiling parameters. Also writes output to the user-specified directory.
#'
#' @importFrom microbiome summarize_probedata
#' @importFrom microbiome read_hitchip
#' @importFrom phyloseq otu_table
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

run.profiling.script <- function (dbuser, dbpwd, dbname, verbose = TRUE, host = NULL, port = NULL, summarization.methods = c("frpa", "sum"), which.projects = NULL) {

  htree.plot <- NULL		     

  # Fetch and preprocess the data		     
  chipdata  <- preprocess.chipdata(dbuser, dbpwd, dbname, 
  	       				verbose = verbose,
					   host = host,
					   port = port, 
		 	  summarization.methods = summarization.methods, 
			         which.projects = which.projects)

  probedata <- chipdata$probedata
  params    <- chipdata$params

  # Phylogeny used for L1/L2/species summarization
  taxonomy <- chipdata$taxonomy

  # Complete phylogeny before melting temperature etc. filters
  taxonomy.full <- chipdata$taxonomy.full

  # Create sample metadata template
  meta <- data.frame(list(index = 1:ncol(probedata), 
       	                  sample = colnames(probedata)), 
			  stringsAsFactors = FALSE)

  # Write preprocessed probe-level data, taxonomy and 
  # metadata template in tab-delimited file
  outd <- WriteChipData(list(oligo = probedata), params$wdir, 
       	  	taxonomy, taxonomy.full, meta, verbose = verbose)

  # Summarize probes into species abundance table
  abundance.tables <- list()
  abundance.tables$oligo <- probedata

  for (method in summarization.methods) {

    output.dir <- params$wdir
    pseq <- read_hitchip(output.dir, method = method, 
    	    			     detection.threshold = 0)$pseq
    spec <- otu_table(pseq)@.Data
    abundance.tables[["species"]][[method]] <- spec

    for (level in setdiff(colnames(taxonomy), c("species", "specimen", "oligoID", "pmTm"))) {
      
      # TODO Switch to this when DB access is running again and possible to 
      # check consistency.
      taxo <- unique(taxonomy[, c(level, "species")])
      rownames(taxo) <- as.character(taxo$species)

      # This includes pseudocount +1 in each cell
      pseq <- hitchip2physeq(t(spec), meta, taxo, detection.limit = 0)
      # This not; compatible with earlier
      # pseq <- hitchip2physeq(t(spec) - 1, meta, taxo, detection.limit = 0)

      tg <- tax_glom(pseq, level)
      ab <- tg@otu_table
      rownames(ab) <- as.character(as.data.frame(tax_table(tg))[[level]])
      ab <- ab[order(rownames(ab)),]
      abundance.tables[[level]][[method]] <- ab

      #ab2 <- species2higher(spec, taxonomy, level, method)
      #abundance.tables[[level]][[method]] <- ab2

    }

  }   

  ## Write preprocessed data in tab delimited file
  outd <- WriteChipData(abundance.tables, params$wdir, taxonomy, taxonomy.full, meta, verbose = verbose)
  
  # Add oligo heatmap into output directory
  # Provide oligodata in the _original (non-log) domain_
  hc.params <- add.heatmap(log10(probedata), 
  	          output.dir = params$wdir, taxonomy = taxonomy)

  # Plot hierachical clustering trees into the output directory
  dat <- abundance.tables[["oligo"]]

  if (ncol(dat) > 2) { 

    if (params$chip == "MITChip") {
      # With MITChip, use the filtered phylogeny for hierarchical clustering
      dat <- dat[unique(taxonomy$oligoID),]
    }

    # Clustering
    # Save into file
    method <- "complete"
    hc <- hclust(as.dist(1 - cor(log10(dat), use = "pairwise.complete.obs", method = "pearson")), method = method)
    pdf(paste(params$wdir, "/hclust_oligo_pearson_", method, "_", nrow(dat), "probes", ".pdf", sep = ""), height = 800, width = 800 * ncol(dat)/20)
    plot(hc, hang = -1, main = "hclust/pearson/oligo/log10/complete", xlab = "Samples", ylab = "1 - Correlation")
    dev.off()

  }

  # Plot hclust trees on screen
  tmp <- htree.plot(dat)

  # Write parameters into log file
  tmp <- WriteLog(chipdata$naHybs, params)
  params$logfilename <- tmp$log.file
  params$paramfilename <- tmp$parameter.file

  params

}

species2higher <- function (species.matrix, taxonomy, level, method) {

  # List all species for the given level (L0 / L1 / L2)")
  phylogroups <- levelmap(NULL, from = level, to = "species", taxonomy)

  summarized.matrix <- matrix(NA, nrow = length(phylogroups), ncol = ncol(species.matrix))
  rownames(summarized.matrix) <- sort(names(phylogroups))
  colnames(summarized.matrix) <- colnames(species.matrix)

  # Go through each phylogroup and summarize from species level
  for (pg in names(phylogroups)) {

    specs <- unique(phylogroups[[pg]])
    mat <- matrix(species.matrix[specs,], nrow = length(specs))

    if (method == "ave") { vec <- colMeans(mat) }
    if (method == "sum") { vec <- colSums(mat)  } 
    if (length(grep("rpa", method)) > 0) { vec <- colSums(mat) } # For RPA, use the sum for L1/L2

    summarized.matrix[pg, ] <- vec

  }

  summarized.matrix

}


#' add.heatmap
#' Description: Add oligprofile heatmap into output directory
#'
#' Arguments:
#'   @param dat oligoprofile data in original (non-log) domain
#'   @param output.dir output data directory
#'   @param output.file output file name
#'   @param taxonomy oligo-phylotype mappings
#'   @param ppcm figure size
#'   @param hclust.method hierarchical clustering method
#'   @param palette color palette ("white/black" / "white/blue" / "black/yellow/white")
#'   @param level taxonomic level to show
#'   @param metric clustering metric
#'   @param figureratio figure ratio
#'   @param fontsize font size
#'   @param tree.display tree.display
#'
#' Returns:
#'   @return Plotting parameters
#'
#' @export
#' @examples # data(peerj32); hc <- add.heatmap(peerj32$microbes[, 1:4])
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

add.heatmap <- function (dat, output.dir, output.file = NULL, taxonomy, ppcm = 150, 
	         hclust.method = "complete", palette = "white/black", level = "L1", metric = "pearson", 
  		 figureratio = 10, fontsize = 40, tree.display = TRUE) {

  # dat <- finaldata[["oligo"]]; output.dir = params$wdir;  output.file = NULL; taxonomy = taxonomy; ppcm = 150; hclust.method = "complete"; palette = "white/blue"; level = "L2"; metric = "pearson"; figureratio = 12; fontsize = 12; tree.display = TRUE
  #output.dir = "~/tmp/";  output.file = NULL; taxonomy = taxonomy; ppcm = 150; hclust.method = "complete"; palette = "white/blue"; level = "L2"; metric = "pearson"; figureratio = 12; fontsize = 12; tree.display = TRUE

  if (is.null(output.file)) {
    output.file <- paste(output.dir,"/", gsub(" ", "", level), "-oligoprofileClustering.pdf",sep="")
  }		 

  hc.params <- list()
  if( ncol(dat) >= 3 ) {

    message(paste("Storing oligo heatmap in", output.file))  
    hc.params$ppcm <- ppcm
    hc.params$output.file <- output.file

    # PLOT THE HEATMAP
    # figure width as a function of the number of the samples
    plotdev <- pdf(output.file, 
  	    width = max(trunc(ppcm*21), trunc(ppcm*21*ncol(dat)/70)), 
	    height = trunc(ppcm*29.7)) 
    try(hc.params <- PlotPhylochipHeatmap(data = dat,
                taxonomy = taxonomy,
                metric = metric,
                level = level,
                tree.display = tree.display,
                palette = palette,
                fontsize = fontsize,
                figureratio = figureratio, 
		hclust.method = hclust.method)) 

    dev.off()
  }

  hc.params

}







#' Description: Default list of removed phylotypes and oligos
#'
#' Arguments:
#'  @param chip Chip name (HIT/MIT/PIT/Chick)Chip
#' Returns:
#'   @return List of removed oligos and phylotypes
#'
#' @export
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

phylotype.rm.list <- function (chip) {

  rm.phylotypes <- list()

  if (chip == "HITChip") {
    
    rm.phylotypes[["oligos"]] <- c("UNI 515", "HIT 5658", "HIT 1503", "HIT 1505", "HIT 1506")
    rm.phylotypes[["species"]] <- c("Victivallis vadensis")
    rm.phylotypes[["L1"]] <- c("Lentisphaerae")
    rm.phylotypes[["L2"]] <- c("Victivallis")

  } else if (chip == "MITChip") {

    rm.phylotypes[["oligos"]] <- c("Bacteria", "DHC_1", "DHC_2", "DHC_3", "DHC_4", "DHC_5", "DHC_6", "Univ_1492")
    rm.phylotypes[["species"]] <- c()
    rm.phylotypes[["L1"]] <- c()
    rm.phylotypes[["L2"]] <- c()

  } else if (chip == "PITChip") {

    # Based on JZ mail 9/2012; LL

    rm.old.oligos <- c("Bacteria", "DHC_1", "DHC_2", "DHC_3", "DHC_4", "DHC_5", "DHC_6", "Univ_1492")
    rm.new.oligos <- c("PIT_1083", "PIT_1022", "PIT_1057", "PIT_1023", "PIT_1118", "PIT_1040", "PIT_1058", "PIT_1119", "PIT_122", "PIT_1221", "PIT_1322", "PIT_1367", "PIT_1489", "PIT_160", "PIT_1628", "PIT_1829", "PIT_1855", "PIT_1963", "PIT_1976", "PIT_1988", "PIT_2002", "PIT_2027", "PIT_2034", "PIT_2101", "PIT_2196", "PIT_2209", "PIT_2281", "PIT_2391", "PIT_2392", "PIT_2418", "PIT_2425", "PIT_2426", "PIT_2498", "PIT_2555", "PIT_2563", "PIT_2651", "PIT_2654", "PIT_2699", "PIT_2741", "PIT_2777", "PIT_2786", "PIT_2936", "PIT_35", "PIT_425", "PIT_427", "PIT_428", "PIT_429", "PIT_435", "PIT_481", "PIT_605", "PIT_7", "PIT_733", "PIT_734", "PIT_892")
    rm.phylotypes[["oligos"]] <- c(rm.old.oligos, rm.new.oligos)
    rm.phylotypes[["species"]] <- c()
    rm.phylotypes[["L0"]] <- c("Nematoda", "Apicomplexa", "Euryarchaeota", "Ascomycota", "Parabasalidea", "Chordata")
    rm.phylotypes[["L1"]] <- c("Chromadorea", "Coccidia", "Methanobacteria", "Saccharomycetales", "Trichomonada", "Mammalia")
    rm.phylotypes[["L2"]] <- c("Ascaris suum et rel.", "Eimeria  et rel.", "Methanobrevibacter et rel.", "Saccharomyces et rel.", "Trichomonas et rel.", "Uncultured Mammalia", "Uncultured methanobacteria")

  } else if (chip == "ChickChip") {
    warning("No universal probes excluded from ChichChip yet!")
  }

  rm.phylotypes

}

