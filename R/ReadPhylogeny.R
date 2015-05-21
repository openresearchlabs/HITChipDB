ReadPhylogeny <- function ( 
	      	 phylogeny, rm.phylotypes, remove.nonspecific.oligos,
	    		     dbuser = dbuser, 
			     dbpwd = dbpwd, 
			     dbname = dbname, 
			     host = host, 
			     port = port,
			     verbose = verbose, 
			     chip) {

    message("Fetching Phylogeny from the database")
    phylogeny.full <- get.phylogeny.info(phylogeny, 
	    		     dbuser = dbuser, 
			     dbpwd = dbpwd, 
			     dbname = dbname, 
			     host = host, 
			     port = port,
			     verbose = verbose, 
			     chip = chip)

    ph <- polish.hitchip.phylogeny(phylogeny.full, chip, rm.phylotypes, 
       	  				remove.nonspecific.oligos)
    phylogeny.full <- ph$full
    phylogeny.filtered <- ph$filtered

    
  list(full = phylogeny.full, filtered = phylogeny.filtered)

}


polish.hitchip.phylogeny <- function (phylogeny.full, chip, rm.phylotypes, remove.nonspecific.oligos) {

    # Fix an issue with PITChip2:	
    #> table(phylogeny.full[grep("Ignatzschineria", phylogeny.full$L2),"L2"])
    #Ignatzschineria et al. Ignatzschineria et rel. 
    #                 64                      61 
    phylogeny.full[grep("Ignatzschineria et al.", phylogeny.full$L2),"L2"] <- "Ignatzschineria et rel."
    phylogeny.full[grep("^Clostridiales$", phylogeny.full$L2),"L2"] <- "Clostridium \\(sensu stricto\\)" 
    phylogeny.full[grep("^Clostridia$", phylogeny.full$L2),"L2"] <- "Clostridium \\(sensu stricto\\)" 

    # This handles also pmTm, complement and mismatch filtering
    # This is the phylogeny used in probe summarization into taxonomic levels
    rm.oligos <- sync.rm.phylotypes(rm.phylotypes, phylogeny.full)$oligos

    phylogeny.filtered <- prune16S(phylogeny.full, pmTm.margin = 2.5, complement = 1, mismatch = 0, rmoligos = rm.phylotypes$oligos, remove.nonspecific.oligos = remove.nonspecific.oligos)

    # Remove certain species from summarization phylogeny.filtered
    rm.species <- c("Victivallis vadensis")  
    phylogeny.filtered <- phylogeny.filtered[!phylogeny.filtered$species %in% rm.species,]

    # Remove probes that target multiple L1 groups
    hits <- table(unique(phylogeny.filtered[, c("oligoID", "L1")])$oligoID)
    oligoID <- NULL
    phylogeny.filtered <- subset(phylogeny.filtered, oligoID %in% names(which(hits <= 1)))

    # The standard database query returns 3631 unique oligoIDs for HITChip after explicitly excluding 
    # 'UNI 515', 'HIT 5658', 'HIT 1503', 'HIT 1505', 'HIT 1506'
    # Then standard filters:
    # pmTm.margin = 2.5 (260 oligoIDs discarded);
    # complement = 1 (1 oligoID discarded);
    # mismatch = 0 (260 oligoIDs discarded; only partially overlapping with other filters);
    # -> The filtered phylogeny is used for species/L1/L2 summarization
    # -> The full phylogeny is still OK for oligo-level analyses, as filtering controls mainly
    #    for mismatches but otherwise the oligos are valid and indeed target some taxa that are
    #    missing from the given phylogeny

    # Keep only relevant cols
    phylogeny.full <- phylogeny.full[, 1:6]
    phylogeny.filtered <- phylogeny.filtered[, 1:6]

    # Remove duplicate rows
    phylogeny.full <- phylogeny.full[!duplicated(phylogeny.full),]
    phylogeny.filtered <- phylogeny.filtered[!duplicated(phylogeny.filtered),]

  list(full = phylogeny.full, filtered = phylogeny.filtered)

}
