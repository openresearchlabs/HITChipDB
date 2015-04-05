### Oligoprofile heatmap

Reproduce and modify the oligo heatmap of the
[run.profiling.script](profiling). Plots heatmap of the oligo profiles
together with phylotype groups and sample clusters. Reload and plot
preprocessed data from the [run.profiling.script](profiling) by using
the [read.profiling](reading) function, assuming you have stored the
output in directory "datadirectory/". Start by reading the data:

    library(HITChipDB)

    ## Loading required package: ade4
    ## Loading required package: microbiome
    ## Loading required package: e1071
    ## Loading required package: vegan
    ## Loading required package: permute
    ## Loading required package: lattice
    ## This is vegan 2.2-1
    ## 
    ## Attaching package: 'vegan'
    ## 
    ## The following object is masked from 'package:ade4':
    ## 
    ##     cca
    ## 
    ## Loading required package: reshape
    ## Loading required package: DBI
    ## Loading required package: AnnotationDbi
    ## Loading required package: BiocGenerics
    ## Loading required package: parallel
    ## 
    ## Attaching package: 'BiocGenerics'
    ## 
    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     xtabs
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, as.vector, cbind,
    ##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
    ##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unlist
    ## 
    ## Loading required package: Biobase
    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.
    ## 
    ## Loading required package: GenomeInfoDb
    ## 
    ## 
    ## microbiome R package (microbiome.github.com)
    ##           
    ## 
    ## 
    ##  Copyright (C) 2011-2015
    ##           Leo Lahti and Jarkko Salojarvi 
    ## 
    ##         
    ##           <microbiome-admin@googlegroups.com>
    ## 
    ## 
    ## Attaching package: 'microbiome'
    ## 
    ## The following object is masked from 'package:lattice':
    ## 
    ##     densityplot
    ## 
    ## The following object is masked from 'package:e1071':
    ## 
    ##     impute
    ## 
    ## Loading required package: RMySQL
    ## 
    ## Attaching package: 'RMySQL'
    ## 
    ## The following object is masked from 'package:RSQLite':
    ## 
    ##     isIdCurrent
    ## 
    ## Loading required package: RPA
    ## Loading required package: affy
    ## Loading required package: affydata

    ##      Package    LibPath                                            
    ## [1,] "affydata" "/home/antagomir/R/x86_64-pc-linux-gnu-library/3.1"
    ##      Item       Title                        
    ## [1,] "Dilution" "AffyBatch instance Dilution"

    ## 
    ## RPA Copyright (C) 2008-2013 Leo Lahti.
    ## This program comes with ABSOLUTELY NO WARRANTY.
    ## This is free software, and you are welcome to redistribute it under the FreeBSD open source license.
    ## 
    ## Loading required package: tcltk
    ## 
    ## HITChipDB R package (microbiome.github.com)
    ## (C) 2011-2015 Leo Lahti and Jarkko Salojarvi <microbiome-admin@googlegroups.com>
    ## This program comes with ABSOLUTELY NO WARRANTY.
    ## This is free software, and you are welcome to redistribute it under the FreeBSD open source license.

    library(microbiome)

    # Define data directory (here: simulated data directory)
    data.directory <- system.file("extdata", package = "microbiome")

    # Read Oligo level data in original domain
    oligo.matrix.log10.simulated <- read.profiling(level = "oligo", 
                              data.dir = data.directory, log10 = TRUE)

    ## Reading /home/antagomir/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata/oligoprofile.tab
    ## Logarithmizing the data

    # Read Oligo-phylogeny mapping table (two methods):
    phylogeny.info <- GetPhylogeny("HITChip", "full")

    ## Reading /home/antagomir/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata/phylogeny.full.tab

Save the image in a file (too large to be opened directly in R). To
prevent ordering of the rows, use hclust.method = NULL in the function
call:

    library(microbiome)
    ppcm <- 150
    png(filename = "oligoprofileClustering.png", 
         width = max(trunc(ppcm*21), 
                     trunc(ppcm*21*ncol(oligo.matrix.log10.simulated)/70)), 
         height = trunc(ppcm*29.7))


    library(HITChipDB)
    tmp <- PlotPhylochipHeatmap(oligo.matrix.log10.simulated, 
                    phylogeny.info, level = "L1", metric = "pearson")
    dev.off()
