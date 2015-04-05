<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{HITChipDB tutorial}
  %\usepackage[utf8]{inputenc}
-->



HITChipDB R package
===================

Installation
------------

### Install/upload and loading the release version

It is advisable to run this before every analysis to use the latest
version:

    install.packages("devtools")
    library(devtools)
    install_github("microbiome/HITChipDB")

### Extracting HITChip data from the MySQL database ('profiling')

[Instructions to extract data from database](wurcomputer)
(HIT/MIT/PIT/ChickChip)

### Oligo heatmap

To reproduce the oligo-level heatmap from the profiling script, modify
[this example](Oligoheatmap.md)

### Read data from SQL database

Read MySQL data for specified projects (request usernames, passwords,
etc. from the admins):

    library(HITChipDB) 

    proj <- c("MetaHIT") # List projects to be extracted; see list.mysql.projects for a complete list
    dbuser = "myusername"; # Request username from the admins
    dbpwd = "mypasswd";    # Request password from the admins
    host <- NULL; # Used with HITChip FTP server; ask details from admins
    port <- NULL; # Used with HITChip FTP server; ask details from admins

    # Get sample information matrix for the selected projects   
    project.info <- fetch.sample.info(proj, dbuser = dbuser, dbpwd = dbpwd, dbname = "Phyloarray", host = host, port = port)

### List projects in MySQL database

    projs <- list.mysql.projects(dbuser, dbpwd, dbname, host = NULL, port = NULL)

Retrieving microarray data
--------------------------

**HITChip** Use virtual machine or WUR database computer through FTP
connection to Helsinki database. Once you have opened the connection,
see the instructions above on how to read the data in R.

To access HITChip FTP database from your own computer, install the
virtual machine. For further instructions, contact the admins.

After installing the virtual machine, see R instructions for HITChip.

**MITChip, PITChip, ChickChip** Use WUR database computer
([instructions](wurcomputer))

### Probe-level operations

[Probe-level preprocessing](Probelevel.md)

### Robust Probabilistic Averaging

[RPA](RPA.md)

### Further SQL functions to be documented

The following functions are also available:

-   get.phylogeny.info: from the MySQL db
-   fetch.projects: Fetch projects from the phyloarray MySQL database
-   fetch.samples: Fetch samples from the phyloarray MySQL database
-   fetch.sample.info: Fetch sample information from HITChip atlas

### Licensing and Citations

This work can be freely used, modified and distributed under the
[Two-clause FreeBSD license](http://en.wikipedia.org/wiki/BSD_licenses).

Kindly cite the work as 'Leo Lahti and Jarkko Salojarvi (2014).
microbiome R package. URL: <http://microbiome.github.com>'.

### Session info

This vignette was created with

    sessionInfo()

    ## R version 3.1.2 (2014-10-31)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] tcltk     parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] HITChipDB_0.5.13     RPA_1.20.01          affydata_1.12.0     
    ##  [4] affy_1.42.3          RMySQL_0.10.1        microbiome_0.99.35  
    ##  [7] AnnotationDbi_1.26.1 GenomeInfoDb_1.0.2   Biobase_2.24.0      
    ## [10] BiocGenerics_0.10.0  RSQLite_1.0.0        DBI_0.3.1           
    ## [13] reshape_0.8.5        vegan_2.2-1          lattice_0.20-29     
    ## [16] permute_0.8-3        e1071_1.6-4          ade4_1.6-2          
    ## [19] rmarkdown_0.5.1     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] acepack_1.3-3.3       affyio_1.32.0         BiocInstaller_1.14.3 
    ##  [4] class_7.3-11          cluster_1.15.3        codetools_0.2-9      
    ##  [7] colorspace_1.2-4      df2json_0.0.2         digest_0.6.8         
    ## [10] doParallel_1.0.8      dynamicTreeCut_1.62   evaluate_0.5.5       
    ## [13] fastcluster_1.1.15    foreach_1.4.2         foreign_0.8-61       
    ## [16] formatR_1.0           Formula_1.2-0         ggplot2_1.0.0        
    ## [19] GO.db_2.14.0          grid_3.1.2            gtable_0.1.2         
    ## [22] Hmisc_3.14-6          htmltools_0.2.6       igraph_0.7.1         
    ## [25] impute_1.38.1         IRanges_1.22.10       iterators_1.0.7      
    ## [28] knitr_1.9             latticeExtra_0.6-26   MASS_7.3-37          
    ## [31] Matrix_1.1-5          matrixStats_0.13.1    mgcv_1.8-3           
    ## [34] mixOmics_5.0-3        munsell_0.4.2         nlme_3.1-119         
    ## [37] nnet_7.3-8            pheatmap_0.7.7        plyr_1.8.1           
    ## [40] preprocessCore_1.26.1 proto_0.3-10          RColorBrewer_1.1-2   
    ## [43] Rcpp_0.11.4           reshape2_1.4.1        RGCCA_2.0            
    ## [46] rgl_0.95.1201         rjson_0.2.15          R.methodsS3_1.6.1    
    ## [49] rpart_4.1-8           scales_0.2.4          splines_3.1.2        
    ## [52] stats4_3.1.2          stringr_0.6.2         survival_2.37-7      
    ## [55] tools_3.1.2           WGCNA_1.43            yaml_2.1.13          
    ## [58] zlibbioc_1.10.0
