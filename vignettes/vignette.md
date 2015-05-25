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

[Instructions to extract data from database](wurcomputer.md)
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

[Robust Probabilistic Averaging
(RPA)](https://github.com/antagomir/RPA/wiki)

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

    ## R version 3.2.0 (2015-04-16)
    ## Platform: x86_64-unknown-linux-gnu (64-bit)
    ## Running under: Ubuntu 15.04
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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] rmarkdown_0.6.1    scimapClient_0.2.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] magrittr_1.5    formatR_1.2     htmltools_0.2.6 tools_3.2.0    
    ##  [5] yaml_2.1.13     RJSONIO_1.3-0   stringi_0.4-1   knitr_1.10.5   
    ##  [9] stringr_1.0.0   digest_0.6.8    evaluate_0.7
