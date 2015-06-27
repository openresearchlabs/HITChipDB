---
title: "HITChipDB vignette"
author: "Leo Lahti and Jarkko Salojarvi"
date: "2015-06-27"
output: md_document
---

<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{HITChipDB tutorial}
  %\usepackage[utf8]{inputenc}
-->



HITChipDB R package
===========

## Installation

### Install/upload and loading the release version

It is advisable to run this before every analysis to use the latest version:


```r
install.packages("devtools")
library(devtools)
install_github("microbiome/HITChipDB")
```

### Extracting HITChip data from the MySQL database ('profiling')

[Instructions to extract data from database](wurcomputer.md) (HIT/MIT/PIT/ChickChip)

### Oligo heatmap

To reproduce the oligo-level heatmap from the profiling script, modify [this example](Oligoheatmap.md)

### Read data from SQL database

Read MySQL data for specified projects (request usernames, passwords, etc. from the admins):


```r
library(HITChipDB) 

proj <- c("MetaHIT") # List projects to be extracted; see list.mysql.projects for a complete list
dbuser = "myusername"; # Request username from the admins
dbpwd = "mypasswd";    # Request password from the admins
host <- NULL; # Used with HITChip FTP server; ask details from admins
port <- NULL; # Used with HITChip FTP server; ask details from admins

# Get sample information matrix for the selected projects	
project.info <- fetch.sample.info(proj, dbuser = dbuser, dbpwd = dbpwd, dbname = "Phyloarray", host = host, port = port)
```

### List projects in MySQL database


```r
projs <- list.mysql.projects(dbuser, dbpwd, dbname, host = NULL, port = NULL)
```

## Retrieving microarray data

**HITChip** Use virtual machine or WUR database computer through FTP
  connection to Helsinki database. Once you have opened the
  connection, see the instructions above on how to read the data in R.

  To access HITChip FTP database from your own computer, install the
  virtual machine. For further instructions, contact the admins.

  After installing the virtual machine, see R instructions for
  HITChip.

**MITChip, PITChip, ChickChip** Use WUR database computer
  ([instructions](wurcomputer))



### Further SQL functions 

* get.phylogeny.info: from the MySQL db
* fetch.projects: Fetch projects from the phyloarray MySQL database
* fetch.samples: Fetch samples from the phyloarray MySQL database
* fetch.sample.info: Fetch sample information from HITChip atlas


### Licensing and Citations

This work can be freely used, modified and distributed under the 
[Two-clause FreeBSD license](http://en.wikipedia.org/wiki/BSD\_licenses).

Kindly cite the work as 'Leo Lahti and Jarkko Salojarvi
(2014). microbiome R package. URL: http://microbiome.github.com'.


### Session info

This vignette was created with


```r
sessionInfo()
```

```
## R version 3.2.1 (2015-06-18)
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
## [1] tcltk     parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] HITChipDB_0.5.25      RMySQL_0.10.3         preprocessCore_1.30.0
##  [4] DBI_0.3.1             microbiome_0.99.55    RPA_1.24.0           
##  [7] affy_1.46.1           Biobase_2.28.0        BiocGenerics_0.14.0  
## [10] phyloseq_1.13.2       rmarkdown_0.7         knitr_1.10.5         
## [13] scimapClient_0.2.1   
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_1.2-6          som_0.3-5                
##   [3] futile.logger_1.4.1       XVector_0.8.0            
##   [5] OAIHarvester_0.1-7        RcppArmadillo_0.5.200.1.0
##   [7] GenomicRanges_1.20.5      affyio_1.36.0            
##   [9] AnnotationDbi_1.30.1      codetools_0.2-11         
##  [11] splines_3.2.1             geneplotter_1.46.0       
##  [13] mixOmics_5.0-4            tgp_2.4-11               
##  [15] ade4_1.7-2                spam_1.0-1               
##  [17] Formula_1.2-1             annotate_1.46.0          
##  [19] cluster_2.0.2             pheatmap_1.0.2           
##  [21] Kendall_2.2               sorvi_0.7.26             
##  [23] assertthat_0.1            Matrix_1.2-1             
##  [25] formatR_1.2               acepack_1.3-3.3          
##  [27] htmltools_0.2.6           tools_3.2.1              
##  [29] igraph_0.7.1              rdryad_0.1.1             
##  [31] gtable_0.1.2              reshape2_1.4.1           
##  [33] dplyr_0.4.2               maps_2.3-9               
##  [35] Rcpp_0.11.6               Biostrings_2.36.1        
##  [37] RJSONIO_1.3-0             multtest_2.24.0          
##  [39] biom_0.3.12               gdata_2.16.1             
##  [41] ape_3.3                   nlme_3.1-120             
##  [43] iterators_1.0.7           lmtest_0.9-34            
##  [45] stringr_1.0.0             proto_0.3-10             
##  [47] gtools_3.5.0              XML_3.98-1.2             
##  [49] zlibbioc_1.14.0           MASS_7.3-41              
##  [51] zoo_1.7-12                scales_0.2.5             
##  [53] BiocInstaller_1.18.3      lambda.r_1.1.7           
##  [55] RColorBrewer_1.1-2        fields_8.2-1             
##  [57] yaml_2.1.13               gridExtra_0.9.1          
##  [59] ggplot2_1.0.1             rpart_4.1-9              
##  [61] latticeExtra_0.6-26       stringi_0.5-2            
##  [63] maptree_1.4-7             RSQLite_1.0.0            
##  [65] genefilter_1.50.0         S4Vectors_0.6.0          
##  [67] tseries_0.10-34           foreach_1.4.2            
##  [69] nortest_1.0-3             permute_0.8-4            
##  [71] boot_1.3-16               BiocParallel_1.2.5       
##  [73] chron_2.3-46              GenomeInfoDb_1.4.1       
##  [75] moments_0.14              bitops_1.0-6             
##  [77] rgl_0.95.1247             evaluate_0.7             
##  [79] lattice_0.20-31           plyr_1.8.3               
##  [81] magrittr_1.5              DESeq2_1.8.1             
##  [83] R6_2.0.1                  IRanges_2.2.4            
##  [85] earlywarnings_1.1.21      Hmisc_3.16-0             
##  [87] foreign_0.8-63            mgcv_1.8-6               
##  [89] survival_2.38-2           RCurl_1.95-4.6           
##  [91] nnet_7.3-9                futile.options_1.0.0     
##  [93] KernSmooth_2.23-14        RGCCA_2.0                
##  [95] locfit_1.5-9.1            grid_3.2.1               
##  [97] data.table_1.9.4          vegan_2.3-0              
##  [99] digest_0.6.8              diptest_0.75-7           
## [101] xtable_1.7-4              stats4_3.2.1             
## [103] munsell_0.4.2             quadprog_1.5-5
```




