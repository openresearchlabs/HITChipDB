# https://support.rstudio.com/hc/en-us/articles/200486518-Customizing-Package-Build-Options
/usr/bin/R CMD BATCH document.R
/usr/bin/R CMD build ../../ --no-build-vignettes --no-examples
/usr/bin/R CMD check HITChipDB_0.6.31.tar.gz --no-build-vignettes --no-examples
/usr/bin/R CMD INSTALL HITChipDB_0.6.31.tar.gz 
#/use/bin/R CMD BiocCheck microbiome_0.99.32.tar.gz 
