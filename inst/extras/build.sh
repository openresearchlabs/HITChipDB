/usr/local/bin/R CMD BATCH document.R
/usr/local/bin/R CMD build ../../
/usr/local/bin/R CMD check HITChipDB_0.5.16.tar.gz
/usr/local/bin/R CMD INSTALL HITChipDB_0.5.16.tar.gz


