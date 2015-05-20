/usr/bin/R CMD BATCH document.R
/usr/bin/R CMD build ../../
/usr/bin/R CMD check HITChipDB_0.5.16.tar.gz
/usr/bin/R CMD INSTALL HITChipDB_0.5.16.tar.gz


