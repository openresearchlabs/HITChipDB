#!/bin/sh
/usr/lib/R CMD BATCH document.R
/usr/lib/R CMD build ../../
/usr/lib/R CMD check HITChipDB_0.6.30.tar.gz
/usr/lib/R CMD INSTALL HITChipDB_0.6.30.tar.gz



