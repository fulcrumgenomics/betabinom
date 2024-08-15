#!/bin/bash

# update md5sums
md5sum \
    DESCRIPTION \
    LICENSE \
    NAMESPACE \
    R/bb.test.R \
    R/ibb.test.R \
    R/fold.change.R \
    R/normalize.R \
    R/utils.R \
    build/vignette.rds \
    inst/doc/betabinom.R \
    inst/doc/betabinom.Rmd \
    inst/doc/betabinom.html \
    man/bb.test.Rd \
    man/fold.change.Rd \
    man/ibb.test.Rd \
    man/normalize.Rd \
    src/Makevars \
    src/bb.c \
    src/ibb.c \
    src/q15.h \
    src/registerDynamicSymbol.c \
    tests/testthat.R \
    tests/testthat/test-normalize.R \
    tests/testthat/test-fold.change.R \
    tests/testthat/test-bb.test.R \
    tests/testthat/test-ibb.test.R \
    tests/testthat/testdata/example-3groups.txt \
    tests/testthat/testdata/example-paired.txt \
    vignettes/betabinom.Rmd \
    | sed -e "s/  / */g" > MD5
cd ..
