#!/bin/bash

# clean previous CRAN build
rm -rf betabinom

# create new CRAN build
mkdir betabinom
cp -r build betabinom
cp -r inst betabinom
cp -r man betabinom
cp -r R betabinom
cp -r src betabinom
rm betabinom/src/*.o
rm betabinom/src/*.so
cp -r vignettes betabinom
cp DESCRIPTION betabinom
cp LICENSE betabinom
cp MD5 betabinom
cp NAMESPACE betabinom
cp cran-comments.md betabinom

# run all CRAN checks
R CMD CHECK --as-cran betabinom
