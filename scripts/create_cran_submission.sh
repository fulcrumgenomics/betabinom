#!/bin/bash

# clean previous CRAN build
rm -rf betabinom

# move new vignette build to expected location
mv Meta/* build
mv doc/* inst/doc

# create new CRAN build
mkdir betabinom
cp -r build betabinom
cp -r inst betabinom
cp -r man betabinom
cp -r R betabinom
cp -r vignettes betabinom
cp -r src betabinom # copy source (without binary blobs)
rm betabinom/src/*.o
rm betabinom/src/*.so
cp metadata/DESCRIPTION betabinom
cp metadata/NAMESPACE betabinom
cp metadata/LICENSE betabinom
cp metadata/MD5 betabinom

# run all CRAN checks
R CMD CHECK --as-cran betabinom
