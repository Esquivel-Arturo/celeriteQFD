#!/bin/bash

tar -xzf Data.tar.gz
tar -xzf CeleriteQFD.tar.gz

mkdir ./res-$1
# make sure the script will use your R installation, 
# and the working directory as its home location
mkdir ./packages
export R_LIBS=$PWD/packages

# run your script
Rscript ./CeleriteQFD/tests/real_data_slicing/QFD-all_CHTC.R $1

tar -czvf QFD-all-$1.tar.gz ./res-$1
