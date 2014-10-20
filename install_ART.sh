#!/bin/bash

# Make and change to lib folder
mkdir lib
cd lib

# Detect OS, correct link and get file
if [ "$(uname)" == "Darwin" ]; then
    filen=http://www.niehs.nih.gov/research/resources/assets/docs/artbinvanillaicecream031114macos64tgz.tgz
    curl -O $filen
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    filen=http://www.niehs.nih.gov/research/resources/assets/docs/artbinvanillaicecream031114linux64tgz.tgz
    wget $filen
fi

# Uncompress
fileb=${filen##*/}
tar zxvf $fileb
