#!/bin/bash

# Detect OS and correct link
if [ "$(uname)" == "Darwin" ]; then
    # Do something under Mac OS X platform
    filen=http://www.niehs.nih.gov/research/resources/assets/docs/artbinvanillaicecream031114macos64tgz.tgz
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # Do something under Linux platform
    filen=http://www.niehs.nih.gov/research/resources/assets/docs/artbinvanillaicecream031114linux64tgz.tgz
fi

# Make and change to lib folder
mkdir lib
cd lib

# Get file and uncompress
curl -O $filen
fileb=${filen##*/}
tar zxvf $fileb
