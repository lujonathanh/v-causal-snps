#!/bin/bash

# Make sure to update package_params.sh first!

source ./package_params.sh

mkdir $FOLDER

cp $DATAFILE $FOLDER
cp $ARGSFILE $FOLDER

while read file; do
    echo Copy $file to $FOLDER
    cp $file $FOLDER
done <package_required_files.txt


