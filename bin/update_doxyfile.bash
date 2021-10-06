#!/bin/bash

# Update all Doxyfiles

# Find 'em
dir_list=`find . -name Doxyfile -exec dirname {} \;`
home_dir=`pwd`
for dir in `echo $dir_list`; do
    cd $dir
    echo "Updating Doxyfile in " `pwd`
    #cp Doxyfile Doxyfile.tmp
    #sed 's/a4wide/a4/g' Doxyfile.tmp > Doxyfile
    #rm Doxyfile.tmp
    doxygen -u
    cd $home_dir
done

