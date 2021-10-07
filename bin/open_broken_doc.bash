#! /bin/bash

# Specify browser
browser=firefox

# Find oomph-lib's bin directory (this file lives in it!)
bin_dir=`dirname "$0"`

# Extract the broken html files from it
rm broken_files.txt
$bin_dir/find_missing_doxygen_hooks.sh  | grep 'html:' | awk 'BEGIN{ FS= ":"}{print $1}' | unique broken_files.txt

for file in `cat broken_files.txt`; do
    echo $file
    txt_file_dir=`dirname $file`/..
    txt_file=`echo $txt_file_dir/*.txt`
    if [ -e $txt_file ] ; then
        firefox $file
        emacs $txt_file
    else
        echo "No txt file in $txt_file_dir; you'll have to find that file yourself."
    fi
done
