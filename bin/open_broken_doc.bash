#! /bin/bash

echo " "
echo "===================================================="
echo "I'm going to use the script "
echo " "
echo "    bin/find_missing_doxygen_hooks.sh "
echo " "
echo "to find potentially broken documentation. "
echo "For each potentially broken html file, I'll display"
echo "the file (using an in-script specificable brower)"
echo "and try to open the *.txt that's likely be the source"
echo "for the documentation. If it can't be found that way,"
echo "you'll have to dig around yourself."
echo " "
echo "If you need the corresponding source code, follow the"
echo "link at the bottom of the html file in the browser"
echo "(for most tutorials)."
echo "======================================================"

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
        echo "I'm assuming that the dodgy-looking html code"
        echo " "
        echo "        "$file
        echo " "
        echo "has been created from "
        echo " "
        echo "        "$txt_file
        echo " "
        firefox $file
        emacs $txt_file
    else
        echo " "
        echo "No *.txt file in $txt_file_dir; you'll have to find"
        echo "the txt file that was used to create this documentation"
        echo "yourself."
        echo " "
    fi
done
