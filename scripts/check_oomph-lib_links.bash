#! /bin/bash


#==================================================================
# Script to check for broken internal links within the oomph-lib
# documentation
#==================================================================


# Home directory to return to
home_dir=`pwd`

# find all html directories
dir_list=`find -name html`

for dir in `echo $dir_list`; do

    cd $dir

    # Get name of directory that contains the txt files
    cd ..
    actual_dir=`pwd`
    cd -

    # Check URLs in index.html
    wget --spider --force-html -i index.html -o .link_check_junk.txt

    # find the local files (i.e. the ones that wget complains about as
    # having incomplete urls and chop off the period at the end of 
    # filenames
    link_list=`grep -v '#' .link_check_junk.txt | grep incomplete  | awk '{l=length($6)-1; print substr($6,1,l)}'`

    for link in `echo $link_list`; do
        if [ ! -e `echo $link` ] ; then
            echo "Warning: Documentation in "$actual_dir" contains link to "$dir"/"$link " which does not exist."
        fi
    done

    rm -f .link_check_junk.txt
    cd $home_dir



#wget --spider --force-html -i html/index.html 

done