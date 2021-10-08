#!/bin/bash

# Update all Doxyfiles

# Find 'em
dir_list=`find . -name Doxyfile -exec dirname {} \;`
home_dir=`pwd`
for dir in `echo $dir_list`; do
    cd $dir
    echo "Updating Doxyfile in " `pwd`
    # Keep this around (commented out) because sed syntax is clunky
    # and I can never remember it... [but see below]
    #cp Doxyfile Doxyfile.tmp
    #sed 's/a4wide/a4/g' Doxyfile.tmp > Doxyfile
    #rm Doxyfile.tmp
    doxygen -u
    cd $home_dir
done

exit

##########################################################################

# Here's a better version proposed by Puneet; don't want to test this
# now. Do it next time we actually use this script in anger

#!/bin/bash

# Control the script behaviour using these variables
enable_sed=false
original_string='a4wide'
replacement_string='a4'

# Helper function to handle the typically-messy escaping of characters
escape_sed_args() {
    echo "$1" | sed 's/[\/&]/\\&/g'
}

# Escape characters
original_string=$(escape_sed_args "$original_string")
replacement_string=$(escape_sed_args "$replacement_string")

# Update all Doxyfiles

# Find 'em
home_dir="$(pwd)"
dir_list=$(find . -name Doxyfile -exec dirname {} \;)
for dir in $(echo $dir_list); do
    cd $dir

    # Update the Doxyfile before using 'sed' otherwise the changes get overwritten
    echo "Updating Doxyfile: ${dir}/Doxyfile"

    # Modify certain strings, if desired
    if [[ "$enable_sed" == true ]]; then
        echo "Running sed with args: 's/${original_string}/${replacement_string}/g'"
        sed "s/${original_string}/${replacement_string}/g" Doxyfile >Doxyfile.tmp
        mv -f Doxyfile.tmp Doxyfile
    fi

    doxygen -u
    cd $home_dir
done
