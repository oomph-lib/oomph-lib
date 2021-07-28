#! /bin/sh

# exit when any command fails
set -e

#-------------------------------------------------------------------------------
# Shell script to perform the following steps (in the described order):
# 1. Remove licence (all lines starting with //LIC//) from all oomph-lib *.h
#    and *.cc files.
# 2. Update end-year of the copyright licence in "cc_and_h_licence_block.txt".
# 3. Apply licence prefix to all oomph-lib *.h and *.cc files.
#-------------------------------------------------------------------------------

# Find *.h and *.cc files in oomph-lib distribution (i.e. exclude external_src)
oomph_lib_h_and_cc_files=$(find demo_drivers src user_src user_drivers self_test bin \( -name '*.h' -o -name '*.cc' \) -exec ls {} \;)

# --------[ REMOVE LICENSE ]----------------------------------------------------

# Remove the licence header from each source file
for file in $oomph_lib_h_and_cc_files; do
    echo "Delicencing:" $file
    gawk -f bin/remove_licence.awk $file >$file.delicenced
    mv -f $file.delicenced $file
done

# --------[ UPDATE COPYRIGHT YEAR ]---------------------------------------------

current_year=$(date +"%Y")
sed --in-place "s/Copyright (C) 2006-.* Matthias Heil and Andrew Hazel/Copyright (C) 2006-$current_year Matthias Heil and Andrew Hazel/g" bin/cc_and_h_licence_block.txt

# --------[ APPLY LICENSE ]-----------------------------------------------------

# Add the licence header back to delicenced files
for file in $oomph_lib_h_and_cc_files; do
    echo "Licencing:" $file
    cat bin/cc_and_h_licence_block.txt $file >$file.licenced
    gawk -f bin/remove_licence.awk $file.licenced >$file.delicenced
    if (test $(diff $file $file.delicenced | wc -w) != 0); then
        echo "Original and licenced/delicenced file don't match!"
        echo "Please run bin/remove_licence.sh first."
        echo "Filename: " $file
        exit 1
    else
        echo "OK"
        mv -f $file.licenced $file
        rm $file.delicenced
    fi
done
