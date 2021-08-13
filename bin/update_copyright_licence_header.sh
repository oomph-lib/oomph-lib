#!/usr/bin/env bash

# exit when any command fails
set -e

#-------------------------------------------------------------------------------
# Shell script to update the copyright notice of all oomph-lib-owned files.
#-------------------------------------------------------------------------------

# --------[ IDENTIFY COPYRIGHTED FILES ]----------------------------------------

# Get the current filename; we need to ignore it from the search rule
current_filename=$(basename "$0")

# Search for all files containing the copyright notice
oomph_lib_copyrighted_files=$(grep . -R -l --exclude="$current_filename" -e 'Copyright (C) 2006-' -e ' Matthias Heil and Andrew Hazel')

# --------[ UPDATE COPYRIGHT YEAR ]---------------------------------------------

# The end-year to update the notice to
current_year=$(date +"%Y")

# Update the licence header for each file. Don't use the "-i" flag for an
# in-place edit as this flag is not supported by the BSD version of "sed" that
# ships with macOS
for file in $oomph_lib_copyrighted_files; do
    echo "Updating copyright notice for file:" $file
    sed "s/Copyright (C) 2006-.* Matthias Heil and Andrew Hazel/Copyright (C) 2006-$current_year Matthias Heil and Andrew Hazel/g" $file >$file.updated
    mv $file.updated $file
done

# ------------------------------------------------------------------------------
