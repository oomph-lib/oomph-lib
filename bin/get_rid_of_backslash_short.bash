#! /bin/bash

# Get rid of all \short annotations in code; no longer needed
# now that we've set
#
#  MULTILINE_CPP_IS_BRIEF = YES
#  ALWAYS_DETAILED_SEC    = YES
#
# in all Doxyfiles.

file_list=`find . -name '*.h' -o -name '*.cc'`
for file in `echo $file_list`; do
    echo "Doing $file"
    cp $file $file.tmp
    sed 's/\\short//g' $file.tmp > $file
    rm $file.tmp
done

