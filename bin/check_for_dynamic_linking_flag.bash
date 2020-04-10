#! /bin/bash

if [ "$#" -ne 1 ]; then
    echo "This script requires the specification of the (temporary)"
    echo "configure options file. You specified $# arguments!"
    exit 1
fi

#echo "Checking if temporary configure file:"
#echo " "
#echo "    " $1
#echo " " 
#echo "requests the dynamic linking for trilinos via specification "
#echo "of --enable-dynamic-linking-for-trilinos flag."
#echo " " 

# Check if flag specifying dynamic linking appears in config options
# (starting at beginning of line: ^). Input argument = temporary
# configure options file. (For some reason this test doesn't work if
# run directly in non-interactive_autogen.sh)
count_match=`grep -c -m 1 '^\-\-enable-dynamic-linking-for-trilinos flag' $1`

#echo $count_match

if [ "$count_match" -ne "0" ]; then
    echo "dynamic"
else
    echo "static"
fi


