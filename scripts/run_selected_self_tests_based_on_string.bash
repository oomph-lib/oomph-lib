#! /bin/bash

echo $#
echo $1
if [ $# -gt 2 -o "$1" == "--help" ] ; then
 echo " "
 echo " "
 echo "Usage: "
 echo " "
 echo "./run_selected_self_tests_based_on_string.bash"
 echo "    Uses hard-coded search string and search root directory (demo_drivers)"
 echo " "
 echo "./run_selected_self_tests_based_on_string.bash search_string"
 echo "    Searches for search string in default search root directory (demo_drivers)" 
echo " "
 echo "./run_selected_self_tests_based_on_string.bash search_string search_root_dir"
 echo "    Searches for search string in specified search root directory"
 echo " "
 echo " "
 exit 
fi

#Search string 
search_string=solve_eigenproblem
if [ $# -ge 1 ]; then
    search_string=$1
fi

#Search root directory
search_root_dir=demo_drivers
if [ $# -eq 2 ]; then
    search_root_dir=$2
fi

echo " "
echo "====================================================="
echo "Running self-tests (make -k) in all directories below"
echo " "
echo "          "$search_root_dir
echo " "
echo "where a demo driver code (*.cc) contains the string"
echo " "
echo "          "$search_string
echo " "
echo "====================================================="
echo " "

if [ ! -d $search_root_dir ]; then
    echo " "
    echo "Sorry; "$search_root_dir" is not a directory"
    echo "bailing out..."
    echo " "
    exit
else
    echo " "
    echo "Yeah; "$search_root_dir" is a directory"
    echo " "
fi

if [ -e partial_validation.log ]; then
    echo " "
    echo "File partial_validation.log already exists."
    echo "Current run will append. Continue (y/n)"
    read continue_flag
    if [ "$continue_flag" != "y" ]; then
        echo "bailing out..."
        exit
    fi
else
    echo " "
    echo "Result of self-tests (concatenated validation.log)"
    echo "will be written to partial_validation.log"
    echo " "
fi

#echo "Hit return to start"
#read continue_flag

#echo "running"
#exit


# Find all demo driver codes that include the search string
#==========================================================

# Clean up
rm -f .tmp_dir_name.bash .tmp_junk.list .tmp_unique_junk.list

# Make script to extract directories where the codes live
find $search_root_dir -name '*.cc' -exec grep -H -m 1 -l  $search_string {} \; | awk '{print "echo `dirname "$1"` >> .tmp_junk.list"}' > .tmp_dir_name.bash


cat .tmp_dir_name.bash

# Run it
source .tmp_dir_name.bash

# Extract unique directories
cat .tmp_junk.list | unique .tmp_unique_junk.list

# Now visit all directories
current_dir=`pwd`
for dir in `cat .tmp_unique_junk.list`; do
    cd $dir
    echo " "
    echo " "
    echo "=========================================================="
    echo " "
    echo " "
    echo "I'm in: " `pwd`
    echo " "
    echo " "
    echo "=========================================================="
    echo " "
    echo " "
    make -k check
    cat Validation/validation.log >> `echo $current_dir`/partial_validation.log
    cd $current_dir
done

# Clean up
rm -f .tmp_dir_name.bash .tmp_junk.list .tmp_unique_junk.list
