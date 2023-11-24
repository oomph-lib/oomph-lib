#! /bin/bash

echo $#
echo $1
if [ $# -gt 2 -o "$1" == "--help" ] ; then
 echo " "
 echo " "
 echo "Usage: "
 echo " "
 echo "./emacs_files_that_contain_string.bash"
 echo "    Uses hard-coded search string and search root directory (demo_drivers)"
 echo " "
 echo "./emacs_files_that_contain_string.bash search_string"
 echo "    Searches for search string in default search root directory (demo_drivers)" 
echo " "
 echo "./emacs_files_that_contain_string.bash search_string search_root_dir"
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
search_root_dir=.
if [ $# -eq 2 ]; then
    search_root_dir=$2
fi

echo " "
echo "====================================================="
echo "Opening all files in all directories below"
echo " "
echo "          "$search_root_dir
echo " "
echo "where a file (*.cc or *.h) contains the string"
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

#echo "Hit return to start"
#read continue_flag

#echo "running"
#exit


# Find all *.h and *.cc files that include the search string
#==========================================================

# Clean up
rm -f .tmp_file_name.list

# Make script to extract directories where the codes live
find $search_root_dir \( -name '*.cc' -o  -name '*.h' \) -exec grep -H -m 1 -l  $search_string {} \; > .tmp_file_name.list


#cat .tmp_file_name.list



# Now visit all files
for file in `cat .tmp_file_name.list`; do
    echo " "
    echo "=========================================================="
    echo " "
    echo " "
    echo "Opening: " $file
    emacs $file
    echo " "
    echo " "
    echo "=========================================================="
    echo " "
    echo " "
done

# Clean up
rm -f .tmp_file_name.list
