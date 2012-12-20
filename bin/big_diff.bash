#! /bin/bash
#----------------------------------------------
# Bash script to document the svn changes
# for all *.cc and *.h files
#----------------------------------------------

# Get list of modified files
modified_list=`svn status  | awk '{if ($1="M")  system("echo " $2)} ' | awk '/.*\.h/{print $1} /.*\.cc/{print $1}'`


# Loop over files
for file in `echo $modified_list`; do

    echo -e "\n\n\n\n##########################################################\n\nFile: " $file"\n\n##########################################################\n\n: "
#    svn diff  -x -w $file
    svn diff $file
done
