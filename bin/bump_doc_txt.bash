#! /bin/bash


files_to_be_bumped="doc/doc.txt doc/the_distribution/the_distribution.txt"
for file in `echo $files_to_be_bumped`; do

    echo " " 
    echo "Bumping revision information in "
    echo " " 
    echo "     $file "
    echo " " 
    echo "by adding dummy comment at end of file."
    echo " "
    
if [ ! -e $file ]; then
    echo "File $file not found. This script should be run "
    echo "from oomph-lib home directory!"
    exit 1
else
    # Check that doc.txt does not have any other changes so far
    bla=`svn diff $file | wc -c`
    if [ $bla -ne 0 ]; then
        echo " " 
        echo "ERROR:" 
        echo " " 
        echo "$file has already been changed; please commmit these"
        echo "changes first. The changes are:"
        echo " " 
        svn diff $file
        echo " " 
        exit 1
    else
        # Bump
        echo "#dummy" > junk.txt; cat $file junk.txt > tmp.txt; mv tmp.txt $file; rm junk.txt
    fi

fi


done

# If we're still going here commit
echo "About to commit..."
svn commit $files_to_be_bumped -m "Bumped revision numbers in $files_to_be_bumped by running bin/bump_doc_txt.bash"
echo "...done"

exit
