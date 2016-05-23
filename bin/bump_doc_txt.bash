#! /bin/bash

echo " " 
echo "Bumping revision information in doc.txt by adding dummy comment "
echo "at end of file..."
echo " "


if [ ! -e doc/doc.txt ]; then
    echo "File doc/doc.txt not found. This script should be run "
    echo "from oomph-lib home directory!"
    exit 1
else
    # Check that doc.txt does not have any other changes so far
    bla=`svn diff doc/doc.txt | wc -c`
    if [ $bla -ne 0 ]; then
        echo " " 
        echo "ERROR:" 
        echo " " 
        echo "doc/doc.txt has already been changed; please commmit these"
        echo "changes first. The changes are:"
        echo " " 
        svn diff doc/doc.txt
        echo " " 
        exit 1
    else
        # Bump and auto-commit
        cd doc; echo "#dummy" > junk.txt; cat doc.txt junk.txt > doc2.txt; mv doc2.txt doc.txt; rm junk.txt
        svn commit doc.txt -m "Bumped doc/doc.txt revision number by running bin/bump_doc_txt.bash"
    fi

fi

echo "...done"

