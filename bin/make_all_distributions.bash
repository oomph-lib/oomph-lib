#! /bin/bash

#======================================================
# Make all four versions of the distribution
#======================================================

#--------------------
# Adjust number
#--------------------
revision=0.85

#--------------------
# No doc, no validata
#--------------------
echo 5 > /tmp/junk.txt
make mydist < /tmp/junk.txt
mv `echo "oomph-lib-"$revision".tar.gz"` `echo " oomph-lib-"$revision"_no_doc_no_validata.tar.gz"`



#--------------------
# With doc, no validata
#--------------------
echo 4 > /tmp/junk.txt
make mydist < /tmp/junk.txt
mv `echo "oomph-lib-"$revision".tar.gz"` `echo " oomph-lib-"$revision"_no_validata.tar.gz"`



#----------------------
# No doc, with validata
#----------------------
echo 3 > /tmp/junk.txt
make mydist < /tmp/junk.txt
mv `echo "oomph-lib-"$revision".tar.gz"` `echo " oomph-lib-"$revision"_no_doc.tar.gz"`



#--------------------
# Full distribution;
# no need to rename
#--------------------
echo 2 > /tmp/junk.txt
make mydist < /tmp/junk.txt

