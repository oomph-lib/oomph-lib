#! /bin/bash

#======================================================
# Make all four versions of the distribution
#======================================================

echo "LIKELY TO BE BROKEN AFTER DAVID/MATTHIAS'S REWRITE; NEED TO UPDATE"
echo "THE RESPONSES"
exit

#--------------------
# Adjust number
#--------------------
revision=0.90

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


exit

#------------------------------------------------------
# Directory for destruct tests
#------------------------------------------------------
dir=test_all_distributions
#rm -rf $dir 
if [ -e $dir   ] 
then
   echo "Please delete directory $dir and try again"
   exit
fi
mkdir $dir
cd $dir
mkdir no_doc_no_validata
cp ../`echo "oomph-lib-"$revision"_no_doc_no_validata.tar.gz"` no_doc_no_validata
mkdir with_doc_no_validata
cp ../`echo "oomph-lib-"$revision"_no_validata.tar.gz"`        with_doc_no_validata
mkdir no_doc_with_validata
cp ../`echo "oomph-lib-"$revision"_no_doc.tar.gz"`             no_doc_with_validata
mkdir with_doc_with_validata
cp ../`echo "oomph-lib-"$revision".tar.gz"`             with_doc_with_validata

#------------------------------------------------------
# Create build script
#------------------------------------------------------
echo "#! /bin/bash" > build_script.bash
chmod a+x build_script.bash
echo "tar xvfz *.tar.gz " >> build_script.bash  
echo "cd "`echo "oomph-lib-"$revision` >> build_script.bash            
echo "./autogen.sh < ../../replies.txt > test_build.log " >> build_script.bash
echo "cd .. " >> build_script.bash

# Replies to queries in autogen.sh:
echo "y" > replies.txt # confirm build dir
echo "y" >> replies.txt # doc size
echo "y" >> replies.txt # config options are ok
echo "y" >> replies.txt # run the self tests
echo " " >> replies.txt # hit return to start the build process




#------------------------------------------------------
# Build and test
#------------------------------------------------------
cd no_doc_no_validata
../build_script.bash &
cd ..

cd with_doc_no_validata
../build_script.bash &
cd ..

cd no_doc_with_validata
../build_script.bash &
cd ..

cd with_doc_with_validata
../build_script.bash &
cd ..
