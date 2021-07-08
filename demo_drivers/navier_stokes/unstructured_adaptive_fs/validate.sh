#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for unstructured drop
#---------------------------------
cd Validation

echo "Running 2D unstructured adaptive drop in channel validation" 
mkdir RESLT

../adaptive_drop_in_channel --validation > OUTPUT_drop
echo "done"
echo " " >> validation.log
echo "2D unstructured adaptive drop in channel validation" >> validation.log
echo "---------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/norm.dat > results_int.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    ../../../../bin/fpdiff.py ../validata/results_int.dat.gz  \
        results_int.dat 0.1 6.0e-8 >> validation.log
fi

mv RESLT RESLT_int


echo "Running 2D unstructured adaptive bubble in channel validation" 
mkdir RESLT
../adaptive_bubble_in_channel --validation > OUTPUT_bubble
echo "done"
echo " " >> validation.log
echo "2D unstructured adaptive bubble in channel validation" >> validation.log
echo "-----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/norm.dat > results_fs.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    ../../../../bin/fpdiff.py ../validata/results_fs.dat.gz  \
        results_fs.dat 0.2 5.0e-8 >> validation.log
fi

mv RESLT RESLT_fs



# Append log to main validation log
cat validation.log >> ../../../../validation.log

cd ..



#######################################################################


#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
