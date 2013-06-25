#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=0

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for circular driven cavity
#-----------------------------------------
cd Validation

echo "Running matrix multiply validation "
mkdir RESLT
../matrix_multiply_test> OUTPUT
echo "done"
echo " " >> validation.log
echo "Matrix multiply validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cd RESLT
#cat .. > results.dat
cd ..



#if test "$1" = "no_fpdiff"; then
#  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
#else
#../../../bin/fpdiff.py ../validata/results.dat.gz RESLT/results.dat >> validation.log
#fi



# Append log to main validation log
cat validation.log >> ../../../validation.log

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
