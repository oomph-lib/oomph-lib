#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for circular driven cavity
#-----------------------------------------
cd Validation

echo "Running circular driven cavity validation "
mkdir RESLT
../circular_driven_cavity > OUTPUT_circular_driven_cavity
echo "done"
echo " " >> validation.log
echo "Circular driven cavity validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cd RESLT
cat  soln0.dat   soln11.dat  soln14.dat  soln17.dat  soln2.dat  soln5.dat  \
     soln8.dat soln1.dat   soln12.dat  soln15.dat  soln18.dat  soln3.dat  \
     soln6.dat  soln9.dat soln10.dat  soln13.dat  soln16.dat  soln19.dat  \
     soln4.dat  soln7.dat > circular_driven_cavity_results.dat
cd ..



if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/circular_driven_cavity_validation.dat.gz  \
         RESLT/circular_driven_cavity_results.dat >> validation.log
fi

mv RESLT RESLT_circular_driven_cavity



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
