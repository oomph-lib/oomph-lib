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

# Validation for rectangular driven cavity
#-----------------------------------------
cd Validation

echo "Running two d biharmonic validation "
mkdir RESLT
../two_d_biharmonic > OUTPUT_two_d_biharmonic
echo "done"
echo " " >> validation.log
echo "Two d biharmonic validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/biharmonic1_results.dat.gz  \
         RESLT/soln_0.dat 0.1 1.0e-12 >> validation.log
../../../../bin/fpdiff.py ../validata/biharmonic2_results.dat.gz  \
         RESLT/soln_1.dat 0.1 1.0e-12 >> validation.log
fi

mv RESLT RESLT_two_d_biharmonic

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
