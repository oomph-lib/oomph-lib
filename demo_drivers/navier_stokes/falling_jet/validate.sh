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

# Validation for single layer free surface Navier Stokes problem
#---------------------------------------------------------------
cd Validation

echo "Running falling jet Navier Stokes validation "
mkdir RESLT_TH RESLT_CR
../falling_jet > OUTPUT_falling_jet
echo "done"
echo " " >> validation.log
echo "Falling jet Navier Stokes validation" >> validation.log
echo "--------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_TH/soln1.dat  RESLT_TH/soln2.dat  > resultsTH.dat
cat  RESLT_CR/soln1.dat  RESLT_CR/soln2.dat  > resultsCR.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "Taylor Hood Elements " >> validation.log
../../../../bin/fpdiff.py ../validata/resultsTH.dat.gz  \
         resultsTH.dat 0.1 1.0e-12 >> validation.log
echo "Crouzeix Raviart Elements " >> validation.log
../../../../bin/fpdiff.py ../validata/resultsCR.dat.gz  \
         resultsCR.dat 0.1 1.0e-12 >> validation.log
fi

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
