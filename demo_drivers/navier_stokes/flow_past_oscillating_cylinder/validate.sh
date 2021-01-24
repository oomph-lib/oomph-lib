#!/bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

# Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
touch Validation
rm -rf Validation
mkdir Validation

# Alias for the location of the fpdiff python script; note that the oomph
# root directory is one level higher because we descended into the Validation
# directory after getting the relative address for the root directory
alias fpdiff="../$OOMPH_ROOT_DIR/bin/fpdiff.py"

# Validation for flow past cylinder
#-----------------------------------------
cd Validation
mkdir RESLT

echo "Running flow_past_oscillating_cylinder"
../flow_past_oscillating_cylinder > OUTPUT
echo "done "
echo " " >> validation.log
echo "Flow past an oscillating cylinder time-integration validation" >> validation.log
echo "-------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  fpdiff ../validata/soln0.dat.gz ./RESLT/soln0.dat 0.1 1.0e-8 >> validation.log
  fpdiff ../validata/unsteady0.dat.gz ./RESLT/unsteady0.dat 0.1 1.0e-8 >> validation.log
fi

# Append log to main validation log
cat validation.log >> ../$OOMPH_ROOT_DIR/validation.log

cd ..

#######################################################################

# Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
