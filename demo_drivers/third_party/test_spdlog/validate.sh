#!/bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1

#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo poisson
#----------------------------
cd Validation

echo "Running nlohmann/json validation "
mkdir RESLT
cd RESLT
../../test_spdlog >../OUTPUT_test_spdlog || {
  echo 'Failed to run ./test_spdlog!'
  exit 1
}
cd ..
echo "done"
echo " " >>validation.log
echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log

# Append output to global validation log file
#--------------------------------------------
cat validation.log >>$OOMPH_ROOT_DIR/validation.log || {
  echo 'Failed to append validation.log to global validation.log!'
  exit 1
}

cd ..

#######################################################################

#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/scripts/validate_ok_count

# Never get here
exit 0
