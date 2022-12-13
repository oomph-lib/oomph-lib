#! /bin/sh

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1

#Set the number of tests to be checked
NUM_TESTS=1

# Append output to global validation log file
#--------------------------------------------
touch Validation
rm -r -f Validation
mkdir Validation
echo "[OK]" >> Validation/validation.log
cat Validation/validation.log >> $OOMPH_ROOT_DIR/validation.log



#######################################################################


#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/scripts/validate_ok_count

