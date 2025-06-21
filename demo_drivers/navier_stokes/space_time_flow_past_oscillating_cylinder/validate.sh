#!/bin/sh

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1

# Set the number of tests to be checked
NUM_TESTS=3

# Setup validation directory
#---------------------------
touch Validation
rm -rf Validation
mkdir Validation

# Alias for the location of the fpdiff python script. To make aliases work
# inside non-interactive shells (i.e. just in a script), we need to enable the
# "expand_aliases" shell option.
shopt -s expand_aliases
alias fpdiff="$OOMPH_ROOT_DIR/scripts/fpdiff.py"

# Validation for flow past cylinder
#-----------------------------------------
cd Validation
mkdir RESLT

echo "Running space_time_oscillating_cylinder"
../space_time_oscillating_cylinder >OUTPUT
echo "done "
echo " " >>validation.log
echo "Space-time flow past oscillating cylinder validation" >>validation.log
echo "----------------------------------------------------" >>validation.log
echo " " >>validation.log
echo "Validation directory: " >>validation.log
echo " " >>validation.log
echo "  " $(pwd) >>validation.log
echo " " >>validation.log

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log
else
  validata_folder="../validata/"
  fpdiff ../validata/soln0.dat.gz ./RESLT/soln0.dat 0.1 1.0e-8 >>validation.log
  fpdiff ../validata/soln1.dat.gz ./RESLT/soln1.dat 0.1 1.0e-8 >>validation.log
  fpdiff ../validata/soln2.dat.gz ./RESLT/soln2.dat 0.1 1.0e-8 >>validation.log
fi

# Append log to main validation log
cat validation.log >>$OOMPH_ROOT_DIR/validation.log

cd ..

#######################################################################

#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/scripts/validate_ok_count

# Never get here
exit 10
