#!/bin/sh

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1

# Set the number of tests to be checked
NUM_TESTS=12

# Setup validation directory
#---------------------------
touch Validation
rm -rf Validation
mkdir Validation
cd Validation

# Alias for the location of the fpdiff python script. To make aliases work
# inside non-interactive shells (i.e. just in a script), we need to enable the
# "expand_aliases" shell option.
shopt -s expand_aliases
alias fpdiff="$OOMPH_ROOT_DIR/scripts/fpdiff.py"

# Validation for equal-order Galerkin discretisation
#---------------------------------------------------
mkdir RESLT

echo "Running test_equal_order_galerkin"
../test_equal_order_galerkin >OUTPUT
echo "done "
echo " " >>validation.log
echo "Equal-order Galerkin discretisation validation" >>validation.log
echo "----------------------------------------------" >>validation.log
echo " " >>validation.log
echo "Validation directory: " >>validation.log
echo " " >>validation.log
echo "  " $(pwd) >>validation.log
echo " " >>validation.log

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log
else
  validata_folder="../validata/eq_order_gal"
  fpdiff $validata_folder/soln_pbc0.dat.gz ./RESLT/soln_pbc0.dat 0.1 1.0e-8 >>validation.log
  fpdiff $validata_folder/soln_pbc1.dat.gz ./RESLT/soln_pbc1.dat 0.1 1.0e-8 >>validation.log
  fpdiff $validata_folder/soln_ic0.dat.gz ./RESLT/soln_ic0.dat 0.1 1.0e-8 >>validation.log
  fpdiff $validata_folder/soln_ic1.dat.gz ./RESLT/soln_ic1.dat 0.1 1.0e-8 >>validation.log
fi

mv RESLT RESLT_equal_order_galerkin

# Validation for equal-order Galerkin-Petrov discretisation
#----------------------------------------------------------
mkdir RESLT

echo "Running test_equal_order_galerkin_petrov"
../test_equal_order_galerkin_petrov >OUTPUT
echo "done "
echo " " >>validation.log
echo "Equal-order Galerkin-Petrov discretisation validation" >>validation.log
echo "-----------------------------------------------------" >>validation.log
echo " " >>validation.log
echo "Validation directory: " >>validation.log
echo " " >>validation.log
echo "  " $(pwd) >>validation.log
echo " " >>validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log
else
  validata_folder="../validata/eq_order_gal_pet"
  fpdiff $validata_folder/soln_pbc0.dat.gz ./RESLT/soln_pbc0.dat 0.1 1.0e-8 >>validation.log
  fpdiff $validata_folder/soln_pbc1.dat.gz ./RESLT/soln_pbc1.dat 0.1 1.0e-8 >>validation.log
  fpdiff $validata_folder/soln_ic0.dat.gz ./RESLT/soln_ic0.dat 0.1 1.0e-8 >>validation.log
  fpdiff $validata_folder/soln_ic1.dat.gz ./RESLT/soln_ic1.dat 0.1 1.0e-8 >>validation.log
fi

mv RESLT RESLT_equal_order_galerkin_petrov

# Validation for mixed-order Galerkin-Petrov discretisation
#----------------------------------------------------------
mkdir RESLT

echo "Running test_mixed_order_galerkin_petrov"
../test_mixed_order_galerkin_petrov >OUTPUT
echo "done "
echo " " >>validation.log
echo "Mixed-order Galerkin-Petrov discretisation validation" >>validation.log
echo "-----------------------------------------------------" >>validation.log
echo " " >>validation.log
echo "Validation directory: " >>validation.log
echo " " >>validation.log
echo "  " $(pwd) >>validation.log
echo " " >>validation.log

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log
else
  validata_folder="../validata/mix_order_gal_pet"
  fpdiff $validata_folder/soln_pbc0.dat.gz ./RESLT/soln_pbc0.dat 0.1 1.0e-8 >>validation.log
  fpdiff $validata_folder/soln_pbc1.dat.gz ./RESLT/soln_pbc1.dat 0.1 1.0e-8 >>validation.log
  fpdiff $validata_folder/soln_ic0.dat.gz ./RESLT/soln_ic0.dat 0.1 1.0e-8 >>validation.log
  fpdiff $validata_folder/soln_ic1.dat.gz ./RESLT/soln_ic1.dat 0.1 1.0e-8 >>validation.log
fi

mv RESLT RESLT_mixed_order_galerkin_petrov

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
