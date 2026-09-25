#! /bin/sh

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1

#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo axisym solid
#----------------------------------
cd Validation

# Non-refineable version
echo "Running non-adaptive axisym solid validation"
mkdir RESLT
../shell >OUTPUT_axisym_solid
echo "done"
echo " " >>validation.log
echo "Axisym solid demo validation (adaptive and non-adaptive)" >>validation.log
echo "--------------------------------------------------------" >>validation.log
echo " " >>validation.log
echo "Validation directory: " >>validation.log
echo " " >>validation.log
echo "  " $(pwd) >>validation.log
echo " " >>validation.log
cat RESLT/bulk_soln_nu_0.4_pres_0.05.dat > shell_results.dat

# Refineable version
echo "Running adaptive axisym solid validation"
../shell_adaptive >OUTPUT_axisym_solid
echo "done"
echo " " >>validation.log
cat RESLT/bulk_soln_nu_0.4_pres_0.05.dat > shell_adaptive_results.dat

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >>validation.log
else
  $OOMPH_ROOT_DIR/scripts/fpdiff.py \
    ../validata/shell_results.dat.gz \
    shell_results.dat >>validation.log

  $OOMPH_ROOT_DIR/scripts/fpdiff.py \
    ../validata/shell_adaptive_results.dat.gz \
    shell_adaptive_results.dat >>validation.log
fi

# Append output to global validation log file
#--------------------------------------------
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
exit 0
