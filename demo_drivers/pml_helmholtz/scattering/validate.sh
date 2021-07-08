#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1


#Set the number of tests to be checked
NUM_TESTS=2


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for non-adaptive scattering
#-----------------------------------
cd Validation

echo "Running unstructured PML scattering"
mkdir RESLT
../unstructured_two_d_helmholtz > OUTPUT_non_adapt
echo "done"
echo " " >> validation.log
echo "Running unstructured PML scattering" >> validation.log
echo "-----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > pml_scattering_results.dat

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/pml_scattering_results.dat.gz   \
    pml_scattering_results.dat  2.0 1.0e-14 >> validation.log
fi
mv RESLT RESLT_non_adapt


# Validation for adaptive scattering
#-----------------------------------

echo "Running adaptive unstructured PML scattering"
mkdir RESLT
../unstructured_two_d_helmholtz_adapt > OUTPUT_adapt
echo "done"
echo " " >> validation.log
echo "Running adaptive unstructured PML scattering" >> validation.log
echo "--------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > pml_adaptive_scattering_results.dat

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/pml_scattering_results.dat.gz   \
    pml_scattering_results.dat  2.0 1.0e-14 >> validation.log
fi
mv RESLT RESLT_adapt




# Append output to global validation log file
#--------------------------------------------
cat validation.log >> $OOMPH_ROOT_DIR/validation.log


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
