#! /bin/sh

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1


#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for axisymmetric single layer free surface Navier-Stokes problem
#----------------------------------------------------------------------------
cd Validation

echo "Running axisymmetric single layer free surface Navier Stokes validation "
mkdir RESLT
../single_layer_free_surface_axisym lalala > OUTPUT_single_layer_free_surface_axisym
echo "done"
echo " " >> validation.log
echo "Axisymmetric single layer free surface Navier-Stokes validation" >> validation.log
echo "---------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0.dat  RESLT/soln1.dat RESLT/soln2.dat > single_layer_free_surface_axisym_results.dat

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/results.dat.gz  \
         single_layer_free_surface_axisym_results.dat >> validation.log
fi

# Append log to main validation log
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
