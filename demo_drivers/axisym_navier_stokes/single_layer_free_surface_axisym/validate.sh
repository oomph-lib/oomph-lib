#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


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

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz  \
         single_layer_free_surface_axisym_results.dat >> validation.log
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
