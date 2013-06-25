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

echo "Running three-dimensional free surface Navier Stokes validation "
mkdir RESLT RESLT_elastic
../3d_free_surface lalala > OUTPUT_3d_free_surface
echo "done"
echo " " >> validation.log
echo "3D free surface Navier Stokes validation" >> validation.log
echo "--------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln2.dat  RESLT/trace.dat > results_spine.dat
cat  RESLT_elastic/soln2.dat  RESLT_elastic/trace.dat > results_elastic.dat

echo " Spine representation " >> validation.log
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_spine.dat.gz  \
         results_spine.dat 0.1 1.0e-12 >> validation.log
fi

echo " Pseudo-elastic representation " >> validation.log
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_elastic.dat.gz  \
         results_elastic.dat 0.1 1.0e-12 >> validation.log
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
