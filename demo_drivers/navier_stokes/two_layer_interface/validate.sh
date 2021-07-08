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

# Validation for spine two layer interface Navier Stokes problem
#---------------------------------------------------------------
cd Validation

echo "Running spine two layer interface Navier Stokes validation "
mkdir RESLT
../spine_two_layer_interface dummy_input > OUTPUT_spine_two_layer_interface
echo "done"
echo " " >> validation.log
echo "Spine two layer interface Navier--Stokes validation" >> validation.log
echo "---------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0.dat  RESLT/soln1.dat RESLT/soln2.dat > results_spine_two_layer_interface.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_spine_two_layer_interface.dat.gz results_spine_two_layer_interface.dat 0.1 2.0e-11 >> validation.log
fi

mv RESLT RESLT_spine_two_layer_interface



# Validation for elastic two layer interface Navier Stokes problem
#---------------------------------------------------------------

echo "Running elastic two layer interface Navier Stokes validation "
mkdir RESLT
../elastic_two_layer_interface --validation > OUTPUT_elastic_two_layer_interface
echo "done"
echo " " >> validation.log
echo "Elastic two layer interface Navier--Stokes validation" >> validation.log
echo "-----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0.dat  RESLT/soln1.dat RESLT/soln2.dat > results_elastic_two_layer_interface.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_elastic_two_layer_interface.dat.gz results_elastic_two_layer_interface.dat 0.1 3.0e-9 >> validation.log
fi

mv RESLT RESLT_elastic_two_layer_interface



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
