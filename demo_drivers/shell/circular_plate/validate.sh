#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

cd Validation



# Validation for flat plate deforming like a beam
#------------------------------------------------

mkdir RESLT

echo "Running circular flat plate "
../unstructured_clamped_circular_plate ../CircularPlate.1.node ../CircularPlate.1.ele ../CircularPlate.1.poly > OUTPUT

echo "done"
echo " " >> validation.log
echo "Unstructured circular flat plate" \
 >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln1.dat RESLT/soln2.dat > circular_plate.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/circular_plate.dat.gz \
 circular_plate.dat 5.0 4.0e-5 >> validation.log
fi



# Append output to global validation log file
#--------------------------------------------
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
