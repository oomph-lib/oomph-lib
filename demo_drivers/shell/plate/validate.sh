#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

cd Validation



# Validation for flat plate deforming like a beam
#------------------------------------------------

mkdir RESLT

echo "Running square flat plate "
../plate > OUTPUT

echo "done"
echo " " >> validation.log
echo "Clamped shell" \
 >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat RESLT/plate0.dat RESLT/plate1.dat  > plate.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/plate.dat.gz \
 plate.dat>> validation.log
fi



echo "Running unstructured square flat plate "
mkdir RESLT_unstructured_plate
../unstructured_clamped_square_plate ../UnitPlate.1.node ../UnitPlate.1.ele ../UnitPlate.1.poly > OUTPUT_unstructured_plate

echo "done"
echo " " >> validation.log
echo "Unstructured Clamped square flat plate" \
 >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_unstructured_plate/soln1.dat RESLT_unstructured_plate/soln2.dat > unstructured_square_plate.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py  ../validata/unstructured_square_plate.dat.gz \
 unstructured_square_plate.dat 0.3 1.0e-14 >> validation.log
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
