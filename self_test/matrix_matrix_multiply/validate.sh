#! /bin/sh

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1


#Set the number of tests to be checked
NUM_TESTS=7

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for circular driven cavity
#-----------------------------------------
cd Validation

echo "Running matrix multiply validation "
mkdir RESLT
../matrix_multiply_test> OUTPUT
echo "done"
echo " " >> validation.log
echo "Matrix multiply validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log 
sort -k 1 -k 2 -n  CC_result1.dat > CC_result1.dat.sorted
sort -k 1 -k 2 -n  CC_result2.dat > CC_result2.dat.sorted
sort -k 1 -k 2 -n  CC_result3.dat > CC_result3.dat.sorted
sort -k 1 -k 2 -n  CR_result1.dat > CR_result1.dat.sorted
sort -k 1 -k 2 -n  CR_result2.dat > CR_result2.dat.sorted
sort -k 1 -k 2 -n  CR_result3.dat > CR_result3.dat.sorted
sort -k 1 -k 2 -n  D_result.dat   > D_result.dat.sorted




if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "CC method 1: ">> validation.log
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/result.dat.gz CC_result1.dat.sorted 0.1 2.0e-14 >> validation.log
echo "CC method 2: ">> validation.log
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/result.dat.gz   CC_result2.dat.sorted  0.1 2.0e-14 >> validation.log
echo "CC method 3: ">> validation.log
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/result.dat.gz   CC_result3.dat.sorted  0.1 2.0e-14 >> validation.log
echo "CR method 1: ">> validation.log
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/result.dat.gz   CR_result1.dat.sorted  0.1 2.0e-14 >> validation.log
echo "CR method 2: ">> validation.log
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/result.dat.gz   CR_result2.dat.sorted  0.1 2.0e-14 >> validation.log
echo "CR method 3: ">> validation.log
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/result.dat.gz   CR_result3.dat.sorted  0.1 2.0e-14 >> validation.log
echo "Dense: ">> validation.log
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/result.dat.gz   D_result.dat.sorted  0.1 2.0e-14 >> validation.log
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
