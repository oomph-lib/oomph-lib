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

# Validation for Elastic bretherton problem
#-------------------------------------
cd Validation

echo "Running elastic Bretherton problem validation "
mkdir RESLT
../elastic_bretherton lalala > OUTPUT
echo "done"
echo " " >> validation.log
echo "Elastic Bretherton validation" >> validation.log
echo "------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/trace.dat  RESLT/soln2.dat > el_breth.dat

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/el_breth.dat.gz  \
         el_breth.dat  0.5 9.0e-9 >> validation.log
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
