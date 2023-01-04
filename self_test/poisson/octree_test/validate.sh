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


cd Validation

echo "Running octree validation "
mkdir RESLT
../octree_test lala > OUTPUT
echo "done"
echo " " >> validation.log
echo "Octree validation" >> validation.log
echo "-----------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
cat trace.dat  RESLT/edge_neighbours0.dat RESLT/neighbours0.dat \
    RESLT/no_true_edge0.dat RESLT/orientation0.dat > results.dat
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/results.dat.gz results.dat >> validation.log
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
