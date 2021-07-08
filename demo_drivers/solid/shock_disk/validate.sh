#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1

#Set the number of tests to be checked
NUM_TESTS=3

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

# Validation for shock disk
#--------------------------

cd Validation
mkdir RESLT0
mkdir RESLT1
mkdir RESLT2

echo "Running shock disk validation "
../shock_disk lalala  > OUTPUT


echo "done"
echo " " >> validation.log
echo "Shock disk validation" >> validation.log
echo "---------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT0/soln3.dat \
    > shock_disk_results0.dat
cat RESLT1/soln3.dat \
    > shock_disk_results1.dat
cat RESLT2/soln3.dat \
    > shock_disk_results2.dat

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/shock_disk_results0.dat.gz \
    shock_disk_results0.dat  >> validation.log
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/shock_disk_results1.dat.gz \
    shock_disk_results1.dat  >> validation.log
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/shock_disk_results2.dat.gz \
    shock_disk_results2.dat  >> validation.log
fi

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
