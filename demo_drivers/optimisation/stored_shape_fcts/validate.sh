#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1


#Set the number of tests to be checked
NUM_TESTS=1


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation



# Validation for demo poisson with stored shape functions
#--------------------------------------------------------
cd Validation
echo "Running 2D poisson with stored shape functions demo "
mkdir RESLT
echo "S" >> tmp_input.dat
../two_d_poisson_stored_shape_fcts 5 < tmp_input.dat > OUTPUT_two_d_poisson_stored_shape_fcts
rm tmp_input.dat
echo "done"
echo " " >> validation.log
echo "2D poisson with stored shape fcts demo" >> validation.log
echo "--------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat  \
    > results.dat

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "Compare results: " >> validation.log
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/results.dat.gz   \
    results.dat  >> validation.log
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
