#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1

#Set the number of tests 
NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation


# Validation for channel with leaflet
#------------------------------------
cd Validation

echo "Running validation for channel with leaflet  "
rm -rf RESLT
mkdir RESLT

# Do validation run
../channel_with_leaflet blabla > OUTPUT 
echo "done"
echo " " >> validation.log

cat RESLT/soln0.dat RESLT/soln3.dat > result.dat

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/result.dat.gz \
    result.dat >> validation.log
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
