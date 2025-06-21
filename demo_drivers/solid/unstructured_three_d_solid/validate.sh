#! /bin/sh

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1

#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

cd Validation


# Validation for 3D unstructured solid
#-------------------------------------

# Get tetgen files
cp ../fsi_bifurcation_solid.* .

mkdir RESLT

echo "Running 3D unstructured solid "
../unstructured_three_d_solid la > OUTPUT

echo "done"
echo " " >> validation.log
echo "3D unstructured solid" >> validation.log
echo "---------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/solid_soln0.dat RESLT/solid_soln1.dat RESLT/solid_soln2.dat > result_three_d.dat

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/result_three_d.dat.gz \
    result_three_d.dat  >> validation.log
fi



# Append output to global validation log file
#--------------------------------------------
cat validation.log >> $OOMPH_ROOT_DIR/validation.log




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
