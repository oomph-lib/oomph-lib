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

# Validation for static fish deformation
#---------------------------------------

cd Validation
mkdir RESLT

echo "Running unstructured adaptive 2D FSI validation "
../unstructured_adaptive_2d_fsi --validation > OUTPUT_adaptive_2d_fsi


echo "done"
echo " " >> validation.log
echo "Unstructured adaptive 2D FSI validation" >> validation.log
echo "----------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > uns_adapt.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/uns_adapt.dat.gz \
    uns_adapt.dat 4.0 1.0e-12 >> validation.log
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
