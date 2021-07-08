#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for acoustic fsi problem
#------------------------------------
cd Validation

echo "Running acoustic fsi validation "
mkdir RESLT
../acoustic_fsi > OUTPUT_structured
echo "done"
echo " " >> validation.log
echo "Acoustic fsi validation" >> validation.log
echo "-----------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/helmholtz_soln1.dat RESLT/elast_soln1.dat > result.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result.dat.gz  \
         result.dat 0.1 1.0e-9 >> validation.log
fi

mv RESLT RESLT_structured

# Validation for acoustic fsi problem
#------------------------------------

echo "Running unstructured acoustic fsi validation "
mkdir RESLT
../unstructured_acoustic_fsi --validation > OUTPUT_unstructured
echo "done"
echo " " >> validation.log
echo "Unstructured acoustic fsi validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > result_unstructured.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_unstructured.dat.gz  \
         result_unstructured.dat  0.1 1.0e-9 >> validation.log
fi

mv RESLT RESLT_unstructured

# Append log to main validation log
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
