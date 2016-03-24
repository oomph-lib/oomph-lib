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

# Validation for unstructured two-d circle Darcy
#-------------------------------------------------------
cd Validation

echo "Running unstructured two-d circle Darcy validation "
mkdir RESLT
../unstructured_two_d_circle --element_area 0.01 --dir RESLT > OUTPUT_unstructured_two_d_circle
echo "done"
echo " " >> validation.log
echo "Unstructured two-d circle Darcy validation" >> validation.log
echo "--------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
touch unstructured_two_d_circle.dat
rm -f unstructured_two_d_circle.dat
head -n 2 RESLT/trace.dat | cut -d ' ' -f 2-3 >> unstructured_two_d_circle.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
 ../../../../bin/fpdiff.py ../validata/unstructured_two_d_circle.dat.gz   \
  unstructured_two_d_circle.dat  >> validation.log
fi


mv RESLT RESLT_non_adapt

# Validation for adaptive unstructured two-d circle Darcy
#-------------------------------------------------------

echo "Running adaptive unstructured two-d circle Darcy validation "
mkdir RESLT
../adaptive_unstructured_two_d_circle --element_area 0.01 --dir RESLT > OUTPUT_adaptive_unstructured_two_d_circle
echo "done"
echo " " >> validation.log
echo "Adaptive unstructured two-d circle Darcy validation" >> validation.log
echo "---------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
touch adaptive_unstructured_two_d_circle.dat
rm -f adaptive_unstructured_two_d_circle.dat
head -n 2 RESLT/trace.dat | cut -d ' ' -f 2-3 >> adaptive_unstructured_two_d_circle.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
 ../../../../bin/fpdiff.py ../validata/adaptive_unstructured_two_d_circle.dat.gz   \
  adaptive_unstructured_two_d_circle.dat 4.0 1.0e-14 >> validation.log
fi


mv RESLT RESLT_adapt




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
