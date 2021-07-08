#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=3

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for acoustic fsi problem
#------------------------------------
cd Validation

echo "Running Fourier-decomposed acoustic fsi validation "
mkdir RESLT
../fourier_decomposed_acoustic_fsi --nstep 2 > OUTPUT_structured
echo "done"
echo " " >> validation.log
echo "Fourier-decomposed acoustic fsi validation" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/helmholtz_soln1.dat RESLT/elast_soln1.dat RESLT/trace.dat > result.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result.dat.gz  \
         result.dat >> validation.log
fi

mv RESLT RESLT_structured



# Validation for unstructured acoustic fsi problem
#-------------------------------------------------

echo "Running unstructured Fourier-decomposed acoustic fsi validation "
mkdir RESLT
../unstructured_fourier_decomposed_acoustic_fsi > OUTPUT_unstructured
echo "done"
echo " " >> validation.log
echo "Unstructured Fourier-decomposed acoustic fsi validation" >> validation.log
echo "-------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/trace.dat > unstructured_result.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_result.dat.gz  \
         unstructured_result.dat >> validation.log
fi

mv RESLT RESLT_unstructured


# Validation for adaptive unstructured acoustic fsi problem
#-----------------------------------------------------------

echo "Running adaptive unstructured Fourier-decomposed acoustic fsi validation "
mkdir RESLT
../adaptive_unstructured_fourier_decomposed_acoustic_fsi > OUTPUT_adaptive_unstructured
echo "done"
echo " " >> validation.log
echo "Unstructured adaptive Fourier-decomposed acoustic fsi validation" >> validation.log
echo "----------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/trace.dat > adaptive_unstructured_result.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_unstructured_result.dat.gz  \
         adaptive_unstructured_result.dat 1.0e-14 0.25 >> validation.log
fi

mv RESLT RESLT_adaptive_unstructured


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
