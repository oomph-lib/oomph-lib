#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=5


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation 
#-----------
cd Validation

echo "Running Fourier decomposed time harmonic lin elast validation "
mkdir RESLT
../cylinder > OUTPUT_structured
echo "done"
echo " " >> validation.log
echo "Fourier decomposed time harmonic lin elast validation" >> validation.log
echo "-----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln.dat  > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz   \
    results.dat  >> validation.log
fi

mv RESLT RESLT_structured


echo "Running unstructured Fourier decomposed time harmonic lin elast validation "
mkdir RESLT
../unstructured_cylinder > OUTPUT_unstructured
echo "done"
echo " " >> validation.log
echo "Unstructured Fourier decomposed time harmonic lin elast validation" >> validation.log
echo "------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/norm.dat  > results_unstructured.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_unstructured.dat.gz   \
    results_unstructured.dat  >> validation.log
fi

mv RESLT RESLT_unstructured



echo "Running unstructured Fourier decomposed time harmonic lin elast validation "
mkdir RESLT
../adaptive_unstructured_cylinder > OUTPUT_adaptive_unstructured
echo "done"
echo " " >> validation.log
echo "Adaptive unstructured Fourier decomposed time harmonic lin elast validation" >> validation.log
echo "---------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/norm.dat  > results_adaptive_unstructured.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_adaptive_unstructured.dat.gz   \
    results_adaptive_unstructured.dat  >> validation.log
fi

mv RESLT RESLT_adaptive_unstructured



echo "Running Fourier decomposed time harmonic lin elast (press load)  validation "
mkdir RESLT
../pressure_loaded_cylinder > OUTPUT_pressure_loaded_unstructured
echo "done"
echo " " >> validation.log
echo "Fourier decomposed time harmonic lin elast (press load) validation" >> validation.log
echo "------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln.dat  > results_press_load.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_press_load.dat.gz   \
    results_press_load.dat  0.1 1.0e-12 >> validation.log
fi

mv RESLT RESLT_press_load



echo "Running unstructured Fourier decomposed time harmonic lin elast (press load) validation "
mkdir RESLT
../adaptive_pressure_loaded_cylinder > OUTPUT_press_load_adaptive_unstructured
echo "done"
echo " " >> validation.log
echo "Adaptive unstructured Fourier decomposed time harmonic lin elast (press load) validation" >> validation.log
echo "----------------------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/norm.dat  > results_press_load_adaptive_unstructured.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_press_load_adaptive_unstructured.dat.gz   \
    results_press_load_adaptive_unstructured.dat 4.0 1.0e-14  >> validation.log
fi

mv RESLT RESLT_press_load_adaptive_unstructured



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
