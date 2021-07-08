#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=4


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for oscillating sphere
#----------------------------------
cd Validation

echo "Running adaptive oscillating sphere validation (N=0). Perfectly matched layers"
mkdir RESLT
../adaptive_oscillating_sphere --validate > OUTPUT_N0_adapt
echo "done"
echo " " >> validation.log
echo "Adaptive oscillating sphere validation N=0 (Perfectly matched layers)" >> validation.log
echo "---------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > adaptive_results_n0.dat
mv RESLT RESLT_N0_adapt

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_results_n0.dat.gz   \
    adaptive_results_n0.dat  >> validation.log
fi

echo "Running adaptive oscillating sphere validation (N=3). Perfectly matched layers"
mkdir RESLT
../adaptive_oscillating_sphere --validate --Fourier_wavenumber 3 > OUTPUT_N3_adapt
echo "done"
echo " " >> validation.log
echo "Adaptive oscillating sphere validation N=3 (Perfectly matched layers)" >> validation.log
echo "---------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > adaptive_results_n3.dat
mv RESLT RESLT_N3_adapt


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_results_n3.dat.gz   \
    adaptive_results_n3.dat  >> validation.log
fi

echo "Running oscillating sphere validation (N=0). Perfectly matched layers"
mkdir RESLT
../oscillating_sphere --validate > OUTPUT_N0
echo "done"
echo " " >> validation.log
echo "Oscillating sphere validation N=0 (Perfectly matched layers)" >> validation.log
echo "-----------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > results_n0.dat
mv RESLT RESLT_N0

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_n0.dat.gz   \
    results_n0.dat  >> validation.log
fi

echo "Running oscillating sphere validation (N=3). Perfectly matched layers"
mkdir RESLT
../oscillating_sphere --validate --Fourier_wavenumber 3 > OUTPUT_N3
echo "done"
echo " " >> validation.log
echo "Oscillating sphere validation N=3 (Perfectly matched layers)" >> validation.log
echo "------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > results_n3.dat
mv RESLT RESLT_N3


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_n3.dat.gz   \
    results_n3.dat 0.25 1.0e-14 >> validation.log
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
