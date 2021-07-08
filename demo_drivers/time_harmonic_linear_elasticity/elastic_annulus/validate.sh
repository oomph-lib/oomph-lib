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

# Validation 
#-----------
cd Validation

echo "Running adaptive time-periodic linear elasticity for annulus"
mkdir RESLT
../adaptive_time_harmonic_elastic_annulus > OUTPUT_adapt
echo "done"
echo " " >> validation.log
echo "Adaptive time-periodic linear elasticity for annulus" >> validation.log
echo "-----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/elast_soln2.dat > adaptive_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_results.dat.gz   \
    adaptive_results.dat 0.1 1.0e-10  >> validation.log
fi
mv RESLT RESLT_adapt



echo "Running non-adaptive time-periodic linear elasticity for annulus"
mkdir RESLT
../time_harmonic_elastic_annulus > OUTPUT_non_adapt
echo "done"
echo " " >> validation.log
echo "Non-adaptive time-periodic linear elasticity for annulus" >> validation.log
echo "--------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/elast_soln2.dat > non_adaptive_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/non_adaptive_results.dat.gz   \
    non_adaptive_results.dat 0.1 1.0e-10 >> validation.log
fi
mv RESLT RESLT_non_adapt




echo "Running adaptive unstructured time-periodic linear elasticity for annulus"
mkdir RESLT
../adaptive_unstructured_time_harmonic_elastic_annulus --validation --alpha 0.0 > OUTPUT_unstructured_adapt
echo "done"
echo " " >> validation.log
echo "Adaptive unstructured time-periodic linear elasticity for annulus" >> validation.log
echo "------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/norm0.dat RESLT/norm1.dat > unstructured_adaptive_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_adaptive_results.dat.gz   \
    unstructured_adaptive_results.dat  4.0 1.0e-12 >> validation.log
fi
mv RESLT RESLT_unstructured_adapt



echo "Running non-adaptive unstructured time-periodic linear elasticity for annulus"
mkdir RESLT
../unstructured_time_harmonic_elastic_annulus --validation --alpha 0.0 > OUTPUT_unstructured_non_adapt
echo "done"
echo " " >> validation.log
echo "Non-adaptive unstructured time-periodic linear elasticity for annulus" >> validation.log
echo "---------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/norm0.dat RESLT/norm1.dat > unstructured_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/unstructured_results.dat.gz   \
    unstructured_results.dat 0.1 1.0e-10 >> validation.log
fi
mv RESLT RESLT_unstructured









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
