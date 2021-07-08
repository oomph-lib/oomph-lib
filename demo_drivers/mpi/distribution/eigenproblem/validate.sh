#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

cd Validation



# Validation for distributed eigenproblem detection
#---------------------------------------------------------------------

if [ -f ../harmonic ]; then

mkdir RESLT_harmonic
cd RESLT_harmonic

echo "Running harmonic eigenproblem (parallel)  "
$MPI_RUN_COMMAND ../../harmonic > ../OUTPUT_harmonic

echo "done"
cd ..
echo " " >> validation.log
echo "Harmonic eigenproblem (parallel)  validation" \
 >> validation.log
echo "---------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_harmonic/eigenvalues1_on_proc0.dat RESLT_harmonic/eigenvalues1_on_proc1.dat RESLT_harmonic/soln1_on_proc0.dat RESLT_harmonic/soln1_on_proc1.dat > harmonic_non_dist.dat
cat RESLT_harmonic/eigenvalues2_on_proc0.dat RESLT_harmonic/eigenvalues2_on_proc1.dat RESLT_harmonic/soln2_on_proc0.dat RESLT_harmonic/soln2_on_proc1.dat > harmonic_dist.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "Non distributed" >> validation.log
../../../../../bin/fpdiff.py ../validata/harmonic_non_dist.dat.gz \
 harmonic_non_dist.dat 0.1 1.0e-13 >> validation.log
echo "Distributed" >> validation.log
../../../../../bin/fpdiff.py ../validata/harmonic_dist.dat.gz \
 harmonic_dist.dat 0.1 1.0e-13 >> validation.log
fi

else
echo "Not runnning test because Trilinos is required"
echo "[OK] (Dummy for non-existent Trilinos)" >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)" >> validation.log
fi


# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../../validation.log




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
