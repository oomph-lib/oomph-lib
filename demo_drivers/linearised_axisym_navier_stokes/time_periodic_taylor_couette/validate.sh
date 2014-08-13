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

# Validation for time-periodic Taylor-Couette
#--------------------------------------------
cd Validation

echo "Running time-periodic Taylor-Couette validation "
mkdir RESLT_CR
mkdir RESLT_TH
../time_periodic_taylor_couette lala > OUTPUT
echo "done"
echo " " >> validation.log
echo "Time-periodic Taylor-Couette validation" >> validation.log
echo "---------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_CR/base_soln_epsilon0.9_Re20.93_soln4.dat \
    RESLT_CR/base_soln_epsilon0.9_Re20.93_soln9.dat \
    RESLT_CR/base_soln_epsilon0.9_Re20.93_soln14.dat \
    RESLT_CR/base_soln_epsilon0.9_Re20.93_soln19.dat \
    RESLT_CR/perturbed_soln_epsilon0.9_Re20.93_soln4.dat \
    RESLT_CR/perturbed_soln_epsilon0.9_Re20.93_soln9.dat \
    RESLT_CR/perturbed_soln_epsilon0.9_Re20.93_soln14.dat \
    RESLT_CR/perturbed_soln_epsilon0.9_Re20.93_soln19.dat \
    RESLT_CR/power_method_trace_epsilon0.9_Re20.93.dat > results_CR.dat
cat RESLT_TH/base_soln_epsilon0.9_Re20.93_soln4.dat \
    RESLT_TH/base_soln_epsilon0.9_Re20.93_soln9.dat \
    RESLT_TH/base_soln_epsilon0.9_Re20.93_soln14.dat \
    RESLT_TH/base_soln_epsilon0.9_Re20.93_soln19.dat \
    RESLT_TH/perturbed_soln_epsilon0.9_Re20.93_soln4.dat \
    RESLT_TH/perturbed_soln_epsilon0.9_Re20.93_soln9.dat \
    RESLT_TH/perturbed_soln_epsilon0.9_Re20.93_soln14.dat \
    RESLT_TH/perturbed_soln_epsilon0.9_Re20.93_soln19.dat \
    RESLT_TH/power_method_trace_epsilon0.9_Re20.93.dat > results_TH.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_CR.dat.gz  \
         results_CR.dat 0.1 4.0e-10 >> validation.log
../../../../bin/fpdiff.py ../validata/results_TH.dat.gz  \
         results_TH.dat 0.1 1.0e-9 >> validation.log
fi

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
