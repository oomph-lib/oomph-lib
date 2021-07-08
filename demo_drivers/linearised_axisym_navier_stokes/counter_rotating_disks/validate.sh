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

# Validation for counter-rotating disks
#--------------------------------------
cd Validation

echo "Running non-refineable counter-rotating disks validation "
mkdir RESLT_CR
mkdir RESLT_TH
../counter_rotating_disks lala > OUTPUT_nonref
echo "done"
echo " " >> validation.log
echo "Counter-rotating disks validation (non-refineable)" >> validation.log
echo "--------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_CR/base_soln_k2_Re300.00_soln0.dat \
    RESLT_CR/perturbed_soln_k2_Re300.00_soln0.dat \
    RESLT_CR/power_method_trace_k2_Re300.00.dat \
 > results_CR.dat
cat RESLT_TH/base_soln_k2_Re300.00_soln0.dat \
    RESLT_TH/perturbed_soln_k2_Re300.00_soln0.dat \
    RESLT_TH/power_method_trace_k2_Re300.00.dat \
 > results_TH.dat

mv RESLT_TH RESLT_TH_nonref
mv RESLT_CR RESLT_CR_nonref

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_CR_nonref.dat.gz  \
         results_CR.dat  0.1 1.5e-9 >> validation.log
../../../../bin/fpdiff.py ../validata/results_TH_nonref.dat.gz  \
         results_TH.dat 0.1 1.0e-11 >> validation.log
fi

rm results_TH.dat
rm results_CR.dat

echo "Running refineable counter-rotating disks validation "
mkdir RESLT_CR
mkdir RESLT_TH
../counter_rotating_disks_ref lala > OUTPUT_ref
echo "done"
echo " " >> validation.log
echo "Counter-rotating disks validation (refineable)" >> validation.log
echo "----------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_CR/base_soln_k2_Re300.00_soln0.dat \
    RESLT_CR/perturbed_soln_k2_Re300.00_soln0.dat \
    RESLT_CR/power_method_trace_k2_Re300.00.dat \
 > results_CR.dat
cat RESLT_TH/base_soln_k2_Re300.00_soln0.dat \
    RESLT_TH/perturbed_soln_k2_Re300.00_soln0.dat \
    RESLT_TH/power_method_trace_k2_Re300.00.dat \
 > results_TH.dat

mv RESLT_TH RESLT_TH_ref
mv RESLT_CR RESLT_CR_ref

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_CR_ref.dat.gz  \
         results_CR.dat 0.1 1.5e-9 >> validation.log
../../../../bin/fpdiff.py ../validata/results_TH_ref.dat.gz  \
         results_TH.dat 0.1 1.0e-13 >> validation.log
fi

rm results_TH.dat
rm results_CR.dat

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
