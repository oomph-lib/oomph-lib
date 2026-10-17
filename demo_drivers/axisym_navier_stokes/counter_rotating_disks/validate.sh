#! /bin/sh

# Get the OOMPH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$1


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
mkdir RESLT_TH
mkdir RESLT_CR

../counter_rotating_disks lala > OUTPUT
echo "done"
echo " " >> validation.log
echo "Counter-rotating disks validation (non-refineable)" >> validation.log
echo "--------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT_TH/base_soln0.dat RESLT_TH/perturbed_soln0.dat \
     RESLT_TH/trace.dat \
 > results_TH_nonref.dat
cat  RESLT_CR/base_soln0.dat RESLT_CR/perturbed_soln0.dat \
     RESLT_CR/trace.dat \
 > results_CR_nonref.dat

mv RESLT_TH RESLT_TH_nonref
mv RESLT_CR RESLT_CR_nonref

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/results_TH_nonref.dat.gz  \
         results_TH_nonref.dat 0.1 1.0e-8 >> validation.log
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/results_CR_nonref.dat.gz  \
         results_CR_nonref.dat 0.1 1.0e-8 >> validation.log
fi


echo "Running refineable counter-rotating disks validation "
mkdir RESLT_TH
mkdir RESLT_CR

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
cat  RESLT_TH/base_soln0.dat RESLT_TH/perturbed_soln0.dat \
     RESLT_TH/trace.dat \
 > results_TH_ref.dat
cat  RESLT_CR/base_soln0.dat RESLT_CR/perturbed_soln0.dat \
     RESLT_CR/trace.dat \
 > results_CR_ref.dat

mv RESLT_TH RESLT_TH_ref
mv RESLT_CR RESLT_CR_ref

if test "$2" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/results_TH_ref.dat.gz  \
         results_TH_ref.dat 0.1 1.0e-11  >> validation.log
$OOMPH_ROOT_DIR/scripts/fpdiff.py ../validata/results_CR_ref.dat.gz  \
         results_CR_ref.dat 0.1 1.0e-11 >> validation.log
fi



# Append log to main validation log
cat validation.log >> $OOMPH_ROOT_DIR/validation.log

cd ..



#######################################################################


#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/scripts/validate_ok_count

# Never get here
exit 10
