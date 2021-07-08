#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for Womersley problems
#----------------------------------
cd Validation

echo "Running Navier-Stokes flux control validation"
mkdir RESLT_flux_control
../flux_control -validation_run -outflow 0 > OUTPUT
echo "done"
echo " " >> validation.log
echo "Navier-Stokes flux control validation" >> validation.log
echo "-------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log



echo "Navier Stokes flux control test: " >> validation.log
cat  RESLT_flux_control/soln0.dat RESLT_flux_control/soln1.dat RESLT_flux_control/soln9.dat RESLT_flux_control/soln10.dat RESLT_flux_control/trace.dat \
> flux_control.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py \
../validata/flux_control.dat.gz  \
flux_control.dat >> validation.log
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
