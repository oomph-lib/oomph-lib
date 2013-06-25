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

# Validation for Womersley problems
#----------------------------------
cd Validation

echo "Running 1D impedance tube validation"
mkdir RESLT_impedance_tube
mkdir RESLT_impedance_tube_with_flux_control
../one_d_womersley -validation_run -outflow 1 > OUTPUT
../one_d_womersley -validation_run -outflow 2 >> OUTPUT
echo "done"
echo " " >> validation.log
echo "1D Womersley validation" >> validation.log
echo "-----------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log



echo "1D impedance tube test: " >> validation.log
cat  RESLT_impedance_tube/soln0.dat RESLT_impedance_tube/soln1.dat RESLT_impedance_tube/soln9.dat RESLT_impedance_tube/soln10.dat RESLT_impedance_tube/trace.dat \
> 1_d_womersley.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py \
../validata/1_d_womersley.dat.gz  \
1_d_womersley.dat >> validation.log
fi


echo "1D Womersley with flux control test: " >> validation.log
cat  RESLT_impedance_tube_with_flux_control/soln0.dat RESLT_impedance_tube_with_flux_control/soln1.dat RESLT_impedance_tube_with_flux_control/soln9.dat RESLT_impedance_tube_with_flux_control/soln10.dat RESLT_impedance_tube_with_flux_control/trace.dat \
> 1_d_womersley_with_flux_control.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py \
../validata/1_d_womersley.dat.gz  \
1_d_womersley_with_flux_control.dat >> validation.log
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
