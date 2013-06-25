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

# Validation for Womersley problems
#----------------------------------
cd Validation

echo "Running 2D Womersley validation "
mkdir RESLT_prescribed_pressure_gradient
mkdir RESLT_prescribed_volume_flux
mkdir RESLT_impedance_tube
mkdir RESLT_navier_stokes
../two_d_womersley > OUTPUT
echo "done"
echo " " >> validation.log
echo "2D Womersley validation" >> validation.log
echo "-----------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log



echo "Prescribed pressure gradient test: " >> validation.log
cat  RESLT_prescribed_pressure_gradient/womersley_soln0.dat RESLT_prescribed_pressure_gradient/womersley_soln1.dat RESLT_prescribed_pressure_gradient/womersley_soln2.dat RESLT_prescribed_pressure_gradient/womersley_soln10.dat RESLT_prescribed_pressure_gradient/trace.dat \
     > presc_pres_grad.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py \
../validata/presc_pres_grad.dat.gz  \
presc_pres_grad.dat >> validation.log
fi



echo "Prescribed volume flux test: " >> validation.log
cat  RESLT_prescribed_volume_flux/womersley_soln0.dat RESLT_prescribed_volume_flux/womersley_soln1.dat RESLT_prescribed_volume_flux/womersley_soln2.dat RESLT_prescribed_volume_flux/womersley_soln10.dat   RESLT_prescribed_volume_flux/trace.dat \
 > prescribed_volume_flux_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py \
../validata/prescribed_volume_flux_results.dat.gz  \
prescribed_volume_flux_results.dat 0.1 2.0e-13 >> validation.log
fi



echo "Impedance tube test: " >> validation.log
cat  RESLT_impedance_tube/womersley_soln0.dat RESLT_impedance_tube/womersley_soln1.dat RESLT_impedance_tube/womersley_soln2.dat  RESLT_impedance_tube/womersley_soln10.dat RESLT_impedance_tube/trace.dat \
 > impedance_tube_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py \
../validata/impedance_tube_results.dat.gz  \
impedance_tube_results.dat 0.1 2.0e-13 >> validation.log
fi


echo "Navier Stokes outflow test: " >> validation.log
cat  RESLT_navier_stokes/womersley_soln0.dat RESLT_navier_stokes/womersley_soln1.dat RESLT_navier_stokes/womersley_soln2.dat  RESLT_navier_stokes/womersley_soln10.dat RESLT_navier_stokes/trace.dat \
 > navier_stokes_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py \
../validata/navier_stokes_results.dat.gz  \
navier_stokes_results.dat 0.1 1.0e-13 >> validation.log
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
