#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=3


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for non-adaptive case - triangular mesh
#---------------------------------------------------
cd Validation

echo "Running PML time-harmonic linear elasticity "
mkdir RESLT
../time_harmonic_elasticity_driver --l_pml 1.6 --n_pml 4 --max_adapt 0 --validation > OUTPUT_non_adapt
echo "done"
echo " " >> validation.log
echo "Generalised time-periodic linear elasticity" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/elast_soln_norm0.dat > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/results.dat.gz   \
    results.dat   >> validation.log
fi
mv RESLT RESLT_non_adapt

# Validation for adaptive case - triangular mesh
#-----------------------------------------------
echo "Running adaptive PML time-harmonic linear elasticity "
mkdir RESLT
../time_harmonic_elasticity_driver --l_pml 1.6 --n_pml 4 --max_adapt 1 --validation > OUTPUT_adapt
echo "done"
echo " " >> validation.log
echo "Adaptive PML time-harmonic linear elasticity" >> validation.log
echo "---------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/elast_soln_norm0.dat > adaptive_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/adaptive_results.dat.gz   \
    adaptive_results.dat   >> validation.log
fi
mv RESLT RESLT_adapt

# Validation for non-adaptive case - rectangular mesh
#---------------------------------------------------
echo "Running quad bulk element PML time-harmonic linear elasticity"
mkdir RESLT
../time_harmonic_elasticity_driver_source --l_pml 4.0 --n_pml 4 --validation > OUTPUT_source
echo "done"
echo " " >> validation.log
echo "Quad bulk element PML time-harmonic linear elasticity" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/elast_soln_norm0.dat > rectangular_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/rectangular_results.dat.gz   \
    rectangular_results.dat   >> validation.log
fi
mv RESLT RESLT_rectangular


# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../validation.log


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
