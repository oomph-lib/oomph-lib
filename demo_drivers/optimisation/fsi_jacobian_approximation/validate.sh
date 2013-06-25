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
cd Validation

# Validation for approximate FSI Jacobian (steady)
#-------------------------------------------------
echo "Running with/without approximate FSI Jacobian (steady)"
mkdir RESLT_APPROX_STEADY RESLT_EXACT_STEADY
../fsi_jacobian_approximation -validation_run -approx_jacobian -steady_run > OUTPUT_APPROX_STEADY
../fsi_jacobian_approximation -validation_run -steady_run> OUTPUT_EXACT_STEADY
echo "done"
echo " " >> validation.log
echo "With/without approximate FSI Jacobian (steady)" >> validation.log
echo "----------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

echo "Approximate FSI Jacobian test (steady): " >> validation.log
cat  RESLT_APPROX_STEADY/soln1.dat RESLT_APPROX_STEADY/soln2.dat \
RESLT_APPROX_STEADY/soln3.dat RESLT_APPROX_STEADY/trace.dat \
RESLT_EXACT_STEADY/soln1.dat RESLT_EXACT_STEADY/soln2.dat \
RESLT_EXACT_STEADY/soln3.dat RESLT_EXACT_STEADY/trace.dat \
> fsi_jacobian_approximation_steady.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py \
../validata/fsi_jacobian_approximation_steady.dat.gz  \
fsi_jacobian_approximation_steady.dat  0.1 1.0e-11 >> validation.log
fi

# Validation for approximate FSI Jacobian (unsteady)
#---------------------------------------------------
echo "Running with/without approximate FSI Jacobian (unsteady)"
mkdir RESLT_APPROX_UNSTEADY RESLT_EXACT_UNSTEADY
../fsi_jacobian_approximation -validation_run -approx_jacobian > OUTPUT_APPROX_UNSTEADY
../fsi_jacobian_approximation -validation_run > OUTPUT_EXACT_UNSTEADY
echo "done"
echo " " >> validation.log
echo "With/without approximate FSI Jacobian (unsteady)" >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

echo "Approximate FSI Jacobian test (unsteady): " >> validation.log
cat  RESLT_APPROX_UNSTEADY/soln1.dat RESLT_APPROX_UNSTEADY/soln2.dat \
RESLT_APPROX_UNSTEADY/soln9.dat RESLT_APPROX_UNSTEADY/soln10.dat \
RESLT_APPROX_UNSTEADY/trace.dat \
RESLT_EXACT_UNSTEADY/soln1.dat RESLT_EXACT_UNSTEADY/soln2.dat \
RESLT_EXACT_UNSTEADY/soln9.dat RESLT_EXACT_UNSTEADY/soln10.dat \
RESLT_EXACT_UNSTEADY/trace.dat \
> fsi_jacobian_approximation_unsteady.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py \
../validata/fsi_jacobian_approximation_unsteady.dat.gz  \
fsi_jacobian_approximation_unsteady.dat 0.1 1.0e-11 >> validation.log
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
