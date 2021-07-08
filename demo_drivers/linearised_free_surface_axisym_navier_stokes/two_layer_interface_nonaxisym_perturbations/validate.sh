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

# Validation for two-layer interface nonaxisymmetric perturbations problem
#-------------------------------------------------------------------------
cd Validation

echo "Running two-layer interface nonaxisymmetric perturbations validation "
mkdir RESLT
../two_layer_interface_nonaxisym_perturbations lala > OUTPUT
echo "done"
echo " " >> validation.log
echo "Two-layer interface nonaxisymmetric perturbations validation" >> validation.log
echo "------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/base_soln0.dat \
    RESLT/base_soln1.dat \
    RESLT/base_soln2.dat \
    RESLT/perturbed_k0_soln0.dat \
    RESLT/perturbed_k0_soln1.dat \
    RESLT/perturbed_k0_soln2.dat \
    RESLT/perturbed_k1_soln0.dat \
    RESLT/perturbed_k1_soln1.dat \
    RESLT/perturbed_k1_soln2.dat \
    RESLT/perturbed_k2_soln0.dat \
    RESLT/perturbed_k2_soln1.dat \
    RESLT/perturbed_k2_soln2.dat \
    RESLT/perturbation_to_interface_k0_soln0.dat \
    RESLT/perturbation_to_interface_k0_soln1.dat \
    RESLT/perturbation_to_interface_k0_soln2.dat \
    RESLT/perturbation_to_interface_k1_soln0.dat \
    RESLT/perturbation_to_interface_k1_soln1.dat \
    RESLT/perturbation_to_interface_k1_soln2.dat \
    RESLT/perturbation_to_interface_k2_soln0.dat \
    RESLT/perturbation_to_interface_k2_soln1.dat \
    RESLT/perturbation_to_interface_k2_soln2.dat > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz  \
         results.dat >> validation.log
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
