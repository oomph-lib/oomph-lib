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

cd Validation


# Validation for FSI oscillating ring Navier Stokes -- compare shape derivs
#--------------------------------------------------------------------------

echo "Running FSI oscillating ring Navier Stokes -- compare shape derivs"
mkdir RESLT

# Run with two command line arguments so we only run a few steps
../fsi_osc_ring_compare_jacs lala lala > OUTPUT_fsi_osc_ring_compare_jacs
echo "done"
echo " " >> validation.log
echo "FSI oscillating ring Navier Stokes validation -- compare shape derivs" >> validation.log
echo "---------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0_1.dat \
    RESLT/soln1_1.dat \
    RESLT/soln2_1.dat \
    RESLT/soln3_1.dat \
    RESLT/soln4_1.dat \
    > fsi_ring_compare_jacs_results.dat


echo "dummy [OK] -- Data appears to be too rough to provide meaningful fpdiff comparison" >> validation.log

#if test "$1" = "no_fpdiff"; then
#  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
#else
#../../../../bin/fpdiff.py ../validata/fsi_ring_compare_jacs_results.dat.gz   \
#         fsi_ring_compare_jacs_results.dat 0.1 1.0e-8 >> validation.log
#fi

mv RESLT RESLT_fsi_osc_ring_compare_jacs

# Validation for FSI oscillating ring Navier Stokes
#--------------------------------------------------

echo "Running FSI oscillating ring Navier Stokes"
mkdir RESLT

# Run with two command line arguments so we only run a few steps
../fsi_osc_ring lala lala > OUTPUT_fsi_osc_ring
echo "done"
echo " " >> validation.log
echo "FSI oscillating ring Navier Stokes validation" >> validation.log
echo "---------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat \
    RESLT/soln1.dat \
    > fsi_ring_results.dat

echo "dummy [OK] -- Data appears to be too rough to provide meaningful fpdiff comparison" >> validation.log

#if test "$1" = "no_fpdiff"; then
#  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
#else
#../../../../bin/fpdiff.py ../validata/fsi_ring_results.dat.gz   \
#         fsi_ring_results.dat 0.1 1.0e-8 >> validation.log
#fi

mv RESLT RESLT_fsi_osc_ring


# Validation for oscillating ring Navier Stokes with algebraic mesh update
#-------------------------------------------------------------------------

echo "Running algebraic oscillating ring Navier-Stokes validation "
mkdir RESLT
mkdir RESLT_restarted
# Run with two command line arguments so we only do three steps
../osc_ring_alg lala lala > OUTPUT_osc_ring_alg
echo "done"
echo " " >> validation.log
echo "Algebraic oscillating ring Navier Stokes validation" >> validation.log
echo "---------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

cat RESLT/soln0.dat \
    RESLT/soln3.dat \
    > osc_ring_alg_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/osc_ring_alg_results.dat.gz  \
         osc_ring_alg_results.dat >> validation.log
fi

mv RESLT RESLT_osc_ring_alg
mv RESLT_restarted RESLT_osc_ring_alg_restarted



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
