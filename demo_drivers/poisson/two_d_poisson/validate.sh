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

# Validation for demo poisson
#----------------------------
cd Validation

echo "Running 2D demo poisson validation "
mkdir RESLT
../two_d_poisson > OUTPUT_two_d_poisson
echo "done"
echo " " >> validation.log
echo "2D demo poisson validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat \
    RESLT/soln3.dat  > two_d_poisson_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/two_d_poisson_results.dat.gz   \
    two_d_poisson_results.dat  >> validation.log
fi

mv RESLT RESLT_noprecompiled_mesh







# Validation for comparison of solvers
#-------------------------------------

echo "Running comparison of solvers for 2D demo poisson validation "

mkdir RESLT_cr
mkdir RESLT_cc
mkdir RESLT_frontal
mkdir RESLT_frontal_reordered
mkdir RESLT_dense_LU
mkdir RESLT_FD_LU

../two_d_poisson_compare_solvers > OUTPUT_compare_solvers
echo "done"
echo " " >> validation.log
echo "Comparison of solvers with 2D demo poisson validation" >> validation.log
echo "-----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

# These output files always exist:
cat RESLT_cr/soln0.dat RESLT_cc/soln0.dat \
    RESLT_dense_LU/soln0.dat RESLT_FD_LU/soln0.dat > \
    compare_solvers_results.dat

# Output from frontal solver only if it exists
HAVE_FRONTAL_RESULTS=false
if test -e RESLT_frontal_reordered/soln0.dat; then
    HAVE_FRONTAL_RESULTS=true
    cat RESLT_frontal_reordered/soln0.dat >> compare_solvers_results.dat
fi
if test -e RESLT_frontal/soln0.dat; then
    HAVE_FRONTAL_RESULTS=true
    cat RESLT_frontal/soln0.dat >> compare_solvers_results.dat
fi


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
if $HAVE_FRONTAL_RESULTS; then
../../../../bin/fpdiff.py ../validata/compare_solvers_results.dat.gz   \
    compare_solvers_results.dat  >> validation.log
else
../../../../bin/fpdiff.py \
    ../validata/compare_solvers_results2.dat.gz   \
    compare_solvers_results.dat  >> validation.log
fi
fi


# Append output to global validation log file
#--------------------------------------------
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
