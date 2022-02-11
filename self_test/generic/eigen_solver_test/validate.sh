#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
if test "$MPI_RUN_COMMAND" = ""; then
NUM_TESTS=3
else
NUM_TESTS=3
fi

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation

# Validation for eigensolver test
#-----------------------------------------
echo "Running eigensolver validation "
mkdir RESLT_lapack
mkdir RESLT_anasazi
../eigen_solver_test > OUTPUT
echo "done"
echo " " >> validation.log
echo "Eigensolver validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_lapack/* > eigen_solver_test_lapack.dat
cat RESLT_anasazi/* > eigen_solver_test_anasazi.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/eigen_solver_test_lapack.dat.gz  \
         eigen_solver_test_lapack.dat >> validation.log
if test "$MPI_RUN_COMMAND" = ""; then
../../../../bin/fpdiff.py ../validata/eigen_solver_test_anasazi.dat.gz  \
         eigen_solver_test_anasazi.dat >> validation.log
fi
fi
rm -rf RESLT_lapack RESLT_anasazi
#-----------------------------------------

# Validation for eigensolver test
#-----------------------------------------
echo "Running find eigenvalues validation "
mkdir RESLT
../lapack_qz_find_eigenvalues_test > OUTPUT
echo "done"
echo " " >> validation.log
echo "Find eigenvalues validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/* > lapack_qz_find_eigenvalues_test.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/lapack_qz_find_eigenvalues_test.dat.gz  \
         lapack_qz_find_eigenvalues_test.dat >> validation.log
fi
rm -rf RESLT
#-----------------------------------------

if test "$MPI_RUN_COMMAND" != ""; then
# Validation for eigensolver test
#-----------------------------------------
echo "Running ANASAZI distributed validation "
mkdir RESLT_anasazi_distributed
mpirun -n 2 ../eigen_solver_test_distributed > OUTPUT
echo "done"
echo " " >> validation.log
echo "Find eigenvalues validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_anasazi_distributed/* > eigen_solver_test_distributed_anasazi.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/eigen_solver_test_distributed_anasazi.dat.gz  \
         eigen_solver_test_distributed_anasazi.dat >> validation.log
fi
rm -rf RESLT
#-----------------------------------------
# Append log to main validation log
cat validation.log >> ../../../../validation.log
fi

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
