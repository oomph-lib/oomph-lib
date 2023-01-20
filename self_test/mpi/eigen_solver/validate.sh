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

cd Validation

# Validation for eigensolver test
#--------------------------------

if [ -f ../eigen_solver_test ]; then
    
echo "Running eigensolver validation "
mkdir RESLT_anasazi
cp ../random_test_matrix.dat .
mpirun -n 1 ../eigen_solver_test > OUTPUT
echo "done"
echo " " >> validation.log
echo "Eigensolver validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_anasazi/* > eigen_solver_test_anasazi.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/eigen_solver_test_anasazi.dat.gz  \
         eigen_solver_test_anasazi.dat >> validation.log
fi
rm -rf RESLT_anasazi


else
    
echo "Not running eigensolver validation; needs trilinos "
echo "Dummy [OK] for missing trilinos " >> validation.log

fi




#-----------------------------------------

#-----------------------------------------
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
