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

# Validation for Fp/PCD preconditioner
#-------------------------------------
cd Validation

echo "Running 2D Fp/PCD preconditioner test "
mkdir RESLT
../two_d_fp bla > OUTPUT_two_d_fp
echo "done"
echo " " >> validation.log
echo "2D Fp/PCD preconditioner test" >> validation.log
echo "-----------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

grep iterations OUTPUT_two_d_fp | awk '{print $NF}' > iter_two_d.dat

cat \
RESLT/soln0.dat   RESLT/soln16.dat RESLT/soln4.dat \
RESLT/soln10.dat  RESLT/soln17.dat RESLT/soln5.dat \
RESLT/soln11.dat  RESLT/soln18.dat RESLT/soln6.dat \
RESLT/soln12.dat  RESLT/soln19.dat RESLT/soln7.dat \
RESLT/soln13.dat  RESLT/soln1.dat  RESLT/soln8.dat \
RESLT/soln14.dat  RESLT/soln2.dat  RESLT/soln9.dat \
RESLT/soln15.dat  RESLT/soln3.dat \
     > two_d_fp.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else

echo "Results: " >> validation.log

../../../../bin/fpdiff.py ../validata/two_d_fp.dat.gz  \
         two_d_fp.dat 0.1 1.0e-12 >> validation.log

echo "Iteration counts (allowing for 20% variation): " >> validation.log
echo " " >> validation.log
echo "[Note: This test fails if  " >> validation.log
echo " " >> validation.log
echo "       (1) the number of linear solver iterations " >> validation.log
echo "           differs by more than 20% " >> validation.log
echo " " >> validation.log
echo "       (2) the linear solver doesn't provide a " >> validation.log
echo "           good enough solution, causing an " >> validation.log
echo "           increase in the number of Newton" >> validation.log
echo "           iterations. " >> validation.log
echo " " >> validation.log
echo "       In the latter case, not even the number of " >> validation.log
echo "       lines in the iter* files will match.] " >> validation.log


../../../../bin/fpdiff.py ../validata/iter_two_d.dat.gz  \
         iter_two_d.dat 20.0 1.0E-10 >> validation.log

fi

mv RESLT RESLT_two_d_fp


#######################################################################



echo "Running 3D Fp/PCD preconditioner test "
mkdir RESLT
../three_d_fp bla > OUTPUT_three_d_fp
echo "done"
echo " " >> validation.log
echo "3D Fp/PCD preconditioner test" >> validation.log
echo "-----------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

grep iterations OUTPUT_three_d_fp | awk '{print $NF}' > iter_three_d.dat

cat \
RESLT/soln0.dat   RESLT/soln4.dat \
RESLT/soln10.dat  RESLT/soln5.dat \
RESLT/soln11.dat  RESLT/soln6.dat \
RESLT/soln1.dat   RESLT/soln7.dat \
RESLT/soln2.dat   RESLT/soln8.dat \
RESLT/soln3.dat   RESLT/soln9.dat \
     > three_d_fp.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else

echo "Results: " >> validation.log

../../../../bin/fpdiff.py ../validata/three_d_fp.dat.gz  \
         three_d_fp.dat 0.1 1.0e-9 >> validation.log

echo "Iteration counts (allowing for 50% variation): " >> validation.log
echo " " >> validation.log
echo "[Note: This test fails if  " >> validation.log
echo " " >> validation.log
echo "       (1) the number of linear solver iterations " >> validation.log
echo "           differs by more than 20% " >> validation.log
echo " " >> validation.log
echo "       (2) the linear solver doesn't provide a " >> validation.log
echo "           good enough solution, causing an " >> validation.log
echo "           increase in the number of Newton" >> validation.log
echo "           iterations. " >> validation.log
echo " " >> validation.log
echo "       In the latter case, not even the number of " >> validation.log
echo "       lines in the iter* files will match.] " >> validation.log


../../../../bin/fpdiff.py ../validata/iter_three_d.dat.gz  \
         iter_three_d.dat 50.0 1.0E-10>> validation.log

fi

mv RESLT RESLT_three_d_fp





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
