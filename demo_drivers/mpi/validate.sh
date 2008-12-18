#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=16


# Doc what we're using to run tests on two processors
echo " " 
echo "Running mpi tests with mpi run command: " $MPI_RUN_COMMAND
echo " " 

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation


# Trilinos test
#==============


# Validation for Trilinos
#------------------------
if [ -f ../TrilinosSolver_test ]; then

echo "Running trilinos tests in parallel on two processors"
mkdir RESLT
$MPI_RUN_COMMAND ../TrilinosSolver_test > OUTPUT_Trilinos
echo "done"
echo " " >> validation.log
echo "Trilinos tests in parallel on two processors" >> validation.log
echo "--------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "Solver/preconditioner combination 0" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln0.dat >> validation.log
echo "Solver/preconditioner combination 1" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln1.dat >> validation.log
echo "Solver/preconditioner combination 2" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln2.dat >> validation.log
echo "Solver/preconditioner combination 3" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln3.dat >> validation.log
echo "Solver/preconditioner combination 4" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln4.dat >> validation.log
echo "Solver/preconditioner combination 5" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln5.dat >> validation.log
echo "Solver/preconditioner combination 6" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln6.dat >> validation.log
echo "Solver/preconditioner combination 7" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln7.dat >> validation.log
echo "Solver/preconditioner combination 8" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln8.dat >> validation.log
echo "Solver/preconditioner combination 9" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln9.dat >> validation.log
echo "Solver/preconditioner combination 10" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln10.dat >> validation.log
echo "Solver/preconditioner combination 11" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln11.dat >> validation.log
echo "Solver/preconditioner combination 12" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln12.dat >> validation.log
echo "Solver/preconditioner combination 13" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln13.dat >> validation.log
echo "Solver/preconditioner combination 14" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln14.dat >> validation.log
echo "Number of Newton iterations for Trilinos solves" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_conv.dat.gz  \
         RESLT/conv.dat >> validation.log

fi

mv RESLT RESLT_Trilinos

else

echo "Not running Trilinos tests as executable doesn't exist"
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
echo "[OK] (Dummy for non-existent Trilinos)"  >> validation.log
fi






# Append log to main validation log
cat validation.log >> ../../../validation.log

cd ..



#######################################################################


#Check that we get the correct number of OKs
OK_COUNT=`grep -c 'OK' Validation/validation.log`
if  [ $OK_COUNT -eq $NUM_TESTS ]; then
 echo " "
 echo "======================================================================"
 echo " " 
 echo "All tests in" 
 echo " " 
 echo "    `pwd`    "
 echo " "
 echo "passed successfully."
 echo " "
 echo "======================================================================"
 echo " " 
else
  if [ $OK_COUNT -lt $NUM_TESTS ]; then
   echo " "
   echo "======================================================================"
   echo " " 
   echo "Only $OK_COUNT of $NUM_TESTS test(s) passed; see"
   echo " " 
   echo "    `pwd`/Validation/validation.log"
   echo " " 
   echo "for details" 
   echo " " 
   echo "======================================================================"
   echo " "
  else 
   echo " "
   echo "======================================================================"
   echo " " 
   echo "More OKs than tests! Need to update NUM_TESTS in"
   echo " " 
   echo "    `pwd`/validate.sh"
   echo " "
   echo "======================================================================"
   echo " "
  fi
fi
