#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=26

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

# Validation for rectangular driven cavity
#-----------------------------------------

echo "Running rectangular driven cavity LSC precond validation "
mkdir RESLT
$MPI_RUN_COMMAND ../driven_cavity > OUTPUT_driven_cavity
echo "done"
echo " " >> validation.log
echo "Rectangular driven cavity validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0.dat RESLT/soln1.dat \
 > driven_cavity_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/driven_cavity_results.dat.gz  \
         driven_cavity_results.dat >> validation.log
fi

# Airy Cantiliver
#================
echo "Running Airy Cantilever with BlockDiagonalPreconditioner"
$MPI_RUN_COMMAND ../airy_cantilever 0 > OUTPUT_airy_0
echo "done"
echo ""
echo "Running Airy Cantilever with BlockDiagonalPreconditioner with two level parallelisation"
$MPI_RUN_COMMAND ../airy_cantilever 1 > OUTPUT_airy_1
echo "done"
echo ""
echo "Running Airy Cantilever with BlockTriangularPreconditioner (Upper)"
$MPI_RUN_COMMAND ../airy_cantilever 2 > OUTPUT_airy_2
echo "done"
echo ""
echo "Running Airy Cantilever with BlockTriangularPreconditioner (Lower)"
$MPI_RUN_COMMAND ../airy_cantilever 3 > OUTPUT_airy_3
echo "done"
echo ""
echo " " >> validation.log
echo "Airy Cantilever block preconditioner tests" >> validation.log
echo "--------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/airy_soln.dat.gz \
    RESLT/airy_soln0.dat 0.1 10e-6 >> validation.log
fi
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/airy_soln.dat.gz \
    RESLT/airy_soln1.dat 0.1 10e-6 >> validation.log
fi
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/airy_soln.dat.gz \
    RESLT/airy_soln2.dat 0.1 10e-6 >> validation.log
fi
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/airy_soln.dat.gz \
    RESLT/airy_soln3.dat 0.1 10e-6 >> validation.log
fi

# Trilinos test
#==============


# Validation for Trilinos
#------------------------
if [ -f ../TrilinosSolver_test ]; then

echo "Running trilinos tests"
$MPI_RUN_COMMAND ../TrilinosSolver_test > OUTPUT_Trilinos
echo "done"
echo " " >> validation.log
echo "Trilinos tests" >> validation.log
echo "--------------" >> validation.log
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
         RESLT/soln_trilinos_0.dat >> validation.log
echo "Solver/preconditioner combination 1" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_1.dat >> validation.log
echo "Solver/preconditioner combination 2" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_2.dat >> validation.log
echo "Solver/preconditioner combination 3" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_3.dat >> validation.log
echo "Solver/preconditioner combination 4" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_4.dat >> validation.log
echo "Solver/preconditioner combination 5" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_5.dat >> validation.log
echo "Solver/preconditioner combination 6" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_6.dat >> validation.log
echo "Solver/preconditioner combination 7" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_7.dat >> validation.log
echo "Solver/preconditioner combination 8" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_8.dat >> validation.log
echo "Solver/preconditioner combination 9" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_9.dat >> validation.log
echo "Solver/preconditioner combination 10" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_10.dat >> validation.log
echo "Solver/preconditioner combination 11" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_11.dat >> validation.log
echo "Solver/preconditioner combination 12" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_12.dat >> validation.log
echo "Solver/preconditioner combination 13" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_13.dat >> validation.log
echo "Solver/preconditioner combination 14" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_14.dat >> validation.log
echo "Number of Newton iterations for Trilinos solves" >> validation.log
../../..//bin/fpdiff.py ../validata/Trilinos_conv.dat.gz  \
         RESLT/conv.dat >> validation.log

fi


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





# Direct Solver Tests
#====================


# Validation for direct solver tests
#-----------------------------------
echo "Running direct solver tests"
$MPI_RUN_COMMAND ../direct_solver_test > OUTPUT_DirectSolver
echo "done"
echo " " >> validation.log
echo "Direct solver tests" >> validation.log
echo "-------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_python"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python" >> validation.log
else
echo "SuperLU_dist matrix based solve w/ global CRDoubleMatrixSolver" >> validation.log
../../..//bin/fpdiff.py ../validata/direct_solver_matrix_solve_result.dat.gz  \
         RESLT/SuperLU_dist_CRDoubleMatrix_global.dat >> validation.log
echo "SuperLU_dist matrix based solve w/ dist CRDoubleMatrixSolver" >> validation.log
../../..//bin/fpdiff.py ../validata/direct_solver_matrix_solve_result.dat.gz  \
         RESLT/SuperLU_dist_CRDoubleMatrix_distributed.dat >> validation.log
echo "SuperLU matrix based solve w/ CCDoubleMatrixSolver" >> validation.log
../../..//bin/fpdiff.py ../validata/direct_solver_matrix_solve_result.dat.gz  \
         RESLT/SuperLU_dist_CCDoubleMatrix.dat >> validation.log
echo "SuperLU_dist global problem based solve" >> validation.log
../../..//bin/fpdiff.py ../validata/direct_solver_problem_solve_result.dat.gz  \
         RESLT/soln_direct_solver_1.dat >> validation.log
echo "SuperLU_dist global problem based solve" >> validation.log
../../..//bin/fpdiff.py ../validata/direct_solver_problem_solve_result.dat.gz  \
         RESLT/soln_direct_solver_2.dat >> validation.log
fi

mv RESLT RESLT_test

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
