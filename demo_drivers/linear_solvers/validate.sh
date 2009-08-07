#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=65

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
../driven_cavity > OUTPUT_driven_cavity
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


mv RESLT RESLT_driven_cavity



# Validation with DenseLU
#------------------------
echo "Running Poisson linear solver test with DenseLU "
mkdir RESLT
../two_d_poisson_dense_lu > OUTPUT_DenseLU
echo "done"
echo " " >> validation.log
echo "Adv diff linear solver test with DenseLU validation" >> validation.log
echo "---------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat  RESLT/soln3.dat > DenseLU_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/DenseLU_results.dat.gz  \
         DenseLU_results.dat >> validation.log
fi


mv RESLT RESLT_DenseLU

# Validation with SuperLU
#------------------------
echo "Running adv diff linear solver test with SuperLU "
mkdir RESLT
../adv_diff_iterative_linear_solver_tester 0 > OUTPUT_SuperLU
echo "done"
echo " " >> validation.log
echo "Adv diff linear solver test with SuperLU validation" >> validation.log
echo "---------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > SuperLU_results.dat
cat RESLT/resolve_error.dat RESLT/la_solve_error.dat > SuperLU_solve_error.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/results_Pe200.dat.gz  \
         SuperLU_results.dat >> validation.log
fi

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/SuperLU_solve_error.dat.gz  \
          SuperLU_solve_error.dat >> validation.log
fi

mv RESLT RESLT_SuperLU


# Threshold for number of iterations in comparison of convergence histories
#===========================================================================
threshold_for_number_of_iterations=3


# Validation with GMRES
#----------------------
echo "Running adv diff linear solver test with GMRES "
mkdir RESLT
../adv_diff_iterative_linear_solver_tester 1 > OUTPUT_GMRES
echo "done"
echo " " >> validation.log
echo "Adv diff linear solver test with GMRES validation" >> validation.log
echo "---------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > GMRES_results.dat
cat RESLT/convergence.dat > GMRES_convergence.dat
cat RESLT/la_solve_convergence.dat > GMRES_la_solve_convergence.dat
cat RESLT/resolve_error.dat RESLT/la_solve_error.dat > GMRES_solve_error.dat


if test "$1" = "no_fpdiff"; then

    echo "dummy [OK] -- Can't run fpdiff.py because we don't have validata" \
        >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have validata" \
        >> validation.log

else

    #Compare number of iterations against reference data and append
    ../../../bin/compare_file_length_with_tolerance.bash GMRES_convergence.dat\
    ../validata/GMRES_convergence.dat $threshold_for_number_of_iterations \
    >>  validation.log

    #Compare number of iterations against reference data and append
    ../../../bin/compare_file_length_with_tolerance.bash \ 
    GMRES_la_solve_convergence.dat \
        ../validata/GMRES_la_solve_convergence.dat \
        $threshold_for_number_of_iterations \
        >>  validation.log
fi


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/results_Pe200.dat.gz  \
         GMRES_results.dat >> validation.log
fi

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/GMRES_solve_error.dat.gz  \
          GMRES_solve_error.dat 0.1 1.0e-12 >> validation.log
fi

mv RESLT RESLT_GMRES



# Validation with BiCGStab
#------------------------
echo "Running adv diff linear solver test with BiCGStab "
mkdir RESLT
../adv_diff_iterative_linear_solver_tester 2 > OUTPUT_BiCGStab
echo "done"
echo " " >> validation.log
echo "Adv diff linear solver test with BiCGStab validation" >> validation.log
echo "---------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > BiCGStab_results.dat
cat RESLT/convergence.dat > BiCGStab_convergence.dat
cat RESLT/la_solve_convergence.dat > BiCGStab_la_solve_convergence.dat
cat RESLT/resolve_error.dat RESLT/la_solve_error.dat > BiCGStab_solve_error.dat

if test "$1" = "no_fpdiff"; then

    echo "dummy [OK] -- Can't run fpdiff.py because we don't have validata" \
        >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have validata" \
        >> validation.log

else

    #Compare number of iterations against reference data and append
    ../../../bin/compare_file_length_with_tolerance.bash \
        BiCGStab_convergence.dat \
        ../validata/BiCGStab_convergence.dat \
        $threshold_for_number_of_iterations \
        >>  validation.log

    #Compare number of iterations against reference data and append
    ../../../bin/compare_file_length_with_tolerance.bash \
        BiCGStab_la_solve_convergence.dat \
        ../validata/BiCGStab_la_solve_convergence.dat \
        $threshold_for_number_of_iterations \
        >>  validation.log

fi 


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/results_Pe200.dat.gz  \
         BiCGStab_results.dat >> validation.log
fi


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/BiCGStab_solve_error.dat.gz  \
          BiCGStab_solve_error.dat 0.1 5.0e-12 >> validation.log
fi

mv RESLT RESLT_BiCGStab





# Validation with CG
#-------------------
echo "Running adv diff linear solver test with CG "
mkdir RESLT
../adv_diff_iterative_linear_solver_tester 3 > OUTPUT_CG
echo "done"
echo " " >> validation.log
echo "Adv diff linear solver test with CG validation" >> validation.log
echo "---------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > CG_results.dat
cat RESLT/convergence.dat > CG_convergence.dat
cat RESLT/la_solve_convergence.dat > CG_la_solve_convergence.dat
cat RESLT/resolve_error.dat RESLT/la_solve_error.dat > CG_solve_error.dat

if test "$1" = "no_fpdiff"; then

    echo "dummy [OK] -- Can't run fpdiff.py because we don't have validata" \
        >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have validata" \
        >> validation.log

else

    #Compare number of iterations against reference data and append
    ../../../bin/compare_file_length_with_tolerance.bash CG_convergence.dat \
        ../validata/CG_convergence.dat $threshold_for_number_of_iterations \
        >>  validation.log
    
    #Compare number of iterations against reference data and append
    ../../../bin/compare_file_length_with_tolerance.bash \
        CG_la_solve_convergence.dat \
        ../validata/CG_la_solve_convergence.dat \
        $threshold_for_number_of_iterations \
        >>  validation.log

fi

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/results_Pe0.dat.gz  \
         CG_results.dat >> validation.log
fi


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/CG_solve_error.dat.gz  \
         CG_solve_error.dat >> validation.log
fi


mv RESLT RESLT_CG




# Validation with GS
#-------------------
echo "Running adv diff linear solver test with GS "
mkdir RESLT
../adv_diff_iterative_linear_solver_tester 4 > OUTPUT_GS
echo "done"
echo " " >> validation.log
echo "Adv diff linear solver test with GS validation" >> validation.log
echo "---------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat > GS_results.dat
cat RESLT/convergence.dat > GS_convergence.dat
cat RESLT/la_solve_convergence.dat > GS_la_solve_convergence.dat
cat RESLT/resolve_error.dat RESLT/la_solve_error.dat > GS_solve_error.dat


    
if test "$1" = "no_fpdiff"; then

    echo "dummy [OK] -- Can't run fpdiff.py because we don't have validata" \
        >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have validata" \
        >> validation.log

else

    #Compare number of iterations against reference data and append
    ../../../bin/compare_file_length_with_tolerance.bash GS_convergence.dat \
        ../validata/GS_convergence.dat $threshold_for_number_of_iterations \
        >>  validation.log
    
    
    #Compare number of iterations against reference data and append
    ../../../bin/compare_file_length_with_tolerance.bash \
        GS_la_solve_convergence.dat \
        ../validata/GS_la_solve_convergence.dat \
        $threshold_for_number_of_iterations \
        >>  validation.log

fi


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/results_Pe0.dat.gz  \
         GS_results.dat >> validation.log
fi


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../bin/fpdiff.py ../validata/GS_solve_error.dat.gz  \
          GS_solve_error.dat >> validation.log
fi


mv RESLT RESLT_GS




# Trilinos test
#==============


# Validation for Trilinos
#------------------------
if [ -f ../TrilinosSolver_test ]; then

echo "Running trilinos tests"
mkdir RESLT
../TrilinosSolver_test > OUTPUT_Trilinos
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
../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln0.dat >> validation.log
echo "Solver/preconditioner combination 1" >> validation.log
../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln1.dat >> validation.log
echo "Solver/preconditioner combination 2" >> validation.log
../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln2.dat >> validation.log
echo "Solver/preconditioner combination 3" >> validation.log
../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln3.dat >> validation.log
echo "Solver/preconditioner combination 4" >> validation.log
../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln4.dat >> validation.log
echo "Solver/preconditioner combination 5" >> validation.log
../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln5.dat >> validation.log
echo "Solver/preconditioner combination 6" >> validation.log
../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln6.dat >> validation.log
echo "Solver/preconditioner combination 7" >> validation.log
../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln7.dat >> validation.log
echo "Solver/preconditioner combination 8" >> validation.log
../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln8.dat >> validation.log
echo "Solver/preconditioner combination 9" >> validation.log
../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln9.dat >> validation.log
echo "Solver/preconditioner combination 10" >> validation.log
../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln10.dat >> validation.log
echo "Solver/preconditioner combination 11" >> validation.log
../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln11.dat >> validation.log
echo "Solver/preconditioner combination 12" >> validation.log
../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln12.dat >> validation.log
echo "Solver/preconditioner combination 13" >> validation.log
../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln13.dat >> validation.log
echo "Solver/preconditioner combination 14" >> validation.log
../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln14.dat >> validation.log
echo "Number of Newton iterations for Trilinos solves" >> validation.log
../../../bin/fpdiff.py ../validata/Trilinos_conv.dat.gz  \
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





# Hypre test
#===========


# Validation for Hypre
#----------------------
if [ -f ../HypreSolver_test ]; then

echo "Running hypre tests"
mkdir RESLT
../HypreSolver_test > OUTPUT_Hypre
echo "done"
echo " " >> validation.log
echo "Hypre tests" >> validation.log
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
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln0.dat >> validation.log
echo "Solver/preconditioner combination 1" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln1.dat >> validation.log
echo "Solver/preconditioner combination 2" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln2.dat >> validation.log
echo "Solver/preconditioner combination 3" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln3.dat >> validation.log
echo "Solver/preconditioner combination 4" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln4.dat >> validation.log
echo "Solver/preconditioner combination 5" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln5.dat >> validation.log
echo "Solver/preconditioner combination 6" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln6.dat >> validation.log
echo "Solver/preconditioner combination 7" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln7.dat >> validation.log
echo "Solver/preconditioner combination 8" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln8.dat >> validation.log
echo "Solver/preconditioner combination 9" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln9.dat >> validation.log
echo "Solver/preconditioner combination 10" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln10.dat >> validation.log
echo "Solver/preconditioner combination 11" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln11.dat >> validation.log
echo "Solver/preconditioner combination 12" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln12.dat >> validation.log
echo "Solver/preconditioner combination 13" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln13.dat >> validation.log
echo "Solver/preconditioner combination 14" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln14.dat >> validation.log
echo "Solver/preconditioner combination 15" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln15.dat >> validation.log
echo "Solver/preconditioner combination 16" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln16.dat >> validation.log
echo "Solver/preconditioner combination 17" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln17.dat >> validation.log
echo "Solver/preconditioner combination 18" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln18.dat >> validation.log
echo "Solver/preconditioner combination 19" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln19.dat >> validation.log
echo "Solver/preconditioner combination 20" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln20.dat >> validation.log
echo "Solver/preconditioner combination 21" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_result.dat.gz  \
         RESLT/soln21.dat >> validation.log

echo "Number of Newton iterations for Hypre solves" >> validation.log
../../../bin/fpdiff.py ../validata/Hypre_conv.dat.gz  \
         RESLT/conv.dat >> validation.log
fi


mv RESLT RESLT_Hypre

else 

echo "Not running Hypre tests as executable doesn't exist"
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log
echo "[OK] (Dummy for non-existent Hypre)"  >> validation.log

fi



# Direct Solver Tests
#====================


# Validation for direct solver tests
#-----------------------------------
echo "Running direct solver tests"
mkdir RESLT
../direct_solver_test > OUTPUT_DirectSolver
echo "done"
echo " " >> validation.log
echo "Direct solver tests" >> validation.log
echo "-------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "DenseLU matrix based solve w/ DenseDoubleMatrix" >> validation.log
../../../bin/fpdiff.py ../validata/direct_solver_matrix_solve_result.dat.gz  \
         RESLT/DenseLU_DenseDoubleMatrix.dat >> validation.log
echo "DenseLU problem based solve" >> validation.log
../../../bin/fpdiff.py ../validata/direct_solver_problem_solve_result.dat.gz  \
         RESLT/soln0.dat >> validation.log
echo "FD_LU problem based solve" >> validation.log
../../../bin/fpdiff.py ../validata/direct_solver_problem_solve_result.dat.gz  \
         RESLT/soln1.dat >> validation.log
echo "SuperLU matrix based solve w/ global CRDoubleMatrix" >> validation.log
../../../bin/fpdiff.py ../validata/direct_solver_matrix_solve_result.dat.gz  \
         RESLT/SuperLU_CRDoubleMatrix.dat >> validation.log
echo "SuperLU matrix based solve w/ global CCDoubleMatrix" >> validation.log
../../../bin/fpdiff.py ../validata/direct_solver_matrix_solve_result.dat.gz  \
         RESLT/SuperLU_CCDoubleMatrix.dat >> validation.log
echo "SuperLU problem based solve" >> validation.log
../../../bin/fpdiff.py ../validata/direct_solver_problem_solve_result.dat.gz  \
         RESLT/soln2.dat >> validation.log
fi

mv RESLT RESLT_DirectSolverTest

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
