#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=85

# Threshold for number of iterations in comparison of convergence histories
#===========================================================================
threshold_for_number_of_iterations=3

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation


# Validation for simple block precond for linear elasticity
#----------------------------------------------------------

echo "Simple block preconditioner for 2D linear elasticity "
mkdir RESLT
../two_d_linear_elasticity_with_simple_block_diagonal_preconditioner > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "2D Linear elasticity simple preconditioner validation" >> validation.log
echo "-----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln.dat > linear_elasticity_simple_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > lin_elast_simple_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/linear_elasticity_simple_prec_results.dat.gz  \
        linear_elasticity_simple_prec_results.dat  1.0e-12 0.1 >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/compare_file_length_with_tolerance.bash \
        lin_elast_simple_iterative_solver_convergence.dat \
        ../validata/lin_elast_simple_iterative_solver_convergence.dat \
        $threshold_for_number_of_iterations \
        >>  validation.log
fi


mv RESLT RESLT_linear_elasticity

# Validation for simple block precond
#------------------------------------

echo "Simple block preconditioner for 'multi-poisson' "
mkdir RESLT
../two_d_multi_poisson --simple > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' simple preconditioner validation" >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0.dat > multi_poisson_simple_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > simple_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/multi_poisson_simple_prec_results.dat.gz  \
        multi_poisson_simple_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/compare_file_length_with_tolerance.bash \
        simple_iterative_solver_convergence.dat \
        ../validata/simple_iterative_solver_convergence.dat \
        $threshold_for_number_of_iterations \
        >>  validation.log
fi


mv RESLT RESLT_multi_poisson_simple_prec



# Validation for two plus three block precond
#--------------------------------------------

echo "Two plus three block preconditioner for 'multi-poisson' "
mkdir RESLT
../two_d_multi_poisson --two_plus_three > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' one plus four preconditioner validation" >> validation.log
echo "-------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0.dat > multi_poisson_two_plus_three_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > two_plus_three_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/multi_poisson_two_plus_three_prec_results.dat.gz  \
        multi_poisson_two_plus_three_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/compare_file_length_with_tolerance.bash \
        two_plus_three_iterative_solver_convergence.dat \
        ../validata/two_plus_three_iterative_solver_convergence.dat \
        $threshold_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_two_plus_three_prec


# Validation for two plus three upper triangular precond
#-------------------------------------------------------

echo "Two plus three upper triangular block preconditioner for 'multi-poisson' "
mkdir RESLT
../two_d_multi_poisson --two_plus_three_upper_triangular > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' two plus three upper triangular validation" >> validation.log
echo "------------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0.dat > multi_poisson_two_plus_three_upper_triangular_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > two_plus_three_upper_triangular_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/multi_poisson_two_plus_three_upper_triangular_prec_results.dat.gz  \
        multi_poisson_two_plus_three_upper_triangular_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/compare_file_length_with_tolerance.bash \
        two_plus_three_upper_triangular_iterative_solver_convergence.dat \
        ../validata/two_plus_three_upper_triangular_iterative_solver_convergence.dat \
        $threshold_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_two_plus_three_upper_triangular_prec


# Validation for two plus three upper triangular with sub block precond
#----------------------------------------------------------------------

echo "Two plus three upper triangular block preconditioner with sub for 'multi-poisson' "
mkdir RESLT
../two_d_multi_poisson --two_plus_three_upper_triangular_with_sub > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' two pus three upper triangular with sub validation" >> validation.log
echo "-----------------------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0.dat > multi_poisson_two_plus_three_upper_triangular_with_sub_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > two_plus_three_upper_triangular_with_sub_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/multi_poisson_two_plus_three_upper_triangular_with_sub_prec_results.dat.gz  \
        multi_poisson_two_plus_three_upper_triangular_with_sub_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/compare_file_length_with_tolerance.bash \
        two_plus_three_upper_triangular_with_sub_iterative_solver_convergence.dat \
        ../validata/two_plus_three_upper_triangular_with_sub_iterative_solver_convergence.dat \
        $threshold_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_two_plus_three_upper_triangular_with_sub_subsidiary_prec

# Validation for two plus three upper triangular with two sub block precond
#--------------------------------------------------------------------------

echo "Two plus three upper triangular block preconditioner with two sub for 'multi-poisson' "
mkdir RESLT
../two_d_multi_poisson --two_plus_three_upper_triangular_with_two_sub > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' two pus three upper triangular with two sub validation" >> validation.log
echo "-----------------------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0.dat > multi_poisson_two_plus_three_upper_triangular_with_two_sub_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > two_plus_three_upper_triangular_with_two_sub_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/multi_poisson_two_plus_three_upper_triangular_with_two_sub_prec_results.dat.gz  \
        multi_poisson_two_plus_three_upper_triangular_with_two_sub_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/compare_file_length_with_tolerance.bash \
        two_plus_three_upper_triangular_with_two_sub_iterative_solver_convergence.dat \
        ../validata/two_plus_three_upper_triangular_with_two_sub_iterative_solver_convergence.dat \
        $threshold_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_two_plus_three_upper_triangular_with_two_sub_subsidiary_prec

# Validation for two plus three upper triangular with replace block precond
#--------------------------------------------------------------------------

echo "Two plus three upper triangular block preconditioner with replace for 'multi-poisson' "
mkdir RESLT
../two_d_multi_poisson --two_plus_three_upper_triangular_with_replace > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' two pus three upper triangular with replace validation" >> validation.log
echo "-----------------------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0.dat > multi_poisson_two_plus_three_upper_triangular_with_replace_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > two_plus_three_upper_triangular_with_replace_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/multi_poisson_two_plus_three_upper_triangular_with_replace_prec_results.dat.gz  \
        multi_poisson_two_plus_three_upper_triangular_with_replace_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/compare_file_length_with_tolerance.bash \
        two_plus_three_upper_triangular_with_replace_iterative_solver_convergence.dat \
        ../validata/two_plus_three_upper_triangular_with_replace_iterative_solver_convergence.dat \
        $threshold_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_two_plus_three_upper_triangular_with_replace_subsidiary_prec

# Validation for coarse two plus two plus two one
#------------------------------------------------

echo "Coarse two plus two plus one block preconditioner for 'multi-poisson' "
mkdir RESLT
../two_d_multi_poisson --coarse_two_plus_two_plus_one > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' coarse two plus two plus one with subsidary and replace preconditioner validation" >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0.dat > multi_poisson_coarse_two_plus_two_plus_one_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > coarse_two_plus_two_plus_one_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/multi_poisson_coarse_two_plus_two_plus_one_prec_results.dat.gz  \
        multi_poisson_coarse_two_plus_two_plus_one_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/compare_file_length_with_tolerance.bash \
        coarse_two_plus_two_plus_one_iterative_solver_convergence.dat \
        ../validata/coarse_two_plus_two_plus_one_iterative_solver_convergence.dat \
        $threshold_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_coarse_two_plus_two_plus_one_prec


# Validation for two level coarsening
#------------------------------------------------

echo "Coarse one plus four with two level coarsening block preconditioner for 'multi-poisson' "
mkdir RESLT
../two_d_multi_poisson --one_plus_four_with_two_coarse > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' two level coarsening preconditioner validation" >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0.dat > multi_poisson_one_plus_four_with_two_coarse_results.dat
cat  RESLT/iterative_solver_convergence.dat > one_plus_four_with_two_coarse_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/multi_poisson_one_plus_four_with_two_coarse_results.dat.gz  \
        multi_poisson_one_plus_four_with_two_coarse_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/compare_file_length_with_tolerance.bash \
        one_plus_four_with_two_coarse_iterative_solver_convergence.dat \
        ../validata/one_plus_four_with_two_coarse_iterative_solver_convergence.dat \
        $threshold_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_one_plus_four_with_two_coarse_prec



# Validation for upper triangular block precond
#--------------------------------------------------------

echo "Upper diag block preconditioner for 'multi-poisson' "
mkdir RESLT
../two_d_multi_poisson --upper_triangular > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' upper_triangular preconditioner validation" >> validation.log
echo "----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0.dat > multi_poisson_upper_triangular_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > upper_triangular_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/multi_poisson_upper_triangular_prec_results.dat.gz  \
        multi_poisson_upper_triangular_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/compare_file_length_with_tolerance.bash \
        upper_triangular_iterative_solver_convergence.dat \
        ../validata/upper_triangular_iterative_solver_convergence.dat \
        $threshold_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_upper_triangular_prec

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
    ../../../bin/compare_file_length_with_tolerance.bash GMRES_convergence.dat \
    ../validata/GMRES_convergence.dat $threshold_for_number_of_iterations \
    >>  validation.log

    #Compare number of iterations against reference data and append
    ../../../bin/compare_file_length_with_tolerance.bash GMRES_la_solve_convergence.dat \
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

echo ""
echo "Not running Trilinos tests as executable doesn't exist"
echo "TRILINOS can't handle serial execution of MPI-enabled code"
echo ""

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

echo ""
echo "Not running Hypre tests as executable doesn't exist"
echo "HYPRE can't handle serial execution of MPI-enabled code"
echo ""

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
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
