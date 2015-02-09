#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=46

# Threshold for number of iterations in comparison of convergence histories
# for multi_poisson
# This is used to compare the number of iterations using fpdiff if it is 
# given as a whole number. Typical number of iterations for 
# two_d_multi_poisson is 10, so allow for up to a difference of 3.
relative_tol_for_number_of_iterations=30

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
cp ../*partition.dat .


# Validation for simple block precond
#------------------------------------
if [ -f ../two_d_multi_poisson ]; then

echo "Simple block preconditioner for 'multi-poisson' "
mkdir RESLT
$MPI_RUN_COMMAND ../two_d_multi_poisson --use_validation_distribution --simple > \
    RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' simple preconditioner validation" >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0_on_proc0.dat RESLT/soln0_on_proc1.dat \
    > multi_poisson_simple_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > \
    simple_iterative_solver_convergence.dat
if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../../../../bin/fpdiff.py \
        ../validata/multi_poisson_simple_prec_results.dat.gz  \
        multi_poisson_simple_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../../../../bin/fpdiff.py \
        simple_iterative_solver_convergence.dat  \
        ../validata/simple_iterative_solver_convergence.dat \
        $relative_tol_for_number_of_iterations \
        >> validation.log
fi

mv RESLT RESLT_multi_poisson_simple_prec


# Validation for two plus three block precond
#--------------------------------------------

echo "Two plus three block preconditioner for 'multi-poisson' "
mkdir RESLT
$MPI_RUN_COMMAND ../two_d_multi_poisson --use_validation_distribution \
    --two_plus_three > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' one plus four preconditioner validation" >> validation.log
echo "-------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0_on_proc0.dat RESLT/soln0_on_proc1.dat \
    > multi_poisson_two_plus_three_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > \
    two_plus_three_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py \
        ../validata/multi_poisson_two_plus_three_prec_results.dat.gz  \
        multi_poisson_two_plus_three_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py \
        two_plus_three_iterative_solver_convergence.dat \
        ../validata/two_plus_three_iterative_solver_convergence.dat \
        $relative_tol_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_two_plus_three_prec


# Validation for two plus three upper triangular precond
#-------------------------------------------------------

echo "Two plus three upper triangular block preconditioner for 'multi-poisson' "
mkdir RESLT
$MPI_RUN_COMMAND ../two_d_multi_poisson --use_validation_distribution \
    --two_plus_three_upper_triangular > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' two plus three upper triangular validation" >> \
    validation.log
echo "------------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0_on_proc0.dat RESLT/soln0_on_proc1.dat \
    > multi_poisson_two_plus_three_upper_triangular_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > \
    two_plus_three_upper_triangular_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py \
        ../validata/multi_poisson_two_plus_three_upper_triangular_prec_results.dat.gz  \
        multi_poisson_two_plus_three_upper_triangular_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py \
        two_plus_three_upper_triangular_iterative_solver_convergence.dat \
        ../validata/two_plus_three_upper_triangular_iterative_solver_convergence.dat \
        $relative_tol_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_two_plus_three_upper_triangular_prec


# Validation for two plus three upper triangular with sub block precond
#----------------------------------------------------------------------

echo "Two plus three upper triangular block preconditioner with sub for 'multi-poisson' "
mkdir RESLT
$MPI_RUN_COMMAND ../two_d_multi_poisson --use_validation_distribution \
    --two_plus_three_upper_triangular_with_sub > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' two pus three upper triangular with sub validation" >> \
    validation.log
echo "-----------------------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0_on_proc0.dat RESLT/soln0_on_proc1.dat \
    > multi_poisson_two_plus_three_upper_triangular_with_sub_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > \
    two_plus_three_upper_triangular_with_sub_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py \
        ../validata/multi_poisson_two_plus_three_upper_triangular_with_sub_prec_results.dat.gz  \
        multi_poisson_two_plus_three_upper_triangular_with_sub_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py \
        two_plus_three_upper_triangular_with_sub_iterative_solver_convergence.dat \
        ../validata/two_plus_three_upper_triangular_with_sub_iterative_solver_convergence.dat \
        $relative_tol_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_two_plus_three_upper_triangular_with_sub_subsidiary_prec

# Validation for two plus three upper triangular with two sub block precond
#--------------------------------------------------------------------------

echo "Two plus three upper triangular block preconditioner with two sub for 'multi-poisson' "
mkdir RESLT
$MPI_RUN_COMMAND ../two_d_multi_poisson --use_validation_distribution --two_plus_three_upper_triangular_with_two_sub > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' two pus three upper triangular with two sub validation" >> validation.log
echo "-----------------------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0_on_proc0.dat RESLT/soln0_on_proc1.dat \
    > multi_poisson_two_plus_three_upper_triangular_with_two_sub_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > \
    two_plus_three_upper_triangular_with_two_sub_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/multi_poisson_two_plus_three_upper_triangular_with_two_sub_prec_results.dat.gz  \
        multi_poisson_two_plus_three_upper_triangular_with_two_sub_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py \
        two_plus_three_upper_triangular_with_two_sub_iterative_solver_convergence.dat \
        ../validata/two_plus_three_upper_triangular_with_two_sub_iterative_solver_convergence.dat \
        $relative_tol_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_two_plus_three_upper_triangular_with_two_sub_subsidiary_prec

# Validation for two plus three upper triangular with replace block precond
#--------------------------------------------------------------------------

echo "Two plus three upper triangular block preconditioner with replace for 'multi-poisson' "
mkdir RESLT
$MPI_RUN_COMMAND ../two_d_multi_poisson --use_validation_distribution \
    --two_plus_three_upper_triangular_with_replace > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' two pus three upper triangular with replace validation" >> validation.log
echo "-----------------------------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0_on_proc0.dat RESLT/soln0_on_proc1.dat \
    > multi_poisson_two_plus_three_upper_triangular_with_replace_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > two_plus_three_upper_triangular_with_replace_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/multi_poisson_two_plus_three_upper_triangular_with_replace_prec_results.dat.gz  \
        multi_poisson_two_plus_three_upper_triangular_with_replace_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py \
        two_plus_three_upper_triangular_with_replace_iterative_solver_convergence.dat \
        ../validata/two_plus_three_upper_triangular_with_replace_iterative_solver_convergence.dat \
        $relative_tol_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_two_plus_three_upper_triangular_with_replace_subsidiary_prec

# Validation for coarse two plus two plus two one
#------------------------------------------------

echo "Coarse two plus two plus one block preconditioner for 'multi-poisson' "
mkdir RESLT
$MPI_RUN_COMMAND ../two_d_multi_poisson --use_validation_distribution --coarse_two_plus_two_plus_one > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' coarse two plus two plus one with subsidary and replace preconditioner validation" >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0_on_proc0.dat RESLT/soln0_on_proc1.dat \
    > multi_poisson_coarse_two_plus_two_plus_one_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > coarse_two_plus_two_plus_one_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/multi_poisson_coarse_two_plus_two_plus_one_prec_results.dat.gz  \
        multi_poisson_coarse_two_plus_two_plus_one_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py \
        coarse_two_plus_two_plus_one_iterative_solver_convergence.dat \
        ../validata/coarse_two_plus_two_plus_one_iterative_solver_convergence.dat \
        $relative_tol_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_coarse_two_plus_two_plus_one_prec


# Validation for one plus four with two level coarsening.
#------------------------------------------------

echo "One plus four with two level coarsening block preconditioner for 'multi-poisson' "
mkdir RESLT
$MPI_RUN_COMMAND ../two_d_multi_poisson --use_validation_distribution --one_plus_four_with_two_coarse > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' one plus four with two level coarsening preconditioner validation" >> validation.log
echo "------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0_on_proc0.dat RESLT/soln0_on_proc1.dat \
    > multi_poisson_one_plus_four_with_two_coarse_results.dat
cat  RESLT/iterative_solver_convergence.dat > one_plus_four_with_two_coarse_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/multi_poisson_one_plus_four_with_two_coarse_results.dat.gz  \
        multi_poisson_one_plus_four_with_two_coarse_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py \
        one_plus_four_with_two_coarse_iterative_solver_convergence.dat \
        ../validata/one_plus_four_with_two_coarse_iterative_solver_convergence.dat.gz \
        $relative_tol_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_one_plus_four_with_two_coarse_prec



# Validation for upper triangular block precond
#--------------------------------------------------------

echo "Upper diag block preconditioner for 'multi-poisson' "
mkdir RESLT
$MPI_RUN_COMMAND ../two_d_multi_poisson --use_validation_distribution --upper_triangular > RESLT/OUTPUT
echo "done"
echo " " >> validation.log
echo "'Multi-Poisson' upper_triangular preconditioner validation" >> validation.log
echo "----------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat  RESLT/soln0_on_proc0.dat RESLT/soln0_on_proc1.dat \
    > multi_poisson_upper_triangular_prec_results.dat
cat  RESLT/iterative_solver_convergence.dat > upper_triangular_iterative_solver_convergence.dat

if test "$1" = "no_fpdiff"; then
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
    # Compare results
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py ../validata/multi_poisson_upper_triangular_prec_results.dat.gz  \
        multi_poisson_upper_triangular_prec_results.dat >> validation.log
    
    #Compare number of iterations against reference data and append
    ../$OOMPH_ROOT_DIR/bin/fpdiff.py \
        upper_triangular_iterative_solver_convergence.dat \
        ../validata/upper_triangular_iterative_solver_convergence.dat \
        $relative_tol_for_number_of_iterations \
        >>  validation.log
fi

mv RESLT RESLT_multi_poisson_upper_triangular_prec


# if we don't have trilinos we can't run the multi poisson self tests...

else
    echo " " 
    echo "Cannot run multi poisson block preconditioners; requires trilinos"
    echo " " 
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
    echo "dummy [OK] -- for multi poisson block preconditioners; can only run with trilinos" >> validation.log
fi


# Validation for rectangular driven cavity
#-----------------------------------------

echo "Running rectangular driven cavity LSC precond validation "
mkdir RESLT

# Wait for a bit to allow parallel file systems to realise
# the existence of the new directory
sleep 5

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
../../../../bin/fpdiff.py ../validata/driven_cavity_results.dat.gz  \
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
../../../../bin/fpdiff.py ../validata/airy_soln.dat.gz \
    RESLT/airy_soln0.dat 0.1 10e-6 >> validation.log
fi
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/airy_soln.dat.gz \
    RESLT/airy_soln1.dat 0.1 10e-6 >> validation.log
fi
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/airy_soln.dat.gz \
    RESLT/airy_soln2.dat 0.1 10e-6 >> validation.log
fi
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/airy_soln.dat.gz \
    RESLT/airy_soln3.dat 0.1 10e-6 >> validation.log
fi

# FSI Preconditioner
#===================
echo "Running FSI Preconditioner on Channel with Leaflet Problem"
$MPI_RUN_COMMAND ../fsi_channel_with_leaflet > OUTPUT_fsi
echo "done"
echo ""
echo " " >> validation.log
echo "FSI Preconditioner Tests" >> validation.log
echo "--------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/fsi_fluid_soln.dat.gz \
    RESLT/fsi_fluid_soln0.dat 0.1 10e-8 >> validation.log
fi
if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/fsi_wall_soln.dat.gz \
    RESLT/fsi_wall_soln0.dat 0.1 10e-6 >> validation.log
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
../../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_0.dat >> validation.log
echo "Solver/preconditioner combination 1" >> validation.log
../../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_1.dat >> validation.log
echo "Solver/preconditioner combination 2" >> validation.log
../../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_2.dat >> validation.log
echo "Solver/preconditioner combination 3" >> validation.log
../../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_3.dat >> validation.log
echo "Solver/preconditioner combination 4" >> validation.log
../../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_4.dat >> validation.log
echo "Solver/preconditioner combination 5" >> validation.log
../../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_5.dat >> validation.log
echo "Solver/preconditioner combination 6" >> validation.log
../../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_6.dat >> validation.log
echo "Solver/preconditioner combination 7" >> validation.log
../../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_7.dat >> validation.log
echo "Solver/preconditioner combination 8" >> validation.log
../../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_8.dat >> validation.log
echo "Solver/preconditioner combination 9" >> validation.log
../../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_9.dat >> validation.log
echo "Solver/preconditioner combination 10" >> validation.log
../../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_10.dat >> validation.log
echo "Solver/preconditioner combination 11" >> validation.log
../../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_11.dat >> validation.log
echo "Solver/preconditioner combination 12" >> validation.log
../../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_12.dat >> validation.log
echo "Solver/preconditioner combination 13" >> validation.log
../../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_13.dat >> validation.log
echo "Solver/preconditioner combination 14" >> validation.log
../../../../bin/fpdiff.py ../validata/Trilinos_result.dat.gz  \
         RESLT/soln_trilinos_14.dat >> validation.log
echo "Number of Newton iterations for Trilinos solves" >> validation.log
../../../../bin/fpdiff.py ../validata/Trilinos_conv.dat.gz  \
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


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "SuperLU_dist matrix based solve w/ global CRDoubleMatrixSolver" >> validation.log
../../../../bin/fpdiff.py ../validata/direct_solver_matrix_solve_result.dat.gz  \
         RESLT/SuperLU_dist_CRDoubleMatrix_global.dat >> validation.log
echo "SuperLU_dist matrix based solve w/ dist CRDoubleMatrixSolver" >> validation.log
../../../../bin/fpdiff.py ../validata/direct_solver_matrix_solve_result.dat.gz  \
         RESLT/SuperLU_dist_CRDoubleMatrix_distributed.dat >> validation.log
echo "SuperLU_dist global problem based solve" >> validation.log
../../../../bin/fpdiff.py ../validata/direct_solver_problem_solve_result.dat.gz  \
         RESLT/soln_direct_solver_1.dat >> validation.log
echo "SuperLU_dist global problem based solve" >> validation.log
../../../../bin/fpdiff.py ../validata/direct_solver_problem_solve_result.dat.gz  \
         RESLT/soln_direct_solver_2.dat >> validation.log

echo "MUMPS-based global problem based solve" >> validation.log
if [ -f RESLT/dummy_mumps.dat ]; then
    echo "Using dummy data for MUMPS self-test (don't have mumps!)"  
    echo "[OK] (Dummy for non-existent MUMPS)"  >> validation.log
else
    ../../../../bin/fpdiff.py ../validata/direct_solver_problem_solve_result.dat.gz  \
        RESLT/soln_direct_solver_3.dat >> validation.log
fi
fi

mv RESLT RESLT_test

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
