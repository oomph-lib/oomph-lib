#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=12

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

# Validation for 2D mesh distribution test
#-----------------------------------------

echo "Running 2D mesh distribution problem "
mkdir RESLT
$MPI_RUN_COMMAND ../two_d_mesh_dist > OUTPUT_two_d_mesh_dist
echo "done"
echo " " >> validation.log
echo "2D mesh distribution validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/mesh0_0.dat RESLT/mesh1_0.dat RESLT/haloed_nodes_on_proc0_0.dat \
    RESLT/haloed_nodes_on_proc1_0.dat RESLT/halo_elements_on_proc0_0.dat \
    RESLT/halo_elements_on_proc1_0.dat > two_d_mesh_dist_results.dat
cat RESLT/soln0_on_proc0.dat RESLT/soln0_on_proc1.dat RESLT/soln1_on_proc0.dat \
    RESLT/soln1_on_proc1.dat RESLT/soln2_on_proc0.dat RESLT/soln2_on_proc1.dat \
    > two_d_soln_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/two_d_mesh_dist_results.dat.gz  \
         two_d_mesh_dist_results.dat >> validation.log
../../../../../bin/fpdiff.py ../validata/two_d_soln_results.dat.gz  \
         two_d_soln_results.dat >> validation.log
fi

mv RESLT RESLT_two_d_mesh_dist

#----------------------------------------------------------------------

# Validation for fish poisson
#----------------------------

echo "Running fish poisson validation "
mkdir RESLT_select_refine
mkdir RESLT_incremental2
mkdir RESLT_fully_automatic
$MPI_RUN_COMMAND ../fish_poisson > OUTPUT_fish_poisson
echo "done"
echo " " >> validation.log
echo "Fish poisson validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_select_refine/soln0_on_proc0.dat RESLT_select_refine/soln0_on_proc1.dat \
    RESLT_select_refine/soln1_on_proc0.dat RESLT_select_refine/soln1_on_proc1.dat \
    RESLT_select_refine/soln2_on_proc0.dat RESLT_select_refine/soln2_on_proc1.dat \
    > fish_poisson_select_results.dat
cat RESLT_incremental2/soln1_on_proc0.dat RESLT_incremental2/soln1_on_proc1.dat \
    RESLT_incremental2/soln6_on_proc0.dat RESLT_incremental2/soln6_on_proc1.dat \
    RESLT_incremental2/soln16_on_proc0.dat RESLT_incremental2/soln16_on_proc1.dat \
    > fish_poisson_incremental_results.dat
cat RESLT_fully_automatic/soln0_on_proc0.dat RESLT_fully_automatic/soln0_on_proc1.dat \
    > fish_poisson_automatic_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/fish_poisson_select_results.dat.gz  \
         fish_poisson_select_results.dat >> validation.log
../../../../../bin/fpdiff.py ../validata/fish_poisson_incremental_results.dat.gz  \
         fish_poisson_incremental_results.dat >> validation.log
../../../../../bin/fpdiff.py ../validata/fish_poisson_automatic_results.dat.gz  \
         fish_poisson_automatic_results.dat >> validation.log
fi

mkdir RESLT_fish_poisson
mv RESLT_select_refine RESLT_fish_poisson
mv RESLT_incremental2 RESLT_fish_poisson
mv RESLT_fully_automatic RESLT_fish_poisson

#----------------------------------------------------------------------

# Validation for adaptive driven cavity (TH & CR)
#------------------------------------------------

echo "Running adaptive rectangular driven cavity (TH & CR) validation "
mkdir RESLT_TH
mkdir RESLT_CR
$MPI_RUN_COMMAND ../adaptive_driven_cavity > OUTPUT_adaptive_driven_cavity
echo "done"
echo " " >> validation.log
echo "Adaptive rectangular driven cavity (TH & CR) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_TH/soln0_on_proc0.dat RESLT_TH/soln0_on_proc1.dat \
    RESLT_TH/soln1_on_proc0.dat RESLT_TH/soln1_on_proc1.dat \
    > adaptive_cavity_TH_results.dat
cat RESLT_CR/soln0_on_proc0.dat RESLT_CR/soln0_on_proc1.dat \
    RESLT_CR/soln1_on_proc0.dat RESLT_CR/soln1_on_proc1.dat \
    > adaptive_cavity_CR_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/adaptive_cavity_TH_results.dat.gz  \
         adaptive_cavity_TH_results.dat >> validation.log
../../../../../bin/fpdiff.py ../validata/adaptive_cavity_CR_results.dat.gz  \
         adaptive_cavity_CR_results.dat >> validation.log
fi

mkdir RESLT_adaptive_driven_cavity
mv RESLT_TH RESLT_adaptive_driven_cavity
mv RESLT_CR RESLT_adaptive_driven_cavity

#----------------------------------------------------------------------

# Validation for circular driven cavity (TH, 2 problems)
#-------------------------------------------------------

echo "Running adaptive circular driven cavity (TH, 2 problems) validation "
mkdir RESLT
$MPI_RUN_COMMAND ../circular_driven_cavity > OUTPUT_circular_driven_cavity
echo "done"
echo " " >> validation.log
echo "Adaptive circular driven cavity validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln1_on_proc0.dat RESLT/soln1_on_proc1.dat > circular_cavity_TH1_results.dat
cat RESLT/soln2_on_proc0.dat RESLT/soln2_on_proc1.dat > circular_cavity_TH2_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/circular_cavity_TH1_results.dat.gz  \
         circular_cavity_TH1_results.dat >> validation.log
../../../../../bin/fpdiff.py ../validata/circular_cavity_TH2_results.dat.gz  \
         circular_cavity_TH2_results.dat >> validation.log
fi

mv RESLT RESLT_circular_driven_cavity

#----------------------------------------------------------------------

# Validation for airy cantilever (2D solid) deformation
#------------------------------------------------------

echo "Running airy cantilever validation "
mkdir RESLT_refine0
mkdir RESLT_refine1
mkdir RESLT_refine2
mkdir RESLT_refine3
mkdir RESLT_refine4
mkdir RESLT_refine5
mkdir RESLT_refine6
mkdir RESLT_refine7
mkdir RESLT_refine8
mkdir RESLT_refine9
$MPI_RUN_COMMAND ../airy_cantilever2_adapt > OUTPUT_airy_cantilever2_adapt
echo "done"
echo " " >> validation.log
echo "Airy cantilever validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_refine0/soln0_on_proc0.dat RESLT_refine1/soln0_on_proc1.dat \
    RESLT_refine2/soln1_on_proc0.dat RESLT_refine3/soln1_on_proc1.dat \
    RESLT_refine4/soln0_on_proc0.dat RESLT_refine5/soln0_on_proc1.dat \
    RESLT_refine6/soln1_on_proc0.dat RESLT_refine7/soln1_on_proc1.dat \
    RESLT_refine8/soln0_on_proc0.dat RESLT_refine8/soln0_on_proc1.dat \
    RESLT_refine9/soln1_on_proc0.dat RESLT_refine9/soln1_on_proc1.dat \
    > airy_cantilever_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/airy_cantilever_results.dat.gz  \
         airy_cantilever_results.dat 0.1 1.0e-7 >> validation.log
fi

mkdir RESLT_airy_cantilever2_adapt
mv RESLT_refine* RESLT_airy_cantilever2_adapt

#----------------------------------------------------------------------

# Validation for 2D poisson with flux b.c.'s
#-----------------------------------------------

echo "Running 2D poisson problem with flux boundary conditions "
mkdir RESLT
$MPI_RUN_COMMAND ../two_d_poisson_flux_bc_adapt > OUTPUT_two_d_poisson_flux_bc_adapt
echo "done"
echo " " >> validation.log
echo "2D poisson with flux b.c.'s validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0_on_proc0.dat RESLT/soln1_on_proc1.dat \
    RESLT/soln2_on_proc0.dat RESLT/soln3_on_proc1.dat \
    > two_d_poisson_flux_bc_adapt_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/two_d_poisson_flux_bc_adapt_results.dat.gz  \
         two_d_poisson_flux_bc_adapt_results.dat >> validation.log
fi

mkdir RESLT_two_d_poisson_flux_bc_adapt
mv RESLT/*.dat RESLT_two_d_poisson_flux_bc_adapt

#------------------------------------------------------

# Validation for prescribed displacement problem with Lagrange multipliers
#-------------------------------------------------------------------------

echo "Running prescribed displ. problem with Lagrange multipliers "
$MPI_RUN_COMMAND ../prescribed_displ_lagr_mult > OUTPUT_prescribed_displ_lagr_mult
echo "done"
echo " " >> validation.log
echo "Prescribed displ. with Lagrange multipliers validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0_on_proc0.dat RESLT/soln1_on_proc1.dat RESLT/soln2_on_proc0.dat \
    RESLT/lagr0_on_proc1.dat RESLT/lagr1_on_proc0.dat RESLT/lagr2_on_proc1.dat \
    > prescribed_displ_lagr_mult_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/prescribed_displ_lagr_mult_results.dat.gz  \
         prescribed_displ_lagr_mult_results.dat >> validation.log
fi

mv RESLT RESLT_prescribed_displ_lagr_mult

# Append log to main validation log
cat validation.log >> ../../../../../validation.log

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
