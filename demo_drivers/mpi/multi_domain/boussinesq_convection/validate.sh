#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=4

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

# Validation for Boussinesq convection problem using multi-domain method
#-----------------------------------------------------------------------

echo "Running Boussinesq convection problem (multi-domain method, analytic) "
mkdir RESLT
$MPI_RUN_COMMAND ../multimesh_boussinesq_convection validate > OUTPUT_multimesh_boussinesq_convection
echo "done"
echo " " >> validation.log
echo "Boussinesq convection (multi-domain, analytic) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0_on_proc0.dat RESLT/soln1_on_proc1.dat RESLT/soln2_on_proc0.dat \
    RESLT/soln3_on_proc1.dat RESLT/soln4_on_proc0.dat RESLT/soln5_on_proc1.dat \
    > multimesh_boussinesq_convection_analytic_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/multimesh_boussinesq_convection_analytic_results.dat.gz  \
         multimesh_boussinesq_convection_analytic_results.dat 0.1 1.5e-7 >> validation.log
fi

mkdir RESLT_multimesh_boussinesq_convection_analytic
mv RESLT/* RESLT_multimesh_boussinesq_convection_analytic

# Validation for Boussinesq convection problem using multi-domain method (FD)
#----------------------------------------------------------------------------

echo "Running Boussinesq convection problem (multi-domain method, FD) "
$MPI_RUN_COMMAND ../multimesh_boussinesq_convection_fd validate > OUTPUT_multimesh_boussinesq_convection_fd
echo "done"
echo " " >> validation.log
echo "Boussinesq convection (multi-domain, FD) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0_on_proc0.dat RESLT/soln1_on_proc1.dat RESLT/soln2_on_proc0.dat \
    RESLT/soln3_on_proc1.dat RESLT/soln4_on_proc0.dat RESLT/soln5_on_proc1.dat \
    > multimesh_boussinesq_convection_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/multimesh_boussinesq_convection_results.dat.gz  \
         multimesh_boussinesq_convection_results.dat 0.1 4e-7 >> validation.log
fi

mkdir RESLT_multimesh_boussinesq_convection_fd
mv RESLT/* RESLT_multimesh_boussinesq_convection_fd

#-----------------------------------------------------------------------

# Validation for refineable Boussinesq convection problem, single domain
#-----------------------------------------------------------------------

echo "Running refineable Boussinesq convection problem (dist. single domain) "
$MPI_RUN_COMMAND ../refineable_b_convection > OUTPUT_refineable_b_convection
echo "done"
echo " " >> validation.log
echo "Refineable Boussinesq convection (single domain) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0_on_proc0.dat RESLT/soln0_on_proc1.dat \
    RESLT/soln1_on_proc0.dat RESLT/soln1_on_proc1.dat \
    > refineable_b_convection_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/refineable_b_convection_results.dat.gz  \
         refineable_b_convection_results.dat 0.1 1.0e-7 >> validation.log
fi

mkdir RESLT_refineable_b_convection
mv RESLT/* RESLT_refineable_b_convection

#----------------------------------------------------------------------

# Validation for refineable Boussinesq convection problem, multi-domain
#----------------------------------------------------------------------

echo "Running refineable Boussinesq convection problem (multi-domain method) "
$MPI_RUN_COMMAND ../multimesh_ind_ref_b_convection > OUTPUT_multimesh_ind_ref_b_convection
echo "done"
echo " " >> validation.log
echo "Refineable Boussinesq convection (multi-domain) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0_on_proc0.dat RESLT/soln0_on_proc1.dat \
    RESLT/soln1_on_proc0.dat RESLT/soln1_on_proc1.dat \
    > multimesh_ind_ref_b_convection_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/multimesh_ind_ref_b_convection_results.dat.gz  \
         multimesh_ind_ref_b_convection_results.dat 0.1 1.0e-7 >> validation.log
fi

mv RESLT RESLT_multimesh_ind_ref_b_convection


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
