#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=6


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo poisson
#----------------------------
cd Validation

echo "Running Boussinesq Convection validation with FD for off-diagonals"
mkdir RESLT
../boussinesq_convection lalala > OUTPUT_boussinesq_convection
echo "done"
echo " " >> validation.log
echo "Boussinesq Convection validation with FD for off-diagonals" >> validation.log
echo "-----------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln5.dat > bous_convection_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz   \
    bous_convection_results.dat 0.1 1.0e-8  >> validation.log
fi

mv RESLT RESLT_non_refineable



echo "Running Boussinesq Convection validation with FD for entire Jacobian"
mkdir RESLT
../boussinesq_convection_fd lalala > OUTPUT_boussinesq_convection_FD
echo "done"
echo " " >> validation.log
echo "Boussinesq Convection validation with FD for entire Jacobian" >> validation.log
echo "------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln5.dat > bous_convection_fd_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz   \
    bous_convection_fd_results.dat 0.1 2.0e-7  >> validation.log
fi

mv RESLT RESLT_non_refineable_fd


echo "Running Boussinesq Convection validation with analytic Jacobian"
mkdir RESLT
../boussinesq_convection_analytic lalala > OUTPUT_boussinesq_convection_analytic
echo "done"
echo " " >> validation.log
echo "Boussinesq Convection validation with analytic Jacobian" >> validation.log
echo "------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln5.dat > bous_convection_anal_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz   \
    bous_convection_anal_results.dat 0.1 5.0e-7  >> validation.log
fi

mv RESLT RESLT_non_refineable_anal


echo "Running Refineable Boussinesq Convection validation "
mkdir RESLT
../refineable_b_convection > OUTPUT_refineable_b_convection
echo "done"
echo " " >> validation.log
echo "Refineable Boussinesq Convection validation" >> validation.log
echo "-------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat > ref_bous_convection_results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/ref_results.dat.gz   \
   ref_bous_convection_results.dat  0.1 1.0e-7 >> validation.log
fi

mv RESLT RESLT_refineable


# Validation for Boussinesq convection problem using multi-domain method
#-----------------------------------------------------------------------

echo "Running Boussinesq convection problem (multi-domain method) "
mkdir RESLT
../multimesh_boussinesq_convection validate > OUTPUT_multimesh_b_convection
echo "done"
echo " " >> validation.log
echo "Boussinesq convection (multi-domain) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln5.dat \
    > multimesh_b_convection.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/multimesh_b_convection.dat.gz  \
         multimesh_b_convection.dat 0.1 1.5e-7 >> validation.log
fi

mv RESLT RESLT_multimesh_boussinesq_convection


# Validation for refineable Boussinesq convection problem, multi-domain
#----------------------------------------------------------------------

echo "Running refineable Boussinesq convection problem (multi-domain method) "
mkdir RESLT
../multimesh_ref_b_convection > OUTPUT_multimesh_ref_b_convection
echo "done"
echo " " >> validation.log
echo "Refineable Boussinesq convection (multi-domain) validation" >> validation.log
echo "------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat \
    > multimesh_ref_b_convection.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/multimesh_ref_b_convection.dat.gz  \
         multimesh_ref_b_convection.dat 0.1 1.0e-7 >> validation.log
fi

mv RESLT RESLT_multimesh_ref_b_convection


# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log


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
