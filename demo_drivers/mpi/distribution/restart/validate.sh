#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=1


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation



# Validation for load-balanced doubly adaptive unsteady heat
#-----------------------------------------------------------
echo "Running distributed load-balanced doubly adaptive 2D unsteady heat validation "
mkdir RESLT
mkdir RESLT_after_load_balance
$MPI_RUN_COMMAND ../two_d_unsteady_heat_2adapt_load_balance --validation_run --partitioning_file doubly_adaptive_partitioning_load_balance.dat > OUTPUT_doubly_adaptive_load_balanced_for_restart
echo "done run for restart"
echo " " >> validation.log
echo "2D distributed load-balanced doubly adaptive unsteady heat validation " >> validation.log
echo "---------------------------------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

sleep 5
cat RESLT/soln0_on_proc0.dat  RESLT/soln0_on_proc1.dat  \
    RESLT/soln9_on_proc0.dat  RESLT/soln9_on_proc1.dat  \
    > result_doubly_adaptive_load_balanced.dat
mv RESLT RESLT_doubly_adaptive_load_balanced_for_restart

mkdir RESLT
$MPI_RUN_COMMAND ../two_d_unsteady_heat_2adapt_load_balance --validation_run --restart_file RESLT_doubly_adaptive_load_balanced_for_restart/restart3 --partitioning_file RESLT_doubly_adaptive_load_balanced_for_restart/partitioning.dat > OUTPUT_doubly_adaptive_load_balanced_restarted
echo "done restarted run"

sleep 5
cat RESLT/soln9_on_proc0.dat  RESLT/soln9_on_proc1.dat  \
    >> result_doubly_adaptive_load_balanced.dat
mv RESLT RESLT_doubly_adaptive_load_balanced_restarted

sleep 5
mkdir RESLT
sleep 20
$MPI_RUN_COMMAND ../two_d_unsteady_heat_2adapt_load_balance --validation_run --restart_file RESLT_doubly_adaptive_load_balanced_for_restart/restart7 --partitioning_file RESLT_doubly_adaptive_load_balanced_for_restart/load_balanced_partitioning.dat > OUTPUT_doubly_adaptive_load_balanced_restarted_from_load_balanced
echo "done restarted run from load balanced solution"

sleep 5
cat RESLT/soln9_on_proc0.dat  RESLT/soln9_on_proc1.dat  \
    metis_input_for_validation.dat \
    >> result_doubly_adaptive_load_balanced.dat
mv RESLT RESLT_doubly_adaptive_load_balanced_restarted_from_load_balanced



if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/result_doubly_adaptive_load_balanced.dat.gz \
    result_doubly_adaptive_load_balanced.dat  >> validation.log
fi







# Append output to global validation log file
#--------------------------------------------
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
