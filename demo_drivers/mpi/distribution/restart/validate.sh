#! /bin/sh


#Set the number of tests to be checked
NUM_TESTS=12

validation_run_flag=--validation_run
#validation_run_flag=
echo $validation_run_flag


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation

# Loop over first pruning / first load balancing
#-----------------------------------------------
first_list="0 1"
for first in `echo $first_list`; do

    prune_first_flag=" " 
    first_flag="prune"
    if [ $first -eq "0" ]; then
        prune_first_flag=" --load_balance_first "
        first_flag="load_balance"
    fi

    
    # Validation for load-balanced doubly adaptive unsteady heat
    #-----------------------------------------------------------
    echo "Running distributed load-balanced doubly adaptive 2D unsteady heat validation; first: $first_flag "
    mkdir RESLT_`echo $first_flag`_first_for_restart
    $MPI_RUN_COMMAND ../two_d_unsteady_heat_2adapt_load_balance $validation_run_flag $prune_first_flag --partitioning_file ../partitioning.dat > OUTPUT_`echo $first_flag`_first_for_restart
    echo "done run for restart"
    echo " " >> validation.log
    echo "2D distributed load-balanced doubly adaptive unsteady heat validation; first: $first_flag " >> validation.log
    echo "---------------------------------------------------------------------------------" >> validation.log
    echo " " >> validation.log
    echo "Validation directory: " >> validation.log
    echo " " >> validation.log
    echo "  " `pwd` >> validation.log
    echo " " >> validation.log
    
    sleep 5
    cat RESLT_`echo $first_flag`_first_for_restart/soln21_on_proc0.dat  RESLT_`echo $first_flag`_first_for_restart/soln21_on_proc1.dat  \
        > result_`echo $first_flag`_first_for_restart.dat
    #mv RESLT RESLT_`echo $first_flag`_first_for_restart
    
    
    if test "$1" = "no_fpdiff"; then
        echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
    else
        ../../../../../bin/fpdiff.py ../validata/result_`echo $first_flag`_first_for_restart.dat.gz  result_`echo $first_flag`_first_for_restart.dat  >> validation.log
    fi

    restart_list="2 7 11 16 20"
    for restart_step in `echo $restart_list`; do
        
        mkdir RESLT_`echo $first_flag`_first_restarted_from_step_restart`echo $restart_step`
        $MPI_RUN_COMMAND ../two_d_unsteady_heat_2adapt_load_balance $validation_run_flag  $prune_first_flag --restart_file RESLT_`echo $first_flag`_first_for_restart/restart$restart_step --partitioning_file RESLT_`echo $first_flag`_first_for_restart/partitioning.dat  > OUTPUT_`echo $first_flag`_first_restarted_from_step_restart$restart_step
        echo "done restarted run from step $restart_step"
        
        sleep 5
        cat RESLT_`echo $first_flag`_first_restarted_from_step_restart`echo $restart_step`/soln21_on_proc0.dat  RESLT_`echo $first_flag`_first_restarted_from_step_restart`echo $restart_step`/soln21_on_proc1.dat  >> result_`echo $first_flag`_first_restarted_from_step_restart`echo $restart_step`.dat
        #mv RESLT RESLT_`echo $first_flag`_first_restarted_from_step_restart`echo $restart_step`
    

        if test "$1" = "no_fpdiff"; then
            echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
        else
            # Note: Should all agree with last file from non-restarted run!
            ../../../../../bin/fpdiff.py ../validata/result_`echo $first_flag`_first_for_restart.dat.gz result_`echo $first_flag`_first_restarted_from_step_restart`echo $restart_step`.dat   >> validation.log
        fi
    
    done
    
done






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
