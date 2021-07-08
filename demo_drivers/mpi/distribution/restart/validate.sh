#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


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
mkdir RESLT

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
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
