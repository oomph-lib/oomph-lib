#!/bin/bash


# Loop over different resolutions
#--------------------------------
resolution_list="50 100 150 200 250 300"
for resolution in `echo $resolution_list`; do

    # Loop over refinement strategies
    strategy_list="1 2 3 4"
    for strategy in `echo $strategy_list`; do

        #Don't dump the jacobians
        dump_jacobians=0

        #Halt code to doc memory usage
        pause_to_doc_memory=0

        # Create flags for executable
        flags=`echo $resolution $resolution $strategy $dump_jacobians $pause_to_doc_memory`

        # Create log file name
        log_file=`echo "log_file_"$resolution"x"$resolution"_method"$strategy`
        
        # Doc which case we're running
        echo " "
        echo "Running executable with flags: " $flags
        echo "Log file: " $log_file
          

        #bsub -o `echo $log_file` ./sparse_assemble_test `echo $flags` 
        ./sparse_assemble_test `echo $flags` 

    done # strategy

done #resolution
