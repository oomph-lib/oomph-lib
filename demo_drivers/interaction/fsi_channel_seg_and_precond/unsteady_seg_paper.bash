#!/bin/bash


#=====================================================
# Shell script to do unsteady runs with simple
# fsi collapsible channel driver. Potentially
# loop over:
# - spatial resolutions, 
# - segregated/monolithic solver
# - Aitken extrapolation (on/off)
# - under-relaxation parameter
# - Irons and Tuck extrapolation (on/off)
# - Q,
# - Re
#
# 
# Procedure: 
#
#    1) Do steady run up to a certain 
#       control displacement, using the monolithic 
#       solver
#
#    2) Do a restart from the steady solution
#       with an incremented external pressure.
#       Run done with the monolithic solver for a small
#       number of timesteps
#
#    3) Restart time-dep. run with segregated solver
#       (if requested) to avoid impulsive start 
#       problems.
#
#======================================================

host_names="\"compute-1-0 compute-1-1 compute-1-2 compute-1-3 compute-1-4 compute-1-5 compute-1-6 compute-1-7 compute-1-8 compute-1-9 compute-1-10 compute-1-11 compute-1-12 compute-1-13 compute-1-14  compute-1-15 compute-1-16\"" 


# ===================================================
# Specify directory in which job is to be run
# ===================================================
RUN_DIR=UNSTEADY_SEGREGATED_COLLAPSIBLE_CHANNEL_RUNS


# Make directory; backup if necessary
if [ -d $RUN_DIR ] ; then
   backed_up_dir=`mktemp -d $RUN_DIR.XXXXXX`
   mv $RUN_DIR $backed_up_dir
   echo " "
   echo "==================================================================="
   echo " "
   echo "Run directory already existed: Backed up into: " $backed_up_dir
   echo " "
   echo "==================================================================="
   echo " "
fi
mkdir $RUN_DIR
cp fsi_collapsible_channel_segregated_driver $RUN_DIR
cd $RUN_DIR




#========================================================
#PARAMETER STUDIES
#========================================================

# Different spatial resolutions
#------------------------------
resolution_list="1 2 3 4"
resolution_list="3"
for resolution in `echo $resolution_list`; do

  # Segregated solver?
  #-------------------
  solver_list="0 1"
  solver_list="0 1"
  for use_segregated_solver in `echo $solver_list`; do

    # Aitken extrapolation
    #---------------------
    aitken_list="0"
    # Aitken only for segregated solver
    if [ $use_segregated_solver -eq 0 ]; then aitken_list="0" ; fi
    for use_aitken in `echo $aitken_list`; do
        
     # Under-relaxation
     #-----------------
     omega_list="1.0e-4 0.25e-3 0.5e-3 0.75e-3 1.0e-3 0.25e-2 0.5e-2 0.75e-2"
     omega_list="1.0e-3 0.5e-3 1.0e-4"
     # Under-relaxation only for segregated solver
     if [ $use_segregated_solver -eq 0 ]; then omega_list="1.0e-0" ; fi
     for omega in `echo $omega_list`; do
        
      # Irons and Tuck extrapolation
      #-----------------------------
      irons_and_tuck_list="0 1"
      # Irons and Tuck only for segregated solver
      if [ $use_segregated_solver -eq 0 ]; then irons_and_tuck_list="0" ; fi
      if [ $use_aitken -eq 1 ]; then irons_and_tuck_list="0" ; fi
      for use_irons_and_tuck in `echo $irons_and_tuck_list`; do

       # Loop over FSI parameters
       #--------------------------
        Q_list="1.0e-2 1.0e-3 1.0e-4"
        Q_list="1.0e-2"
        for Q in `echo $Q_list`; do
            
         # Loop over Reynolds numbers
         #---------------------------
         Re_list="0 250 500"
         Re_list="500"
         for Re in `echo $Re_list` ; do
                
           # Always use displacement control for steady run
           #-----------------------------------------------
           use_displ_ctrl=1

           # Max. prescribed y_ctrl for displacement control
           #------------------------------------------------
           y_prescr_max=0.65

           # Timestep
           #---------
           dt=0.1

           # Convergence crit
           #----------------------
           convergence_crit=1


           # Loop over convergence tol
           #--------------------------------
           convergence_tol_list="1.0e-2 1.0e-4 1.0e-6 1.0e-8"
           convergence_tol_list="1.0e-8"
           if [ $use_segregated_solver -eq 0 ]; then convergence_tol_list="0" ; fi
           for convergence_tol in `echo $convergence_tol_list` ; do

           echo " " 
           echo " " 
           echo "##########################################################" 
           echo "##########################################################" 
           echo "##########################################################" 
           echo " " 
           echo " " 
           
           
           
           # Create flags for steady run (always do it with Newton)
           #-------------------------------------------------------
           steady_flag=1
           nsteps=5
           dummy_use_segregated_solver=0
           steady_use_aitken=0

           steady_flags=`echo $resolution $dummy_use_segregated_solver $steady_use_aitken $omega $use_irons_and_tuck $use_displ_ctrl $y_prescr_max $steady_flag $convergence_crit $convergence_tol $nsteps $Re $Q $dt`
           
           
          # Doc which case we're running
           echo " "
           echo "Doing steady run with flags  : " $steady_flags
           
           
           # Create flags for first unsteady run (with Newton)
           #--------------------------------------------------
           restart_directory="RESLT_steady"
           restart_file_number=$nsteps
           steady_flag=0
           first_nsteps=5
           first_pressure_jump=0.8333 # as in paper # hierher add this to loop later
           first_use_segregated_solver=0
           first_aitken_flag=0

           #if [ $Q = "1.0e-3" ]; then first_pressure_jump=0.4333 ; fi
           #if [ $Q = "1.0e-4" ]; then first_pressure_jump=0.2333 ; fi

           echo "First unsteady run with pressure_jump = "  $first_pressure_jump

           first_unsteady_flags=`echo $resolution $first_use_segregated_solver $first_aitken_flag $omega $use_irons_and_tuck $use_displ_ctrl $y_prescr_max $steady_flag $convergence_crit $convergence_tol $first_nsteps $Re $Q $dt $restart_directory"/restart"$restart_file_number".dat" $first_pressure_jump`
           
           
          # Doc which case we're running
           echo " "
           echo "Doing first unsteady run with flags: " $first_unsteady_flags
           


           # Create flags for second unsteady run (with whatever solver)
           #------------------------------------------------------------
           restart_directory="RESLT_unsteady"
           overlap=3
           restart_file_number=`echo "$first_nsteps - $overlap" | bc`
           steady_flag=0
           nsteps=200
           pressure_jump=0.0 # we're not doing any further impulsive increase

           echo "Running with pressure_jump = "  $pressure_jump

           unsteady_flags=`echo $resolution $use_segregated_solver $use_aitken $omega $use_irons_and_tuck $use_displ_ctrl $y_prescr_max $steady_flag $convergence_crit $convergence_tol $nsteps $Re $Q $dt $restart_directory"/restart"$restart_file_number".dat" $pressure_jump`
           
           
          # Doc which case we're running
           echo " "
           echo "Doing continued unsteady run with flags: " $unsteady_flags
           

         
          # Remember where we are: The directory where all these jobs
          # are spawned
           SPAWN_DIR=`pwd`
           
          # Make a temporary directory
          BATCH_DIR=`mktemp -d $SPAWN_DIR/tmp_dir.XXXXXX`
          echo "Temporarily, results for this run are in: " $BATCH_DIR

          # Copy executable in to exec directory
          cp fsi_collapsible_channel_segregated_driver $BATCH_DIR


          # Create shell script for this job
          #---------------------------------
          # Note: Script needs to be
          # in local dir not in /tmp because /tmp differs from
          # compute node to compute node!
          RUN_SCRIPT=`mktemp run_script.XXXXXX`
          echo "Run script: " $RUN_SCRIPT

          echo "#!/bin/bash" > $RUN_SCRIPT
          echo "cd $BATCH_DIR" >> $RUN_SCRIPT
          echo "rm -rf RESLT"  >> $RUN_SCRIPT

          echo "mkdir RESLT"  >> $RUN_SCRIPT
          echo "echo $steady_flags > command_line_flags" >> $RUN_SCRIPT
          echo "cp command_line_flags  RESLT"  >> $RUN_SCRIPT
          echo "./fsi_collapsible_channel_segregated_driver \`cat command_line_flags\` > OUTPUT_steady" >> $RUN_SCRIPT
          echo "mv RESLT RESLT_steady"  >> $RUN_SCRIPT

          echo "mkdir RESLT"  >> $RUN_SCRIPT
          echo "echo $first_unsteady_flags > command_line_flags" >> $RUN_SCRIPT
          echo "cp command_line_flags  RESLT"  >> $RUN_SCRIPT
          echo "./fsi_collapsible_channel_segregated_driver \`cat command_line_flags\` > OUTPUT_unsteady" >> $RUN_SCRIPT
          echo "mv RESLT RESLT_unsteady"  >> $RUN_SCRIPT

          echo "mkdir RESLT"  >> $RUN_SCRIPT
          echo "echo $unsteady_flags > command_line_flags" >> $RUN_SCRIPT
          echo "cp command_line_flags  RESLT"  >> $RUN_SCRIPT
          echo "./fsi_collapsible_channel_segregated_driver \`cat command_line_flags\` > OUTPUT_unsteady2" >> $RUN_SCRIPT
          echo "mv RESLT RESLT_unsteady2"  >> $RUN_SCRIPT

          echo "cd $SPAWN_DIR"  >> $RUN_SCRIPT
          solver_flag="newton"
          if [ $use_segregated_solver -eq 1 ]; then solver_flag="picard" ; fi 
          echo "post_fix="`echo "res"$resolution"_aitken"$use_aitken"_omega"$omega"_irons_and_tuck"$use_irons_and_tuck"_conv_crit"$convergence_crit"_conv_tol"$convergence_tol"_Re"$Re"_Q"$Q"_solver"$solver_flag`  >> $RUN_SCRIPT
          echo "cat $BATCH_DIR/RESLT_unsteady2/soln2.dat >>  $SPAWN_DIR/unsteady.dat" >> $RUN_SCRIPT
        echo "mv $BATCH_DIR/RESLT_steady $SPAWN_DIR/RESLT_steady_\$post_fix"  >> $RUN_SCRIPT
          echo "mv $BATCH_DIR/RESLT_unsteady $SPAWN_DIR/RESLT_unsteady_\$post_fix"  >> $RUN_SCRIPT
          echo "mv $BATCH_DIR/RESLT_unsteady2 $SPAWN_DIR/RESLT_unsteady2_\$post_fix"  >> $RUN_SCRIPT
          echo "mv $BATCH_DIR/OUTPUT_steady $SPAWN_DIR/OUTPUT_steady_\$post_fix"  >> $RUN_SCRIPT
          echo "mv $BATCH_DIR/OUTPUT_unsteady $SPAWN_DIR/OUTPUT_unsteady_\$post_fix"  >> $RUN_SCRIPT
          echo "mv $BATCH_DIR/OUTPUT_unsteady2 $SPAWN_DIR/OUTPUT_unsteady2_\$post_fix"  >> $RUN_SCRIPT

          echo "rm -rf $BATCH_DIR" >> $RUN_SCRIPT
          chmod a+x  $RUN_SCRIPT

          #echo " " 
          #echo "-----------FYI: Here's the run script--------------" 
          #echo " " 
          #cat  $RUN_SCRIPT
          #echo " " 
          #echo "---------------------------------------------------" 
          #echo " " 
          mv $RUN_SCRIPT $BATCH_DIR

          # Run it via bsub
          echo "Running job....."

          command="bsub -m `echo $host_names` $BATCH_DIR/$RUN_SCRIPT -o $BATCH_DIR/../$RUN_SCRIPT.log" 

          echo "COMMAND: " $command 
          echo $command > tmp_command
          source tmp_command & 
          rm tmp_command

          # $BATCH_DIR/$RUN_SCRIPT 


          done # convergence tol

         done # Re

       done #Q

      done # irons and tuck

     done # omega

    done # aitken

 done # solver

done # resolution













