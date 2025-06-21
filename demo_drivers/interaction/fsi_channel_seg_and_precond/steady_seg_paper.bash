#!/bin/bash


#=====================================================
# Shell script to do steady runs with simple
# fsi collapsible channel driver. Potentially
# loop over:
# - spatial resolutions,
# - segregated/monolithic solver
# - Aitken extrapolation (on/off)
# - Irons and Tuck extrapolation (on/off)
# - displacement control (on/off)
# - Q,
# - Re
#
#======================================================

host_names="\"compute-1-0 compute-1-1 compute-1-2 compute-1-3 compute-1-4 compute-1-5 compute-1-6 compute-1-7 compute-1-8 compute-1-9 compute-1-10 compute-1-11 compute-1-12 compute-1-13 compute-1-14 compute-1-16\"" 

# ===================================================
# Specify directory in which job is to be run
# ===================================================
RUN_DIR=STEADY_SEGREGATED_COLLAPSIBLE_CHANNEL_RUNS



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
resolution_list="2 3 4"
for resolution in `echo $resolution_list`; do

  # Segregated solver?
  #-------------------
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
     omega_list="1.0e-0 0.9e-0 0.8e-0 0.7e-0"
     # Under-relaxation only for segregated solver
     if [ $use_segregated_solver -eq 0 ]; then omega_list="1.0e-0" ; fi
     for omega in `echo $omega_list`; do
        
     # Irons and Tuck extrapolation
     #-----------------------------
     irons_and_tuck_list="0 1"
     # Irons and Tuck only for segregated solver
     if [ $use_segregated_solver -eq 0 ]; then irons_and_tuck_list="0" ; fi
     for use_irons_and_tuck in `echo $irons_and_tuck_list`; do


      # Displacement control
      #---------------------
      displ_ctrl_list="1"
      for use_displ_ctrl in `echo $displ_ctrl_list`; do

      # Max. prescribed y_ctrl for displacement control
      #------------------------------------------------
      y_prescr_max=0.65

      # Only steady run
      #-----------------
      steady_flag=1

       # Loop over FSI parameters
       #--------------------------
       Q_list="1.0e-2 1.0e-3 1.0e-4"
       for Q in `echo $Q_list`; do

         # Loop over Reynolds numbers
         #---------------------------
         Re_list="0 250 500"
         for Re in `echo $Re_list` ; do

          # Convergence criterion (1 = max. change in displacement)
          #--------------------------------------------------------
          convergence_crit=1

           # Loop over convergence tol
           #--------------------------------
           convergence_tol_list="1.0e-8"
           if [ $use_segregated_solver -eq 0 ]; then convergence_tol_list="0" ; fi
           for convergence_tol in `echo $convergence_tol_list` ; do


          # Number of steps in parameter study
          nsteps=6 #20
          
          # Dummy timestep
          dt=0.1

          # Create flags for executable
          flags=`echo $resolution $use_segregated_solver $use_aitken $omega $use_irons_and_tuck $use_displ_ctrl $y_prescr_max $steady_flag $convergence_crit $convergence_tol $nsteps $Re $Q $dt`

#          flags=`echo $resolution $use_segregated_solver $use_aitken $omega $use_irons_and_tuck $use_displ_ctrl $y_prescr_max $steady_flag $convergence_criterion $nsteps $Re $Q $dt`


          # Doc which case we're running
          echo " "
          echo "Running executable with flags: " $flags


          # Remember where we are: The directory where all these jobs
          # are spawned
          SPAWN_DIR=`pwd`

          # Make a temporary directory
          BATCH_DIR=`mktemp -d $SPAWN_DIR/tmp_dir.XXXXXX`
          #echo "Temporarily, results for this run are in: " $BATCH_DIR

          # Copy executable in to exec directory
          cp fsi_collapsible_channel_segregated_driver $BATCH_DIR


          # Create shell script for this job
          #---------------------------------
          # Note: Script needs to be
          # in local dir not in /tmp because /tmp differs from
          # compute node to compute node!
          RUN_SCRIPT=`mktemp run_script.XXXXXX`
          #echo "Run script: " $RUN_SCRIPT

          echo "#!/bin/bash" > $RUN_SCRIPT
          echo "cd $BATCH_DIR" >> $RUN_SCRIPT
          echo "rm -rf RESLT"  >> $RUN_SCRIPT
          echo "mkdir RESLT"  >> $RUN_SCRIPT
          echo "echo $flags > command_line_flags" >> $RUN_SCRIPT
          echo "cp command_line_flags  RESLT"  >> $RUN_SCRIPT
          echo "./fsi_collapsible_channel_segregated_driver \`cat command_line_flags\` > OUTPUT" >> $RUN_SCRIPT
          echo "cd $SPAWN_DIR"  >> $RUN_SCRIPT
          solver_flag="newton"
          if [ $use_segregated_solver -eq 1 ]; then solver_flag="picard" ; fi 
          echo "post_fix="`echo "res"$resolution"_aitken"$use_aitken"_omega"$omega"_irons_and_tuck"$use_irons_and_tuck"_conv_crit"$convergence_crit"_conv_tol"$convergence_tol"_Re"$Re"_Q"$Q"_solver"$solver_flag`  >> $RUN_SCRIPT

#          echo "post_fix="`echo $resolution"_"$use_displ_ctrl"_"$Re"_"$Q"_"$solver_flag"_"$use_aitken"_"$use_irons_and_tuck`  >> $RUN_SCRIPT

          echo "cat $BATCH_DIR/RESLT/soln5.dat >> $SPAWN_DIR/steady.dat" >> $RUN_SCRIPT
          echo "mv $BATCH_DIR/RESLT $SPAWN_DIR/RESLT_\$post_fix"  >> $RUN_SCRIPT
          echo "mv $BATCH_DIR/OUTPUT $SPAWN_DIR/OUTPUT_\$post_fix"  >> $RUN_SCRIPT
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

          #$BATCH_DIR/$RUN_SCRIPT

          done #conv. tolerance

         done # Re

       done #Q

      done #displ control

      done # irons and tuck

     done # omega

    done # aitken

 done # solver

done # resolution

