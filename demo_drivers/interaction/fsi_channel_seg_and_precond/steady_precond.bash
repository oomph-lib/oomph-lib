#!/bin/bash

#=====================================================
# Shell script to do steady runs with simple
# fsi collapsible channel driver. Potentially
# loop over:
#
# - spatial resolutions, 
# - displacement control (on/off)
# - Solver: direct or iterative with various preconditioners
# - preconditioner options (where appropriate)
# - Q,
# - Re
#
# Solve is always done with monolithic Newton solver
# but with different solvers/preconditioners.
#======================================================

#host_names="\"compute-1-0 compute-1-1 compute-1-2 compute-1-3 compute-1-4 compute-1-5 compute-1-6 compute-1-7 compute-1-8 compute-1-9 compute-1-10 compute-1-11 compute-1-12 compute-1-13 compute-1-14 compute-1-16\"" 


# ===================================================
# Specify directory in which job is to be run
# ===================================================
RUN_DIR=STEADY_PRECONDITIONED_COLLAPSIBLE_CHANNEL_RUNS


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
cp fsi_collapsible_channel_precond_driver $RUN_DIR
cd $RUN_DIR



#========================================================
#PARAMETER STUDIES
#========================================================

# Different spatial resolutions
#------------------------------
resolution_list="1 2 3 4 5 6 7 8 9"
resolution_list="2"
for resolution in `echo $resolution_list`; do
    
  echo "resolution: " $resolution
    
      # Displacement control
      #---------------------
      displ_ctrl_list="0 1"
      displ_ctrl_list="1"
      for use_displ_ctrl in `echo $displ_ctrl_list`; do


      # Solver/preconditioner
      #---------------------
      solver_flag_list="0 1 2 3"
      for solver_flag in `echo $solver_flag_list`; do


      # Preconditioner settings (where appropriate -- not for direct solver
      #--------------------------------------------------------------------
      # or exact preconditioner)
      #-------------------------
      solver_sub_flag_list="0 1 2"
      if [ $solver_flag -eq 0 ]; then solver_sub_flag_list="0" ; fi
      if [ $solver_flag -eq 1 ]; then solver_sub_flag_list="0" ; fi
      for solver_sub_flag in `echo $solver_sub_flag_list`; do

          echo " "
          echo " "
          echo " "
          echo "Solver_flag " $solver_flag
          echo "Solver_sub_flag_list " $solver_sub_flag_list
          echo " "
          echo " "
          echo " "

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

          # Only steady run
          #-----------------
          steady_flag=1


          # Max. prescribed y_ctrl for displacement control
          #------------------------------------------------
          y_prescr_max=0.65


          # Number of steps in parameter study
          nsteps=6

          echo " " 
          echo " " 
          echo "##########################################################" 
          echo "##########################################################" 
          echo "##########################################################" 
          echo " " 
          echo " " 

          # Create flags for executable
          flags=`echo $resolution $use_displ_ctrl $y_prescr_max $steady_flag $nsteps $solver_flag  $solver_sub_flag $Re $Q `

          # Doc which case we're running
          echo " "
          echo "Running executable with flags: " $flags

          # Remember where we are: The directory where all these jobs
          # are spawned
          SPAWN_DIR=`pwd`

          # Make a temporary directory
          BATCH_DIR=`mktemp -d $SPAWN_DIR/tmp_dir.XXXXXX`
          echo "Temporarily, results for this run are in: " $BATCH_DIR

          # Copy executable in to exec directory
          cp fsi_collapsible_channel_precond_driver $BATCH_DIR


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
          echo "echo $flags > command_line_flags" >> $RUN_SCRIPT
          echo "cp command_line_flags  RESLT"  >> $RUN_SCRIPT
          echo "./fsi_collapsible_channel_precond_driver \`cat command_line_flags\` > OUTPUT" >> $RUN_SCRIPT
          echo "cd $SPAWN_DIR"  >> $RUN_SCRIPT
          echo "post_fix="`echo "res"$resolution"_displ_ctrl"$use_displ_ctrl"_solver"$solver_flag"_subsolver"$solver_sub_flag"_Re"$Re"_Q"$Q`  >> $RUN_SCRIPT
          echo "mv $BATCH_DIR/RESLT $SPAWN_DIR/RESLT_\$post_fix"  >> $RUN_SCRIPT
          echo "mv $BATCH_DIR/OUTPUT $SPAWN_DIR/OUTPUT_\$post_fix"  >> $RUN_SCRIPT
          echo "rm -rf $BATCH_DIR/*" >> $RUN_SCRIPT
          chmod a+x  $RUN_SCRIPT

          echo " " 
          echo "-----------FYI: Here's the run script--------------" 
          echo " " 
          cat  $RUN_SCRIPT
          echo " " 
          echo "---------------------------------------------------" 
          echo " " 
          mv $RUN_SCRIPT $BATCH_DIR

          # Run it via bsub
          echo "Running job....."


          #command="bsub -m `echo $host_names` $BATCH_DIR/$RUN_SCRIPT -o $BATCH_DIR/../$RUN_SCRIPT.log" 
          command="$BATCH_DIR/$RUN_SCRIPT" 

          echo "COMMAND: " $command 
          echo $command > tmp_command
          source tmp_command 
          rm tmp_command

         done # Re

       done #Q

       done # solver sub flag

       done # solver flag

      done #displ control

done # resolution

