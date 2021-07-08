#! /bin/bash

# Rebuild executable
make stefan_boltzmann_melt

# Directory to store data
dir=./NEW_RUNS_nintpt
if [ -e `echo $dir` ];  then
    echo " "
    echo "ERROR: Please delete directory $dir and try again"
    echo " "
    exit
fi
mkdir $dir

# Move to working directory
cp stefan_boltzmann_melt       $dir
cp stefan_boltzmann_melt.cc    $dir
cd $dir




# Loop over timesteps
#--------------------
dt_list="0.005" # 0.001"
for dt in `echo $dt_list`; do
    
    
# Loop over timesteps
#--------------------
nintpt_list="10 20"
for nintpt_aux in `echo $nintpt_list`; do
    
    
    nintpt_flag=""
    nintpt_dirflag="_gauss"
    if [ $nintpt_aux -gt 0 ]; then
        nintpt_flag="--nintpt  "$nintpt_aux
        nintpt_dirflag="_nintpt"$nintpt_aux
    fi

    # Make and move to actual run directory
    run_dir="RESLT_dt"`echo $dt``echo $nintpt_dirflag`
    echo "Dir name: " $run_dir
    mkdir $run_dir
    
    # Assemble command line flag
    command_line_flag=`echo " --dt $dt --dir $run_dir $nintpt_flag"`
    
    # Run it
    ./stefan_boltzmann_melt `echo $command_line_flag` > `echo $run_dir`/OUTPUT &
    
done
done

cd ..
