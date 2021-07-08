#! /bin/bash

# Rebuild executable
make solid_contact

# Directory to store data
dir=./NEW_RUNS
if [ -e `echo $dir` ];  then
    echo " "
    echo "ERROR: Please delete directory $dir and try again"
    echo " "
    exit
fi
mkdir $dir

# Move to working directory
cp solid_contact       $dir
cp solid_contact.cc    $dir
cp plot_solid_contact.bash    $dir
cp solid_contact.pvsm    $dir
cd $dir


# Loop over the element areas
#----------------------------
el_area_list="0.02 0.002" ## hierher  0.0002 0.00002 0.000002" # adjust
for el_area in `echo $el_area_list`; do

    
    # Make and move to actual run directory
    mkdir RESLT
    run_dir="RESLT_el_area"`echo $el_area`
    echo "Dir name: " $run_dir

    
    # Assemble command line flag
    #command_line_flag=`echo " --proper_elasticity  --no_adapt --el_area  $el_area"` # hierher
    command_line_flag=`echo " --proper_elasticity  --el_area  $el_area"`
    
    # Run it
    ./solid_contact `echo $command_line_flag` > RESLT/OUTPUT
   
    # Postprocess
    cd RESLT
    ../plot_solid_contact.bash
    #paraview --state=../solid_contact.pvsm
    cd ..

    # Move
    mv RESLT $run_dir
    
    
done
cd ..