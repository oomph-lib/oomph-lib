
#! /bin/bash

# Make the thing
make circular_couette
make circular_couette_without_pseudo_singularity_subtracted_off

# Setup result directory
dir=NEW_RUNS_CIRCULAR_COUETTE
if [ -e $dir ]; then
    echo " "
    echo " "
    echo "Directory " $dir "already exists -- please move it out of the way."
    echo " "
    echo " "
    exit
fi
mkdir $dir

# Copy driver code in for future reference
cp convergence_and_cpu_couette.lay driven_cavity.cc circular_couette circular_couette_without_pseudo_singularity_subtracted_off $dir

# Move to new directory
cd $dir

#With/without singularity removal
sing_list="0 1"
for sing in `echo $sing_list`; do

    executable=circular_couette
    if [ $sing -eq 0 ]; then
        executable=circular_couette_without_pseudo_singularity_subtracted_off
    fi
        
    # Specify the number of elements
    element_multiplier_list="1 2 3 4 5"
    for element in `echo $element_multiplier_list`; do
        
        n_y=$((10*$element))
        
        echo "Running for element multiplier: " $element
        
        mkdir RESLT
# hierher blend and Re=10.0 --blend 
        ./`echo $executable` --element_multiplier $element --re 10.0 > RESLT/OUTPUT
        cpu=`grep "Total time for Newton solver" RESLT/OUTPUT | awk '{print $NF}'`
        error_with_pressure=`grep "error (with pressure)" RESLT/OUTPUT | awk '{print $4}'`
        error_without_pressure=`grep "error (without pressure)" RESLT/OUTPUT | awk '{print $4}'`
        
        
        sing_identifier="_with_singularity"
        if [ $sing -eq 0 ]; then
            sing_identifier="_without_singularity"
        fi
        new_result_dir="RESLT"$sing_identifier"_re10_el_multiplier"$element

        echo "Moving results to: " $new_result_dir
        mv RESLT $new_result_dir
        
        echo $n_y " " $error_with_pressure" " $error_without_pressure " " $cpu >> convergence$sing_identifier.dat
        
    done
    
done
echo  " "
echo  "Visualise with "
echo  " "
echo  "    tecplot convergence_and_cpu_couette.lay"
echo  " "
