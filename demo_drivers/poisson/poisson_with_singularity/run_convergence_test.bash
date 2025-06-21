
#! /bin/bash

# Make the thing
make two_d_poisson two_d_poisson_without_singularity_subtracted_off

# Setup result directory
dir=NEW_RUNS
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
cp two_d_poisson.cc two_d_poisson two_d_poisson_without_singularity_subtracted_off poisson_elements_with_singularity.h $dir

# Move to new directory
cd $dir

# Specify the number of elements
element_multiplier_list="1 2 3 4 5" # 6 7 8 9 10"

for element in `echo $element_multiplier_list`; do
        
    # With singularity
    #-----------------
    echo "Running with singul for element multiplier: " $element
    mkdir RESLT
    ./two_d_poisson --element_multiplier $element > RESLT/OUTPUT
    # Postprocess results
    n_y=$((4*$element))
    cpu=`grep "Total time for Newton solver" RESLT/OUTPUT | awk '{print $8}'`
    error=`grep "Norm of error" RESLT/OUTPUT | awk '{print $5}'`
    new_result_dir="RESLT_el_multiplier"$element
    echo "Writing results to: " $new_result_dir
    mv RESLT $new_result_dir
    
    echo $n_y " " $error " " $cpu >> convergence_with_sing.dat



    # Without singularity
    #--------------------
    echo "Running without singul for element multiplier: " $element
    mkdir RESLT
    ./two_d_poisson_without_singularity_subtracted_off --element_multiplier $element > RESLT/OUTPUT
    # Postprocess results
    n_y=$((4*$element))
    cpu=`grep "Total time for Newton solver" RESLT/OUTPUT | awk '{print $8}'`
    error=`grep "Norm of error" RESLT/OUTPUT | awk '{print $5}'`
    new_result_dir="RESLT_without_singularity_subtracted_off_el_multiplier"$element
    echo "Writing results to: " $new_result_dir
    mv RESLT $new_result_dir
    
    echo $n_y " " $error " " $cpu >> convergence_without_sing.dat

    
done

exit
