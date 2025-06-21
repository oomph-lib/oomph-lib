#! /bin/bash


root_dir=`pwd`
other_root_dir=/home/mheil/version_reference_against_louis_stuff/demo_drivers/mpi/distribution/two_d_unstructured_adaptive_poisson

dir_names_list=`ls -d Validation/RES*NP*`
for dir in `echo $dir_names_list`; do

    echo $dir
    echo " "
    echo "--------------------------------"
    echo " "
    cd $root_dir/$dir
    pwd
    ls -l soln2*.dat
    cat soln2*.dat > combined.dat    
    /home/mheil/bin/oomph-convert combined.dat
    echo " " 
    cd $other_root_dir/$dir
    ls -l soln2*.dat   
    cat soln2*.dat > combined.dat
    /home/mheil/bin/oomph-convert combined.dat
    pwd
    echo " "
    echo "--------------------------------"
    echo " "
    cd $root_dir
    rm combined_mine.vtu
    rm  combined_other.vtu
    ln -sf $root_dir/$dir"/combined.vtu" combined_mine.vtu
    ln -sf $other_root_dir/$dir"/combined.vtu" combined_other.vtu
    paraview --state=compare.pvsm

done
