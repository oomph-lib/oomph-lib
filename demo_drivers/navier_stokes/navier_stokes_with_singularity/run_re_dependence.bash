#! /bin/bash


# Adaptive
adaptive=1

# Make the thing
executable=driven_cavity
post_fix=""
if [ $adaptive == 1 ]; then
    make adaptive_driven_cavity
    executable=adaptive_driven_cavity
    dir_postfix="_ADAPTIVE"
else
    make driven_cavity
    executable=driven_cavity
fi





# Setup result directory
dir=NEW_RUNS_RE_DEPENDENCE$dir_postfix
if [ -e $dir ]; then
    echo " "
    echo " "
    echo "Directory " $dir "already exists -- please move it out of the way."
    echo " "
    echo " "
    exit
fi
mkdir $dir

# Copy driver codes and stuff in for future reference
files_to_be_copied="compare_driven_cavity.pvsm compare_driven_cavity_different_re.pvsm  navier_stokes_elements_with_singularity.h adaptive_driven_cavity adaptive_driven_cavity.cc driven_cavity.cc driven_cavity driven_cavity_without_singularity_subtracted_off finite_re_perturbation.h"
cp $files_to_be_copied $dir

# Move to new directory
cd $dir

soln_file_list=""

# Reynolds number
Re_list="0.0 5.0 10.0 15.0 20.0"
for Re in `echo $Re_list`; do

    # Make new case director
    case_dir="Case_Re"$Re
    mkdir $case_dir
    cp $files_to_be_copied $case_dir
    cd $case_dir
    mkdir RESLT

    if [ $adaptive == 1 ]; then
        ./`echo $executable` --re $Re > RESLT/OUTPUT
    else

        # Specify the number of elements
        element_multiplier=21
        n_y=$((10*$element_multiplier))
        echo "Running for element multiplier: " $element_multiplier
        ./`echo $executable` --element_multiplier $element_multiplier --re $Re > RESLT/OUTPUT
    fi

    cd RESLT
    /home/mheil/bin/oomph-convert soln0.dat
    soln_file_list=`echo $soln_file_list`" "`pwd`"/soln0.vtu"
    cd ..

    #vtu_file_list=`ls RESLT$sing_identifier*/coarse_soln0_filtered.vtu`
    #/home/mheil/bin/make_combined_pvd.bash $vtu_file_list
    #mv combined.pvd coarse_combined$sing_identifier".pvd"
    #vtu_file_list=`ls RESLT$sing_identifier*/soln0_filtered.vtu`
    #/home/mheil/bin/make_combined_pvd.bash $vtu_file_list
    #mv combined.pvd combined$sing_identifier".pvd"
    
    # Ready for next case
    cd ..
    
done

echo "soln file list: " $soln_file_list
/home/mheil/bin/make_combined_pvd.bash  $soln_file_list


exit

#################################################



make_combined_pvd.bash Case_Re*/RESLT_with_singularity_re*_el_multiplier5/soln0_filtered.vtu; mv combined.pvd combined_with_singularity.pvd
make_combined_pvd.bash Case_Re*/RESLT_without_singularity_re*_el_multiplier5/coarse_soln0_filtered.vtu; mv combined.pvd combined_without_singularity.pvd
make_combined_pvd.bash Case_Re*/RESLT_with_singularity_re*_el_multiplier5/coarse_soln0_filtered.vtu; mv combined.pvd coarse_combined_with_singularity.pvd
make_combined_pvd.bash Case_Re*/RESLT_ad*/soln9.vtu; mv combined.pvd combined_adaptive.pvd
make_combined_pvd.bash Case_Re*/RESLT_ad*/coarse_soln9.vtu; mv combined.pvd coarse_combined_adaptive.pvd

echo  " "
echo  "Visualise with "
echo  " "
echo  "-- in Case* directory: "
echo  " "
echo  "      paraview --state=compare_driven_cavity.pvsm"
echo  " "
echo  "and "
echo  " "
echo  "-- in top level directory: "
echo  " "
echo  "      paraview --state=compare_driven_cavity_different_re.pvsm"
echo  " "
