#! /bin/bash

# Make the thing
make driven_cavity
make driven_cavity_without_singularity_subtracted_off
make adaptive_driven_cavity

# Setup result directory
dir=NEW_RUNS_DRIVEN_CAVITY
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
files_to_be_copied="compare_driven_cavity.pvsm compare_driven_cavity_different_re.pvsm  navier_stokes_elements_with_singularity.h adaptive_driven_cavity adaptive_driven_cavity.cc driven_cavity.cc driven_cavity driven_cavity_without_singularity_subtracted_off"
cp $files_to_be_copied $dir

# Move to new directory
cd $dir

# Reynolds number
Re_list="0.0 50.0 100.0 150.0 200.0"
for Re in `echo $Re_list`; do

    # Make new case director
    case_dir="Case_Re"$Re
    mkdir $case_dir
    cp $files_to_be_copied $case_dir
    cd $case_dir
    
    # Adaptive driven cavity
    mkdir RESLT 
    ./adaptive_driven_cavity --re $Re > RESLT/OUTPUT
    cd RESLT
    oomph-convert soln9.dat
    oomph-convert coarse_soln9.dat
    cd ..
    mv RESLT RESLT_adaptive_driven_cavity
    
    
    #With/without singularity removal
    sing_list="0 1"
    for sing in `echo $sing_list`; do

        executable=driven_cavity
        if [ $sing -eq 0 ]; then
            executable=driven_cavity_without_singularity_subtracted_off
        fi
        
        # Specify the number of elements
        element_multiplier_list="1 2 3 4 5"
        for element in `echo $element_multiplier_list`; do
            
            n_y=$((10*$element))
            
            echo "Running for element multiplier: " $element
            
            mkdir RESLT
            ./`echo $executable` --element_multiplier $element --re $Re > RESLT/OUTPUT
            
            
            cd RESLT
            sed 's/inf/50000/g' soln0.dat >soln0_filtered.dat; oomph-convert soln0_filtered.dat
            sed 's/inf/50000/g' coarse_soln0.dat > coarse_soln0_filtered.dat; oomph-convert coarse_soln0_filtered.dat
            cd ..
            
            cpu=`grep "Total time for Newton solver" RESLT/OUTPUT | awk '{print $8}'`
            error_with_pressure=`grep "error (with pressure)" RESLT/OUTPUT | awk '{print $4}'`
            error_without_pressure=`grep "error (without pressure)" RESLT/OUTPUT | awk '{print $4}'`
            
            
            sing_identifier="_with_singularity"
            if [ $sing -eq 0 ]; then
                sing_identifier="_without_singularity"
            fi
            new_result_dir="RESLT"$sing_identifier"_re"$Re"_el_multiplier"$element
            
            echo "Moving results to: " $new_result_dir
            mv RESLT $new_result_dir
            
            #echo $n_y " " $error_with_pressure" " $error_without_pressure " " $cpu >> convergence$sing_identifier.dat
            
        done
        
        vtu_file_list=`ls RESLT$sing_identifier*/coarse_soln0_filtered.vtu`
        /home/mheil/bin/make_combined_pvd.bash $vtu_file_list
        mv combined.pvd coarse_combined$sing_identifier".pvd"
        vtu_file_list=`ls RESLT$sing_identifier*/soln0_filtered.vtu`
        /home/mheil/bin/make_combined_pvd.bash $vtu_file_list
        mv combined.pvd combined$sing_identifier".pvd"
        
    done

    # Ready for next case
    cd ..

done

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
