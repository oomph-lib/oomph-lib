
#------------------------------------------------------------------------
# Script to run all versions of the "sliding, clamped circular FvK plate"
#------------------------------------------------------------------------


# Setup directories
main_dir=NEW_RUNS
if [ -e $main_dir ]; then
    echo "Directory " $main_dir "already exists -- please move it out of the way."
    exit
fi
mkdir $main_dir

driver_list="circular_disk displacement_based_circular_disk axisym_fvk axisym_displ_based_fvk"

for driver in `echo $driver_list`; do

    make $driver
    cp $driver $main_dir
    cd $main_dir
    mkdir RESLT
    ./$driver > RESLT/OUTPUT
    mv RESLT RESLT_$driver
    cd ..

done
