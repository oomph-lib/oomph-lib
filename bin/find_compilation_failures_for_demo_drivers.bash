#! /bin/bash

#===================================================================
# Script to search through all directories in demo_drivers and self_test
# directories that contain validate.sh scripts. In all these
# directories check if all executables required for the self tests
# have been built succesfully. Used to spot compilation failures that
# could go unobserved in a "make -k check" based self test.
#===================================================================


# Where are we?
current_dir=`pwd`

# Check if we're running mpi self-tests
mpi_run_command_length_aux=`make do-we-have-mpi`

#echo "mpi_run_command_length_aux: "$mpi_run_command_length_aux


mpi_run_command_length=`echo $mpi_run_command_length_aux | awk '{first=match($0, "START_BLA")+9; second=match($0, "END_BLA"); lngth=second-first; print substr($0, first, lngth)}'`

echo "mpi_run_command_length: "$mpi_run_command_length

have_mpi=0
if [ "$mpi_run_command_length" -ne 0 ]; then
    have_mpi=1
    echo "Have mpi compiler..."
else
    echo "Don't have mpi compiler..."
fi

# Find all directories with validate.sh in them:
validate_sh_list=`find demo_drivers self_test -name validate.sh`

# Prepare summary file
rm -f .failed_compilation_list
failed_file=`echo $current_dir"/.failed_compilation_list"`

for file in `echo $validate_sh_list`; do
    dir=`dirname $file`
    do_it=1
    if [ $have_mpi -eq 0 ]; then
        if [[ "$dir" =~ "/mpi/" ]]; then
            #echo "Dir "$dir" contains mpi but we don't have mpi compiler..."
            do_it=0
        fi
    fi
    if [ $do_it -eq 1 ]; then
        cd $dir
        #echo "executing make check-if-...-exist in:"
        #pwd
        make check-if-executables-exist | grep 'LIKELY' >> $failed_file
        cd $current_dir
    fi
done


echo " "
echo "==========================================================="
echo "Compilation of demo driver/self test codes failed in the"
echo "following directories: "
echo "==========================================================="
echo " "
failed_list=`awk '{print $NF}' $failed_file `
for failed in `echo $failed_list`; do dirname $failed; done | sort --unique
n_failed=`echo $failed_list | wc -w`
return_flag=0
if [ $n_failed -eq 0 ]; then
    echo "No failures!"
    return_flag=0
else
    echo " "
    echo "Error -- check these out!"
    return_flag=1
fi
echo " "

# Clean up
rm -f .failed_compilation_list

exit $return_flag
