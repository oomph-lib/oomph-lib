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

# Find all directories with validate.sh in them:
validate_sh_list=`find demo_drivers self_test -name validate.sh`

# Prepare summary file
rm -f .failed_compilation_list
failed_file=`echo $current_dir"/.failed_compilation_list"`

for file in `echo $validate_sh_list`; do
    dir=`dirname $file`
    cd $dir
    make check-if-executables-exist | grep 'LIKELY' >> $failed_file
    cd $current_dir
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
    echo "Check these out!"
    return_flag=1
fi
echo " "

# Clean up
rm -f .failed_compilation_list

exit $return_flag
