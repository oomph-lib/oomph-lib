#! /bin/bash

string='CPU for projection of solution onto new mesh'
string='CPU for (one way) aux_setup_multi_domain_interaction'
echo $string

file_names_list=`ls Validation/OUTPUT*`
for file in `echo $file_names_list`; do
    grep "$string" $file | awk '{print $NF}' > .junk_new.dat
    grep "$string" /home/mheil/version_reference_against_louis_stuff/demo_drivers/mpi/distribution/two_d_unstructured_adaptive_poisson/$file| awk '{print $NF}' > .junk_old.dat
    paste .junk_new.dat .junk_old.dat
done
