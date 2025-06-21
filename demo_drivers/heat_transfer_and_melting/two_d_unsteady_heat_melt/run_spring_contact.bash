#! /bin/bash
case_list="--old_version --single_kink --completely_smooth --kuhn_tucker"

n_primary_list="1" # 40 80 120" 


for n_primary in `echo $n_primary_list`; do 
    for case in `echo $case_list`; do 
        dir=RESLT_`echo $case`_`echo $n_primary`
        mkdir $dir
        cd $dir
        echo "Doing "$dir
        ../spring_contact $case --n_primary $n_primary > OUTPUT 
        ../plot_spring_contact_landscape.bash > plot.out
        #cat newton_iter?.dat newton_iter??.dat | awk '{print $4}' > all_newton_iter.dat
        cd ..
    done
done