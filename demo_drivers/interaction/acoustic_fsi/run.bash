#! /bin/bash

make acoustic_fsi

case_postfix="_refine"
case_postfix="_outer_radius"

dir=runs`echo $case_postfix`
# Check existence of working dir
if [ -e `echo $dir`   ] 
then
    echo " "
    echo "ERROR: Please delete directory $dir and try again"
    echo " "
    exit
fi
mkdir $dir
cd $dir


if [ $case_postfix == "_refine" ]; then

    el_mult_list="1 4 8 16"
    for el_mult in `echo $el_mult_list`; do
        
        mkdir RESLT
        ../acoustic_fsi --el_multiplier $el_mult --max_adapt 0 --nstep 10 > RESLT/OUTPUT 
        mv RESLT RESLT_`echo $el_mult`
        
    done

elif [ $case_postfix == "_outer_radius" ]; then

    r_list="4.0 3.0 2.0 1.5"
    for r in `echo $r_list`; do
        
        mkdir RESLT
        ../acoustic_fsi --el_multiplier 3 --outer_radius $r --max_adapt 0 --nstep 10 > RESLT/OUTPUT 
        mv RESLT RESLT_r`echo $r`
        
    done

fi