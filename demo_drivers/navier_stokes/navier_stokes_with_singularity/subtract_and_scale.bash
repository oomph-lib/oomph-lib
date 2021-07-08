#! /bin/bash


step_number=9

Re_list="5.0 10.0 15.0 20.0"
for Re in `echo $Re_list`; do
    paste Case_Re$Re/RESLT/soln$step_number.dat  Case_Re0.0/RESLT/soln$step_number.dat | awk -v Re=$Re '{if ($1!="ZONE"){print $1" "$2" "($14-$3)/Re" "($15-$4)/Re" "($16-$5)/Re" "}else{print $1" "$2" "$3}}' > scaled_diff_from_stokes_$Re.dat
    paste Case_Re$Re/RESLT/line_plot$step_number.dat  Case_Re0.0/RESLT/line_plot$step_number.dat | awk -v Re=$Re '{if ($1!="ZONE"){print $1" "$2" "$3" "($10-$4)/Re" "($11-$5)/Re" "($12-$6)/Re" "}else{for (i=1;i<=NF/2;i++){printf "%s", $i" "}; print""}}' > scaled_line_plot_diff_from_stokes_$Re.dat
    echo "done scaled_diff_from_stokes_$Re.dat"
done
