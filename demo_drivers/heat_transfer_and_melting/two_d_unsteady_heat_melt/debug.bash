#! /bin/bash

compare_jacs.bash cr_jac.dat fd_jac.dat
../../../bin/fpdiff.py in_both_file*.dat 1.0 1.0e-13 > diff.dat

rm -f offending_rows.dat
rm -f offending_row_dofs.dat
awk '{if (substr($0,1,1) ~ /^[0-9][0-9]*$/ ){print $1}}' diff.dat | unique  offending_rows.dat
list=`cat offending_rows.dat`
for row in `echo $list`; do 
    echo "Finding dof type associated with row "$row
    grep eqn_number typescript | awk -v eqn=$row '{if ($10==eqn){print $0}}' >> offending_row_dofs.dat
done



rm -f offending_cols.dat
rm -f offending_col_dofs.dat
awk '{if (substr($0,1,1) ~ /^[0-9][0-9]*$/ ){print $2}}' diff.dat | unique  offending_cols.dat
list=`cat offending_cols.dat`
for col in `echo $list`; do 
    echo "Finding dof type associated with col "$col
    grep eqn_number typescript | awk -v eqn=$col '{if ($10==eqn){print $0}}'>> offending_col_dofs.dat
done
