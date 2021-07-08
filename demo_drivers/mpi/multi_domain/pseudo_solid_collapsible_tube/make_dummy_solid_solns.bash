#! /bin/bash

file_list=`ls RESLT*/solid_soln*.dat`

for file in `echo $file_list`; do
    
    non_empty=`wc $file | awk '{print $1}'`

    if [ $non_empty -eq 0 ]; then
        echo "ZONE"  > $file
        echo "0 0 0 0 0 0 0" >> $file
    fi

done

