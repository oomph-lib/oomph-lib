#! /bin/bash

file_list=`ls animate_displacement*.png`

for file in `echo $file_list`; do
    base=`basename $file png`
    new=`echo $base`gif
    echo $new
    convert -resize 20% $file $new
done


gifmerge -l0 animate_displacement?.gif animate_displacement??.gif > anim.gif
 