list=`ls animate*png`
for file in `echo $list`; do
    convert -resize 20% $file `basename $file png`gif 
done
gifmerge -l0 animate_displacement?.gif  animate_displacement??.gif > all.gif
