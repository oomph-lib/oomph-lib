file_list=`ls animate_displacement*png`
for file in `echo $file_list`; do
    convert -resize 20%  $file `basename $file .png`.gif 
done
gifmerge -l0 animate_displacement?.gif animate_displacement??.gif > anim.gif
