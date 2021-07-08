file_list=`ls comp*png`
for file in `echo $file_list`; do
    convert -resize 20%  $file ../figures/`basename $file .png`.gif 
    convert -resize 20%  $file ../figures/`basename $file .png`.eps
done
