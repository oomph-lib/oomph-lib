#stem="original_contact"
#stem="single_kink_contact"
#stem="kuhn_tucker_contact"
stem="sb_rays"

mkdir cropped
list=`ls $stem.*.png`
count=0
for file in `echo $list`; do
  convert -trim -transparent white $file cropped/$stem$count.png
  let count=count+1
done

