#stem="original_contact"
#stem="single_kink_contact"
#stem="pseudo_melt_half"
#stem="proper_melt"
stem="melt_without_refreeze"

mkdir re_enumerated
list=`ls $stem.*.png`
count=0
for file in `echo $list`; do
  cp $file re_enumerated/$stem$count.png
  let count=count+1
done

