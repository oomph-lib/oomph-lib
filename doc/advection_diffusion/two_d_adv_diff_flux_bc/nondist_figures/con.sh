#! /bin/sh

rm -f *.gif

convert adv_diff_source1.png adv_diff_source1.gif
convert adv_diff_source2.png adv_diff_source2.gif

gifmerge -l0 -100 adv_diff_source?.gif > adv_diff_source.gif


cp adv_diff_source1.img.eps adv_diff_source.eps

rm -rf adv_diff_source1.img.eps adv_diff_source1.gif
rm -rf adv_diff_source2.img.eps adv_diff_source2.gif





