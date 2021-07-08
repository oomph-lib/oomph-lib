#! /bin/bash

clean_checkout_dir=/media/mheil/9e93acef-2d9b-499c-ae86-480fff25e3b1/scratch/version_clean_checkout_can_go



base_dir=src

changed_file_list=`svn status | grep 'M '  | awk '{print $2}'`
for changed_file in `echo $changed_file_list`; do
    orig_file=$clean_checkout_dir"/"$base_dir"/"$changed_file
    chmod u-w $orig_file
    meld $changed_file $orig_file
    chmod u+w $orig_file
done



