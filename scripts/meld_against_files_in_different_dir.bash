#! /bin/bash

other_dir=/home/mheil/junk_tar_file_from_puneet/oomph-lib-1.0.1306M





base_dir=.

#h_and_cc_file_list=`find -name '*.h' -o -name '*.cc' -o -name '*.sh -o -name '*.bash' `

h_and_cc_file_list=`find -name '*.sh' `

echo $h_and_cc_file_list


for file in `echo $h_and_cc_file_list`; do
    orig_file=$other_dir"/"$base_dir"/"$file
    if [ $orig_file  ] && [ -e $file ]; then
        ls $orig_file $file
        if cmp -s "$orig_file" "$file"
        then
            echo "The files match"
        else
            echo "The files are different"
            chmod u-w $orig_file
            meld $file $orig_file
            chmod u+w $orig_file
        fi
    fi
done



