#! /bin/bash

if [ -e Compare_sample_point_containers ]; then
    echo "Directory Compare_sample_point_containers exists. Please delete."
    exit
fi
mkdir Compare_sample_point_containers

#make check
#mv Validation Compare_sample_point_containers/Validation_default

../../../../bin/wrapper_for_validate.sh --ref_bin
mv Validation Compare_sample_point_containers/Validation_ref_bin

../../../../bin/wrapper_for_validate.sh --non_ref_bin
mv Validation Compare_sample_point_containers/Validation_non_ref_bin

../../../../bin/wrapper_for_validate.sh --cgal
mv Validation Compare_sample_point_containers/Validation_cgal



