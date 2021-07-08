#! /bin/sh

# clean up directories
rm -rf *~
cd ../
make clean
cd ../../../demo_drivers/beam/unsteady_ring/
make clean
cd ../../../
# remove excluded-files file
rm -rf excluded-files
# create new excluded-files file
find demo_drivers/beam/unsteady_ring/*/.svn \
     demo_drivers/beam/unsteady_ring/.svn \
     doc/beam/unsteady_ring/*/.svn \
     doc/beam/unsteady_ring/.svn \
     -maxdepth 12 -print > excluded-files

tar -cv -X excluded-files -f unsteady_ring.tar \
    doc/beam/unsteady_ring/ \
    demo_drivers/beam/unsteady_ring \
    bin/txt2h.sh \
    doc/oomph-lib_header.html \
    doc/figures/oomph_original.png \

gzip -9 unsteady_ring.tar

rm -rf excluded-files


    
