#! /bin/sh

# clean up directories
rm -rf *~
cd ../
make clean
cd ../../../demo_drivers/navier_stokes/osc_plate/
make clean
cd ../../../
# remove excluded-files file
rm -rf excluded-files
# create new excluded-files file
find demo_drivers/navier_stokes/osc_plate/*/.svn \
     demo_drivers/navier_stokes/osc_plate/.svn \
     doc/navier_stokes/osc_plate/*/.svn \
     doc/navier_stokes/osc_plate/.svn \
     -maxdepth 12 -print > excluded-files

tar -cv -X excluded-files -f osc_plate.tar \
    doc/navier_stokes/osc_plate/ \
    demo_drivers/navier_stokes/osc_plate \
    bin/txt2h.sh \
    doc/oomph-lib_header.html \
    doc/figures/oomph_original.png \

gzip -9 osc_plate.tar

rm -rf excluded-files


    
