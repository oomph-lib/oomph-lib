#! /bin/sh

# clean up directories
rm -rf *~
cd ../
make clean
cd ../../../demo_drivers/navier_stokes/osc_quarter_circle/
make clean
cd ../../../
# remove excluded-files file
rm -rf excluded-files
# create new excluded-files file
find demo_drivers/navier_stokes/osc_quarter_circle/*/.svn \
     demo_drivers/navier_stokes/osc_quarter_circle/.svn \
     doc/navier_stokes/osc_quarter_circle/*/.svn \
     doc/navier_stokes/osc_quarter_circle/.svn \
     -maxdepth 12 -print > excluded-files

tar -cv -X excluded-files -f osc_quarter_circle.tar \
    doc/navier_stokes/osc_quarter_circle/ \
    demo_drivers/navier_stokes/osc_quarter_circle \
    bin/txt2h.sh \
    doc/oomph-lib_header.html \
    doc/figures/oomph_original.png \

gzip -9 osc_quarter_circle.tar

rm -rf excluded-files


    
