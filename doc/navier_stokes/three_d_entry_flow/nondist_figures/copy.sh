#! /bin/sh

# clean up directories
cd ../
make clean
cd ../../../demo_drivers/navier_stokes/three_d_entry_flow/
make clean
cd ../../../
# remove excluded-files file
rm -rf excluded-files
# create new excluded-files file
find demo_drivers/navier_stokes/three_d_entry_flow/*/.svn \
     demo_drivers/navier_stokes/three_d_entry_flow/.svn \
     doc/navier_stokes/three_d_entry_flow/*/.svn \
     doc/navier_stokes/three_d_entry_flow/.svn \
     -maxdepth 12 -print > excluded-files

tar -cv -X excluded-files -f three_d_entry_flow.tar \
    doc/navier_stokes/three_d_entry_flow/ \
    demo_drivers/navier_stokes/three_d_entry_flow \
    bin/txt2h.sh \
    doc/oomph-lib_header.html \
    doc/figures/oomph_original.png \
#    config/makefile_templates/doc \
#    configure.ac 

gzip -9 three_d_entry_flow.tar

rm -rf excluded-files


    
