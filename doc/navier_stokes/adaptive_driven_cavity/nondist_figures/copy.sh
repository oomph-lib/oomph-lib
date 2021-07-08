#! /bin/sh

# clean up directories
cd ../
make clean
cd ../../../demo_drivers/navier_stokes/adaptive_driven_cavity/
make clean
cd ../../../
# remove excluded-files file
rm -rf excluded-files
# create new excluded-files file
find demo_drivers/navier_stokes/adaptive_driven_cavity/*/.svn \
     demo_drivers/navier_stokes/adaptive_driven_cavity/.svn \
     doc/navier_stokes/adaptive_driven_cavity/*/.svn \
     doc/navier_stokes/adaptive_driven_cavity/.svn \
     -maxdepth 12 -print > excluded-files

tar -cv -X excluded-files -f adaptive_driven_cavity.tar \
    doc/navier_stokes/adaptive_driven_cavity/ \
    demo_drivers/navier_stokes/adaptive_driven_cavity \
    bin/txt2h.sh \
    doc/oomph-lib_header.html \
    doc/figures/oomph_original.png \
#    config/makefile_templates/doc \
#    configure.ac 

gzip -9 adaptive_driven_cavity.tar

rm -rf excluded-files


    
