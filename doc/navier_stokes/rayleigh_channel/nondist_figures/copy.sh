#! /bin/sh

# clean up directories
rm -rf *~
cd ../
make clean
cd ../../../demo_drivers/navier_stokes/rayleigh_channel/
make clean
cd ../../../
# remove excluded-files file
rm -rf excluded-files
# create new excluded-files file
find demo_drivers/navier_stokes/rayleigh_channel/*/.svn \
     demo_drivers/navier_stokes/rayleigh_channel/.svn \
     doc/navier_stokes/rayleigh_channel/*/.svn \
     doc/navier_stokes/rayleigh_channel/.svn \
     -maxdepth 12 -print > excluded-files

tar -cv -X excluded-files -f rayleigh_channel.tar \
    doc/navier_stokes/rayleigh_channel/ \
    demo_drivers/navier_stokes/rayleigh_channel \
    bin/txt2h.sh \
    doc/oomph-lib_header.html \
    doc/figures/oomph_original.png \

gzip -9 rayleigh_channel.tar

rm -rf excluded-files


    
