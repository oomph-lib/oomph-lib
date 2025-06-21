#! /bin/sh

# clean up directories
rm -rf *~
cd ../
make clean
cd ../../../demo_drivers/wave/two_d_wave/
make clean
cd ../../../
# remove excluded-files file
rm -rf excluded-files
# create new excluded-files file
find demo_drivers/wave/two_d_wave/*/.svn \
     demo_drivers/wave/two_d_wave/.svn \
     doc/wave/two_d_wave/*/.svn \
     doc/wave/two_d_wave/.svn \
     -maxdepth 12 -print > excluded-files

tar -cv -X excluded-files -f two_d_wave.tar \
    doc/wave/two_d_wave/ \
    demo_drivers/wave/two_d_wave \
    bin/txt2h.sh \
    doc/oomph-lib_header.html \
    doc/figures/oomph_original.png \

gzip -9 two_d_wave.tar

rm -rf excluded-files


    
