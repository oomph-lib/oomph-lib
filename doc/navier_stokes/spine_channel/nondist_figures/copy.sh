#! /bin/sh

# clean up directories
cd ../
make clean
cd ../../../demo_drivers/navier_stokes/spine_channel_flow/
make clean
cd ../../../
# remove excluded-files file
rm -rf excluded-files
# create new excluded-files file
find demo_drivers/navier_stokes/spine_channel_flow/*/.svn \
     demo_drivers/navier_stokes/spine_channel_flow/.svn \
     doc/navier_stokes/spine_channel_flow/*/.svn \
     doc/navier_stokes/spine_channel_flow/.svn \
     -maxdepth 12 -print > excluded-files

tar -cv -X excluded-files -f spine_channel_flow.tar \
    doc/navier_stokes/spine_channel_flow/ \
    demo_drivers/navier_stokes/spine_channel_flow \
    bin/txt2h.sh \
    doc/oomph-lib_header.html \
    doc/figures/oomph_original.png \
    src/meshes/channel_spine_mesh*.h \
    src/meshes/channel_spine_mesh*.cc \


gzip -9 spine_channel_flow.tar

rm -rf excluded-files


    
