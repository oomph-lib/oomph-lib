#! /bin/sh

#--------------------------------------------------------------
# Shell script to remove licence (all lines starting 
# with //LIC//) from all oomph-lib *.h *.cc files.
#--------------------------------------------------------------


#--------------------------------------------------------------
# Find *.h and *.cc files in oomph-lib distribution (i.e. exclude
# external_src)
#--------------------------------------------------------------
oomph_lib_h_and_cc_files=`find demo_drivers src user_src user_drivers self_test bin \( -name '*.h' -o -name '*.cc' \) -exec ls  {} \; `


#--------------------------------------------------------------
# Loop over all of these files
#--------------------------------------------------------------
for file in $oomph_lib_h_and_cc_files; do
    echo "Delicencing:" $file
    gawk -f bin/remove_licence.awk $file > $file.delicenced
    mv -f $file.delicenced $file
done
