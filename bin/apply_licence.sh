#! /bin/sh

#--------------------------------------------------------------
# Shell script to apply licence prefix to all oomph-lib
# *.h *.cc and  files. Also activate the svn keywords flag.
#--------------------------------------------------------------


#--------------------------------------------------------------
# Find *.h and *.cc files in oomph-lib distribution (i.e. exclude
# external_src)
#--------------------------------------------------------------
oomph_lib_h_and_cc_files=`find demo_drivers src user_src user_drivers \( -name '*.h' -o -name '*.cc' \) -exec ls  {} \; `


#--------------------------------------------------------------
# Loop over all of these files
#--------------------------------------------------------------
for file in $oomph_lib_h_and_cc_files; do
    echo "Licencing:" $file
    cat bin/cc_and_h_licence_block.txt $file > $file.licenced 
    gawk -f bin/remove_licence.awk $file.licenced > $file.delicenced
    if (test `diff $file $file.delicenced | wc -w` != 0); then
        echo "Original and licenced/delicenced file don't match!"
        echo "Please run bin/remove_licence.sh first."
        echo "Filename: " $file
        exit 1
    else
        echo "OK"
        mv -f $file.licenced $file
        rm  $file.delicenced
        svn propset svn:keywords "LastChangedDate LastChangedRevision" $file
    fi
#    echo " "
done
