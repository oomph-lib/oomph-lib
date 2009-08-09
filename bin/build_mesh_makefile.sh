#! /bin/sh

#-------------------------------------------------------------------
# Shell script to generate input files for mesh makefile from the
# list of filename stems in mesh_names.list
# The root directory must be passed as the first argument
#-------------------------------------------------------------------

# Test that we do have (only one) command line argument
if (test $# -ne 1); then
 echo "Path to working directory must be the only argument\
 to build_mesh_makefile.sh"
 exit 1
fi

#Set the mesh directory
mesh_dir=$1/src/meshes

#Check that the list file exists
if !(test -e $mesh_dir/mesh_names.list); then
 echo "File of mesh names (" $mesh_dir/mesh_names.list ") does not exist"
 echo "Have you passed the correct working directory\
 to build_mesh_makesfile.sh"
 exit 1
fi

# Get mesh filename stems from list file
#---------------------------------------
mesh_stem=`cat $mesh_dir/mesh_names.list`

#Write definitions of names
#--------------------------
echo "#Automatically generated file with filename definitions"\
 > $mesh_dir/mesh_names.aux
echo "#-------------------------------------------------------"\
 >> $mesh_dir/mesh_names.aux
echo " ">> $mesh_dir/mesh_names.aux

#old sources=`for i in $mesh_stem; do echo -n $i.template.cc " "; done`
sources=`for i in $mesh_stem; do printf "$i.template.cc "; done`
echo "#Define the templated sources" >> $mesh_dir/mesh_names.aux
echo "#----------------------------" >> $mesh_dir/mesh_names.aux
echo "sources=" $sources >> $mesh_dir/mesh_names.aux
echo " " >> $mesh_dir/mesh_names.aux

#old templated_headers=`for i in $mesh_stem; do echo -n $i.template.h " "; done`
templated_headers=`for i in $mesh_stem; do printf "$i.template.h "; done`
echo "#Define the templated headers" >> $mesh_dir/mesh_names.aux
echo "#----------------------------" >> $mesh_dir/mesh_names.aux
echo "templated_headers=" $templated_headers >> $mesh_dir/mesh_names.aux
echo " " >> $mesh_dir/mesh_names.aux

#old headers=`for i in $mesh_stem; do echo -n $i.h " "; done`
headers=`for i in $mesh_stem; do printf "$i.h "; done`
echo "#Define the combined header" >> $mesh_dir/mesh_names.aux
echo "#--------------------------" >> $mesh_dir/mesh_names.aux
echo "headers=" $headers >> $mesh_dir/mesh_names.aux
echo " " >> $mesh_dir/mesh_names.aux


# Write instructions how to generate the combined mesh
#-----------------------------------------------------
echo "# How to build the combined meshes " > $mesh_dir/mesh_instructions.aux
echo "#----------------------------------" >> $mesh_dir/mesh_instructions.aux
for i in $mesh_stem; do
  echo   $i.h": " $i.dummy  >> $mesh_dir/mesh_instructions.aux
  echo   $i.dummy": " $i.template.h " " $i.template.cc \
      >> $mesh_dir/mesh_instructions.aux
  echo  \
   "	 echo \"// Automatically generated, combined header file\" > $i.h" \
      >> $mesh_dir/mesh_instructions.aux
  echo  "	 echo \"#include \\\"$i.template.h\\\" \" >> $i.h " \
      >> $mesh_dir/mesh_instructions.aux
  echo "	 echo \"#include \\\"$i.template.cc\\\" \" >> $i.h " \
      >> $mesh_dir/mesh_instructions.aux
  echo " ">> $mesh_dir/mesh_instructions.aux
done

# Write instructions how to delete the combined mesh header files
#-----------------------------------------------------
echo "# How to delete the combined mesh header files "\
> $mesh_dir/mesh_clean.aux
echo "#----------------------------------" >> $mesh_dir/mesh_clean.aux
echo "clean-local:" >> $mesh_dir/mesh_clean.aux
#old echo -n "	 rm -f " >> $mesh_dir/mesh_clean.aux
printf "	 rm -f " >> $mesh_dir/mesh_clean.aux
for i in $mesh_stem; do
  #old echo -n $i.h" meshes.h " >> $mesh_dir/mesh_clean.aux
 printf "$i.h meshes.h " >> $mesh_dir/mesh_clean.aux
done
